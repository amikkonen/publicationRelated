/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "wallHeatTransfer.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "fvCFD.H"

#include "wallFvPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallHeatTransfer, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallHeatTransfer::wallHeatTransfer
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    phiName_("phi"),
    heName_("h"),    
    TRef_(dimensionedScalar("zero", dimTemperature, 0.0)),
    obr_(obr),
    active_(true),
    log_(true),
//    writeHeatFlux_(true),
//    writeHeatTransferCoefficient_(true),
    patchSet_(),
    heatFlux_
    (
        IOobject
        (
            "heatFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("heatFlux", heatFlux_.dimensions(), 0.0)
    ),
    heatTransferCoefficient_
    (
        IOobject
        (
            "heatTransferCoefficient",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("hHeat",dimPower / sqr(dimLength) / dimTemperature, 0)
    ),
    meanHeatTransferCoefficient_(SMALL),
    convergenceCriteria_(0)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "heatTransfer::wallHeatTransfer"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallHeatTransfer::~wallHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallHeatTransfer::read(const dictionary& dict)
{
    if (active_)
    {

        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
        heName_ = dict.lookupOrDefault<word>("heName", "h");
        log_ = dict.lookupOrDefault<Switch>("log", true);
        TRef_ = dict.lookup("Tref");
        convergenceCriteria_ = dict.lookupOrDefault<scalar>("converged", SMALL);
        
        
//        writeHeatFlux_ = dict.lookupOrDefault<Switch>("writeHeatFlux", true);
//        writeHeatTransferCoefficient_ = 
//            dict.lookupOrDefault<Switch>("writeHeatTransferCoefficient", true);


        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        patchSet_ = pbm.patchSet(wordReList(dict.lookup("patches")));
    }
}


void Foam::wallHeatTransfer::writeFileHeader(const label i)
{
//    OFstream& file = this->file();

//    writeHeader(file, "Wall heat transfer");
//    writeCommented(file, "Time");

//    writeTabbed(file, "wallHeatFlux");

//    file<< endl;
}


void Foam::wallHeatTransfer::execute()
{
    // Do nothing - only valid on write
}


void Foam::wallHeatTransfer::end()
{
    // Do nothing - only valid on write
}


void Foam::wallHeatTransfer::timeSet()
{
    // Do nothing - only valid on write
}

void Foam::wallHeatTransfer::compressibleHeatTransfer()
{
    // Turbulence model.
    typedef compressible::turbulenceModel cmpModel;
    const cmpModel& turbulence = mesh_.lookupObject<cmpModel>
    (
        turbulenceModel::propertiesName
    );

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
//        const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
//        volScalarField Cp = thermo.Cp();

//        //    const volScalarField& mu = thermo.mu();
//        //    const volScalarField& mut = turbulence.mut();
//        const volScalarField Pr = mesh_.lookupObject<volScalarField>
//        (
//            fluidThermo::Pr
//        );

//        //    volScalarField k = turbulence.alphaEff()() * rho * Cp;


//        const volScalarField& T =
//            mesh_.lookupObject<volScalarField>("T");

//        volScalarField  thermalConductivity = 
//            (turbulence.nu()/Pr+turbulence.nut()/Prt)*Cp*rho;
//        
//        heatFlux_.boundaryField()[patchI] =
//        		 thermalConductivity.boundaryField()[patchI] 
//        		 * T.boundaryField()[patchI].snGrad();




        // Enthalpy / internal energy field.
        const volScalarField& he =
                mesh_.lookupObject<volScalarField>(heName_);

        heatFlux_.boundaryField()[patchI] =
	         turbulence.alphaEff()().boundaryField()[patchI] 
	         * he.boundaryField()[patchI].snGrad();
    }
}

void Foam::wallHeatTransfer::incompressibleHeatTransfer()
{
    const volScalarField& T =
            mesh_.lookupObject<volScalarField>("T"); 

    const dictionary& transportProperties = mesh_.lookupObject<IOdictionary>
    (
       "transportProperties"
    );
    
    // Prandtl number.
    dimensionedScalar Pr(transportProperties.lookup("Pr"));

    // Turbulent Prandtl number.
    dimensionedScalar Prt(transportProperties.lookup("Prt"));

    // Heat capacity.
    dimensionedScalar Cp(transportProperties.lookup("Cp"));

    // Fluid density.
    dimensionedScalar rho(transportProperties.lookup("rho"));
    
    // Turbulence model.
    typedef incompressible::turbulenceModel icoModel;
    const icoModel& turbulence = mesh_.lookupObject<icoModel>
    (
        turbulenceModel::propertiesName
    );

    // Thermal conductivity.
    volScalarField  thermalConductivity = 
        (turbulence.nu()/Pr+turbulence.nut()/Prt)*Cp*rho;

    // Loop all given patches and calculate wall heat flux.
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
	    heatFlux_.boundaryField()[patchI] =
        		 thermalConductivity.boundaryField()[patchI] 
        		 * T.boundaryField()[patchI].snGrad();    
    }
}

void Foam::wallHeatTransfer::calcHeatTransferCoefficient()
{    
    // Temperature field.    
    const volScalarField& T =
            mesh_.lookupObject<volScalarField>("T");   

    // Loop all given patches and calculate wall heat transfer coefficient.
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
	    heatTransferCoefficient_.boundaryField()[patchI] =
            heatFlux_.boundaryField()[patchI]			    
            / (T.boundaryField()[patchI] - scalar(TRef_.value()));
    }
    
    // Write results.
    mesh_.time().write();
}

void Foam::wallHeatTransfer::writeLog()
{

    // Area
    const surfaceScalarField::GeometricBoundaryField& magSf =
        mesh_.magSf().boundaryField();
       
    // Loop all given patches and calculate mean heat flux
    // and mean heat transfer coefficient for log.    
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
        scalar totalPatchHeatFlux = 
            gSum(magSf[patchI]*heatFlux_.boundaryField()[patchI]);

        scalar old_meanHeatTransferCoefficient = meanHeatTransferCoefficient_;
            
        meanHeatTransferCoefficient_ = 
            gSum(magSf[patchI]
            * heatTransferCoefficient_.boundaryField()[patchI])
            / gSum(magSf[patchI]);

        scalar diff = meanHeatTransferCoefficient_ -            
                        old_meanHeatTransferCoefficient;

        scalar maxHeatTransferCoefficient = 
            max(heatTransferCoefficient_.boundaryField()[patchI]);


        Info<< mesh_.boundary()[patchI].name() << endl
             << "    Total heat power               (W)      " 
             << totalPatchHeatFlux << endl
             << "    h - h_old                      (W/m^2K) "              
             << diff << endl
             << "    Max heat transfer coefficient  (W/m^2K) " 
             << maxHeatTransferCoefficient << endl
             << "    Mean heat transfer coefficient (W/m^2K) " 
             << meanHeatTransferCoefficient_  << endl;

        
//        if (mag(diff) < convergenceCriteria_) {
//            obr_.time().stopAt(Time::saWriteNow);
//            Info << endl
//                 << "Mean heat transfer coefficient converged." << endl
//                 << "   h - h_old < converged" << endl;              
//        }
  

             
    }
    Info << endl;
}

void Foam::wallHeatTransfer::write()
{

    if (active_)
    {
//        functionObjectFile::write();
        
        // Needed for choosing between compressible and incompressible case.
        const surfaceScalarField& phi =
                mesh_.lookupObject<surfaceScalarField>(phiName_);

        // Call correct model for heat flux calculation
        // (compressible/incompressible).
        if (phi.dimensions() == dimMass/dimTime)
        {
            compressibleHeatTransfer();
        }
        else if (phi.dimensions() == dimVolume/dimTime)
        {
            incompressibleHeatTransfer();
        }
        else
        {
            NotImplemented;
        }
        
        // Calculate heat transfer coefficient.        
        calcHeatTransferCoefficient();
        
        // Write log if log switch is on.
        if (log_)
        {
            writeLog();
        }
    }
}

// ************************************************************************* //
