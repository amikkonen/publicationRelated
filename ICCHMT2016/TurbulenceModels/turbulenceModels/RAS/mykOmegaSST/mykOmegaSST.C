/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "mykOmegaSST.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// my mod
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> mykOmegaSST<BasicTurbulenceModel>::mykOmegaSST::F1
(
    const volScalarField& CDkOmega
) const
{

    tmp<volScalarField> CDkOmegaPlus = CDkOmega*this->rho_;
    if (this->rho_.dimensions() == dimDensity)
    {
        CDkOmegaPlus = max
        (
            CDkOmegaPlus,
            dimensionedScalar("1.0e-10", dimless/sqr(dimTime)*dimDensity, 1.0e-10)
        );    

    } else {
        CDkOmegaPlus = max
        (
            CDkOmegaPlus, // value of rho_?
            dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
        );    
    }
    
//    tmp<volScalarField> CDkOmegaPlus = max
//    (
//        CDkOmega*this->rho_,
//        dimensionedScalar("1.0e-10", dimless/sqr(dimTime)*dimDensity, 1.0e-10)
//    ); 


    
    // added rho
    tmp<volScalarField> arg1 = min
    (
//        min
//        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_*this->rho_)*k_/(CDkOmegaPlus*sqr(y_))
        );//,
//        scalar(10)
//    );

    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykOmegaSST<BasicTurbulenceModel>::mykOmegaSST::F2() const
{
    tmp<volScalarField> arg2 = //min
    //(
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        );//,
//        scalar(100)
//    );

    return tanh(sqr(arg2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykOmegaSST<BasicTurbulenceModel>::mykOmegaSST::F3() const
{
//    tmp<volScalarField> arg3 = min
//    (
//        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
//        scalar(10)
//    );

    tmp<volScalarField> arg3 = 150*(this->mu()/this->rho_)/(omega_*sqr(y_));

    return 1 - tanh(pow4(arg3));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> mykOmegaSST<BasicTurbulenceModel>::mykOmegaSST::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23() *= F3();
    }

    return f23;
}


template<class BasicTurbulenceModel>
void mykOmegaSST<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void mykOmegaSST<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> mykOmegaSST<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> mykOmegaSST<BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> mykOmegaSST<BasicTurbulenceModel>::Qsas
(
    const volScalarField& S2,
    const volScalarField& gamma,
    const volScalarField& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
mykOmegaSST<BasicTurbulenceModel>::mykOmegaSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),
    // my mod
    SSTV_
    (
        Switch::lookupOrAddToDict
        (
            "SSTV",
            this->coeffDict_,
            false
        )
    ),
    // my mod end    
    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool mykOmegaSST<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());
        // my mod
        SSTV_.readIfPresent("SSTV", this->coeffDict());
        // my mod end        
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void mykOmegaSST<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
//    tmp<volScalarField> tOmega2 = 2*magSqr(skew(tgradU()));
    volScalarField Omega2(2*magSqr(skew(tgradU())));
    
  
//    volScalarField GbyNu((tgradU() && dev(twoSymm(tgradU()))));
//    volScalarField G(this->GName(), nut*GbyNu);

    volScalarField GbyNu((tgradU() && twoSymm(tgradU()))); // consistent SST 2003

//    volScalarField GbyNu();
    if ( SSTV_ ) {
    
//        Info<< "    KatoLaunder" << endl;
        // GbyNu = sqrt(S2) && sqrt(2*magSqr(skew(tgradU())));
 
        // GbyNu = tgradU() && 2*skew(tgradU());
//        GbyNu = tOmega2 - (2.0/3.0) * this->omega_ * divU; // SST-V
        GbyNu = Omega2 - (2.0/3.0) * this->omega_ * divU; // SST-V        
    }

// muuttuuko k omegaksi kerronnassa
//    GbyNu = sqrt(S2) && sqrt(Omega2); // kato launder
    volScalarField G(this->GName(), nut*GbyNu);
    
    
//    volScalarField P((nut*rho * tgradU() && twoSymm(tgradU())));
//    volScalarField P("P", nut*rho*PbyMut);



    // RC-Hellsten
//    scalar Crc = 1.4;
//    volScalarField Ri  = sqrt(Omega2)/sqrt(S2) 
//                         * ( sqrt(Omega2)/sqrt(S2) - scalar(1) );
//    volScalarField F4 = (1 / (1 + Crc*Ri));





    // RC-
//    scalar Cr1 = 1.0;
//    scalar Cr2 = 2.0;
//    scalar Cr3 = 1.0;
    
//    volScalarField D2 = max(S2,0.09*sqr(omega_) );    
//    volTensorField S =  symm( tgradU() ) ;

    // TODO
//    volTensorField S = scalar(0.5) * (tgradU() + tgradU().T());
//    volTensorField W = scalar(0.5) * (tgradU() - tgradU().T());

//    volTensorField W(skew(tgradU()));

    // note that reference frame rotation term is missing (not implemented)
    
    
//    volVectorField gradS = fvc::laplacian(U);
//    
//    volScalarField rhat = 
//            scalar(2) * W && S / (sqrt(Omega2) && pow(D2,3/2)) & gradS;
            
            //* (U & fvc::grad(S));
            
//    volScalarField rStar = sqrt(S2)/sqrt(Omega2);
    
    
    // half of the terms missing
//    volScalarField frot = (1+Cr1) * 2 * rStar / (1+rStar) * 0.5 - Cr1;
    
//    frot = fvc::smooth(frot, 0.5);
    
    
    
    
//    volScalarField fr1 = max(
//                                min(frot, 1.25), 0.1
//                             );
    
    
    // write
//    volScalarField rho
//    (
//        IOobject
//        (
//            "rho",
//            runTime.timeName(),
//            mesh
//        ),
//        thermo.rho()
//    );
    
//    volScalarField fr1_write
//    (
//        IOobject
//        (
//            IOobject::groupName("fr1", U.group()),
//            this->runTime_.timeName(),
//            this->mesh_,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        fr1
//    );
//    fr1_write.write();
    


    
    
    
//    GbyNu *= fr1;
//    G *= fr1;
    
    
    tgradU.clear();
//    tOmega2.clear();

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));

    {
        volScalarField gamma(this->gamma(F1));
        volScalarField beta(this->beta(F1));

        // TODO
        dimensionedScalar smallNut(
                dimensionedScalar("smallNut", nut.dimensions(), ROOTVSMALL));



        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
//            fr1*alpha*rho*gamma 
            alpha*rho*gamma 
           *
//           GbyNu
           min
            (
                GbyNu, 
                //P,
                (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, b1_*F23()*sqrt(S2))
                //c1_*betaStar_*omega_*omega_
                //c1_*betaStar_*omega_*k_ / (nut + smallNut)
            )
          //- fvm::SuSp((2.0/3.0)*alpha*rho*gamma*divU, omega_)
          - fvm::Sp(alpha*rho*beta*omega_, omega_)
//          - fvm::Sp(F4*alpha*rho*beta*omega_, omega_)
          - fvm::SuSp
            (
                alpha*rho*(F1 - scalar(1))*CDkOmega/omega_,
                omega_
            )
          + Qsas(S2, gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn().relax();
        fvOptions.constrain(omegaEqn());
        omegaEqn().boundaryManipulate(omega_.boundaryField());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
//        fr1*min(alpha*rho*G, (c1_*betaStar_)*alpha*rho*k_*omega_)
        min(alpha*rho*G, (c1_*betaStar_)*alpha*rho*k_*omega_)
//      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
//       min(alpha*P, (c1_*betaStar_)*alpha*rho*k_*omega_) 
      - fvm::Sp(alpha*rho*betaStar_*omega_, k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn().relax();
    fvOptions.constrain(kEqn());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(S2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
