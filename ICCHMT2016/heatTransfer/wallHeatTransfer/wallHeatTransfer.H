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

Class
    Foam::wallHeatTransfer

Group
    grpHeatTransferFunctionObjects

Description
    Short description...
    

    Example of function object specification:
    \verbatim
    wallHeatTransfer1
    {
        type        wallHeatTransfer;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
    \endtable

    Output data is written to the files \<timeDir\>/wallHeatFlux
                                        \<timeDir\>/heatTransferCoefficient

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    wallHeatTransfer.C
    IOwallHeatTransfer.H

\*---------------------------------------------------------------------------*/

#ifndef wallHeatTransfer_H
#define wallHeatTransfer_H


#include "volFields.H"
#include "functionObjectFile.H"
#include "Switch.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class wallHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class wallHeatTransfer
:
    public functionObjectFile
{
protected:

    // Protected data

        //- Name of this set of field min/max
        //  Also used as the name of the output directory
        word name_;
        
        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Name of flux field (optional)
        word phiName_;
        
        //- Name of entalpy/internal energy field (optional)
        word heName_;

        // Reference value for heat transfer coefficient (h = q / (Twall - Tref))
        dimensionedScalar TRef_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well
        Switch log_;

//        //- Write wallHeatFlux
//        Switch writeHeatFlux_;

//        //- Write heatFlux
//        Switch writeHeatTransferCoefficient_;
//            
        //- Patches to to calculated heat transfer
        labelHashSet patchSet_;
        

        //- Heat flux
        volScalarField heatFlux_;
        
        //- Heat transfer coeffiecint
        volScalarField heatTransferCoefficient_;


        scalar meanHeatTransferCoefficient_;
        scalar convergenceCriteria_;
    
    

    // Private Member Functions
        //- Disallow default bitwise copy construct
        wallHeatTransfer(const wallHeatTransfer&);

        //- Disallow default bitwise assignment
        void operator=(const wallHeatTransfer&);

        //- Output file header information
        virtual void writeFileHeader(const label i);

        void compressibleHeatTransfer();
        void incompressibleHeatTransfer();

        void calcHeatTransferCoefficient();
        void writeLog();


public:

    //- Runtime type information
    TypeName("wallHeatTransfer");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        wallHeatTransfer
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
            //,const bool readFields = true
        );


    //- Destructor
    virtual ~wallHeatTransfer();


    // Member Functions

        //- Return name of the set of field min/max
        virtual const word& name() const
        {
            return name_;
        }

        //- Read options, default from controlDict
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the wallHeatTransfer
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
