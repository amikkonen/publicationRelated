#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "RASModel.H"

namespace Foam
{
    typedef IncompressibleTurbulenceModel<transportModel> 
        transportModelIncompressibleTurbulenceModel;
    typedef RASModel<transportModelIncompressibleTurbulenceModel>
        RAStransportModelIncompressibleTurbulenceModel;
}

#include "kOmegaSSTSASnew.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    kOmegaSSTSASnew
);

#include "mykOmega.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    mykOmega
);

#include "mykOmegaSST.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    mykOmegaSST
);

//#include "myv2f.H"
//makeTemplatedTurbulenceModel
//(
//    transportModelIncompressibleTurbulenceModel,
//    RAS,
//    myv2f
//);
