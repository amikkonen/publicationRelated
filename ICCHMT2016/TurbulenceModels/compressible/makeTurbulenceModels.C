#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

namespace Foam
{

typedef ThermalDiffusivity<CompressibleTurbulenceModel<fluidThermo> >
    fluidThermoCompressibleTurbulenceModel;

typedef RASModel<EddyDiffusivity<fluidThermoCompressibleTurbulenceModel> >
    RASfluidThermoCompressibleTurbulenceModel;

}

#include "kOmegaSSTSASnew.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    kOmegaSSTSASnew
);

#include "mykOmega.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    mykOmega
);

#include "mykOmegaSST.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    mykOmegaSST
);


//#include "myv2f.H"
//makeTemplatedTurbulenceModel
//(
//    fluidThermoCompressibleTurbulenceModel,
//    RAS,
//    myv2f
//);
