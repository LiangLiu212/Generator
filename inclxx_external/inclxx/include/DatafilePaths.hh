#ifndef DATAFILEPATHS_HH
#define DATAFILEPATHS_HH

#include <string>

#ifdef USE_INSTALL_PATH
const std::string dataPath = "/exp/sbnd/app/users/liangliu/sbnd_genie/inclxx_external/work/stage/inclxx/share";
#else
const std::string dataPath = "/exp/sbnd/app/users/liangliu/sbnd_genie/inclxx_external/work/src";
#endif

const std::string defaultINCLXXDatafilePath = dataPath + "/data/";
#ifdef INCL_DEEXCITATION_ABLA07
const std::string defaultABLA07DatafilePath = dataPath + "/de-excitation/abla07/upstream/tables/";
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
const std::string defaultABLAXXDatafilePath = dataPath + "/data/";
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
const std::string defaultGEMINIXXDatafilePath = dataPath + "/de-excitation/geminixx/upstream/";
#endif

#endif // DATAFILEPATHS_HH
