#ifndef CONVERTER_HH
#define CONVERTER_HH

#include "TLeaf.h"
#include "ROOTWriter.hh"
#ifdef HAS_THBOOK
#include "THbookFile.h"
#include "THbookTree.h"
#endif
#include <map>
#include <vector>

class Converter {
  public:
    Converter(std::vector<std::string> const &inFilename, std::string const &outFilename);
    ~Converter();

    int doIt();

  private:
    bool parseTextFile(std::string const &textFilename);
    void initialiseVariableMap();
    void initialiseLeafMaps();
    void initialiseLeafMapsFromVarMap(std::map<std::string,std::string> const &theVarMap);
    int fillGlobalTree(Long64_t nEntries);
    int projectileTypeToAp(const int t);
    int projectileTypeToZp(const int t);

    int exitStatus;
    Int_t nShots;
    Int_t totalNShots;
    Float_t geometricalCrossSection;
    Float_t largestGeometricalCrossSection;
    Short_t Ap;
    Short_t Zp;
    Float_t Ep;
    Short_t At;
    Short_t Zt;

    TFile *theInFile;
#ifdef HAS_THBOOK
    THbookFile *theInHbookFile;
#endif
    TTree *theInTree;
    TTree *theOutTree;
    TTree *theOutGlobalTree;
    std::vector<std::string> const &theInFilenames;
    G4INCL::IO::ROOTWriter theROOTWriter;
    std::map<std::string, void *> theInLeafMap;
    std::map<std::string, void *> theOutLeafMap;
    std::map<std::string, std::string> theShortMap;
    std::map<std::string, std::string> theShortArrayMap;
    std::map<std::string, std::string> theIntMap;
    std::map<std::string, std::string> theFloatMap;
    std::map<std::string, std::string> theFloatArrayMap;

    typedef std::map<std::string, std::string>::const_iterator VariableIterator;
    typedef std::map<std::string, void *> LeafIterator;
};

#endif // CONVERTER_HH
