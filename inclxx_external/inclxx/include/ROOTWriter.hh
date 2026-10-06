#ifndef G4INCL_IO_ROOTWriter_hh
#define G4INCL_IO_ROOTWriter_hh 1

#ifdef INCL_ROOT_USE

//#include "TRint.h"
#include <TFile.h>
#include <TTree.h>
#include <string>
#include "G4INCLEventInfo.hh"
#include "G4INCLGlobalInfo.hh"
#include "IWriter.hh"

namespace G4INCL {
  namespace IO {
    class ROOTWriter : public IWriter {
    public:
      ROOTWriter(std::string const &aSelection="", const bool c=false, const bool invKin=false);
      ~ROOTWriter();

      bool openFile(std::string filename);
      bool writeEventInfo(G4INCL::EventInfo const &e);
      bool writeGlobalInfo(G4INCL::GlobalInfo const &e);
      bool closeFile();

      void flush();

      void selectEvents();

      TTree *getEventTree() const { return theEventTree; }
      TTree *getGlobalTree() const { return theGlobalTree; }

      void setCompressionLevel(const int l) { compressionLevel=l; }
      int getCompressionLevel() const { return compressionLevel; }
    private:
      TFile *outputStream;
      TTree *theGlobalTree, *theEventTree, *theOutputEventTree;
      std::string selection;
      bool outputStreamOpen;
      G4INCL::GlobalInfo theGlobalInfo;
      G4INCL::EventInfo theEventInfo;
      int compressionLevel;
      const bool concise;
      const bool inverseKinematics;
    };
  }
}

#endif
#endif
