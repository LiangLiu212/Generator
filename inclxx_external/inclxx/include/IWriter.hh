#ifndef G4INCL_IO_IWriter_hh
#define G4INCL_IO_IWriter_hh 1

#include <string>

#include "G4INCLGlobalInfo.hh"
#include "G4INCLEventInfo.hh"

namespace G4INCL {
  namespace IO {

    class IWriter {
    public:
      IWriter();
      virtual ~IWriter();

      virtual bool openFile(std::string filename) = 0;
      virtual bool writeGlobalInfo(G4INCL::GlobalInfo const &globalInfo) = 0;
      virtual bool writeEventInfo(G4INCL::EventInfo const &eventInfo) = 0;
      virtual bool closeFile() = 0;
      virtual void flush() = 0;
    };
  }
}

#endif
