#include "IWriter.hh"

#ifndef G4INCL_IO_INCLTree_hh
#define G4INCL_IO_INCLTree_hh 1

namespace G4INCL {
  namespace IO {

    class INCLTree {
    public:
      INCLTree(IWriter *aWriter, std::string file) {
        theWriter = aWriter;
        theWriter->openFile(file);
      };
      ~INCLTree() {
        theWriter->closeFile();
        delete theWriter;
      };

      bool fillEventInfo(G4INCL::EventInfo const &eventInfo) {
        sanityChecks(eventInfo);
        theWriter->writeEventInfo(eventInfo);
        return true;
      }

      bool fillGlobalInfo(G4INCL::GlobalInfo const &globalInfo) {
        theWriter->writeGlobalInfo(globalInfo);
        return true;
      }

      void flush() {
        theWriter->flush();
      }

      void sanityChecks(G4INCL::EventInfo const &eventInfo) const;

    private:
      IWriter *theWriter;
    };
  }
}
#endif
