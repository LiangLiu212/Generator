#ifndef G4INCL_IO_ASCIIWriter_hh
#define G4INCL_IO_ASCIIWriter_hh 1

#include <iostream>
#include <fstream>
#include "IWriter.hh"

namespace G4INCL {
  namespace IO {

    class ASCIIWriter : public IWriter {
    public:
      ASCIIWriter(const bool invKin);
      ~ASCIIWriter();

      bool openFile(std::string filename);
      /** \brief Writes event information in ASCII format
       *
       * Uses a simple markup language.
       * An "e" character denotes the start of a new event. On the next line,
       * generic information about the event such as projectile and target
       * characteristics, impact parameter, etc.
       * A "p" character denotes an emitted particle. Information about the
       * particle follows.
       * An "r" character denotes a remnant. Information about the remnant
       * follows.
       *
       * \param eventInfo event information to be written
       * \return true
       */
      bool writeEventInfo(G4INCL::EventInfo const &eventInfo);

      /** \brief Writes global information in ASCII format
       *
       * Uses a simple markup language.
       * A "g" character denotes the start of the global calculation information.
       * Information follows on the same line.
       *
       * \param globalInfo global information to be written
       * \return true
       */
      bool writeGlobalInfo(G4INCL::GlobalInfo const &globalInfo);
      bool closeFile();
      void flush();

    private:
      std::ofstream *outputStream;
      bool outputStreamOpen;
      const bool inverseKinematics;
    };
  }
}

#endif
