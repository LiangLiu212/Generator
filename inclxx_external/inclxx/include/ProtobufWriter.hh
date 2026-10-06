#ifndef G4INCL_IO_ProtobufWriter_hh
#define G4INCL_IO_ProtobufWriter_hh 1

#ifdef INCL_PROTOBUF_USE

#include <iostream>
#include <fstream>
#include "IWriter.hh"
#include <fcntl.h>
#include <sys/stat.h>

#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include "eventio.pb.h"

namespace G4INCL {
  namespace IO {

    class ProtobufWriter : public IWriter {
    public:
      ProtobufWriter();
      ~ProtobufWriter();

      bool openFile(std::string filename);
      /** \brief Writes event information in Protobuf format
       *
       * \param eventInfo event information to be written
       * \return true
       */
      bool writeEventInfo(G4INCL::EventInfo const &eventInfo);

      /** \brief Writes global information in Protobuf format
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
      int fd;
      google::protobuf::io::FileOutputStream *file_stream;
      google::protobuf::io::GzipOutputStream *gzip_stream;
      std::ofstream *outputStream;
      bool outputStreamOpen;
    };
  }
}

#endif
#endif
