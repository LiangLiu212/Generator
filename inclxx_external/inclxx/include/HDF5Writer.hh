#ifndef G4INCL_IO_HDF5Writer_hh
#define G4INCL_IO_HDF5Writer_hh 1

#ifdef INCL_HDF5_USE

#include <iostream>
#include <fstream>
#include "IWriter.hh"

#include "H5LT.h"
#include "hdf5.h"

namespace G4INCL {
  namespace IO {

    class HDF5Writer : public IWriter {
    public:
      HDF5Writer();
      ~HDF5Writer();

      bool openFile(std::string filename);
      bool writeParticle(G4INCL::Particle *p);
      bool closeFile();
      bool flush() {}

    private:
      hid_t fileID, space, dset;
      bool outputStreamOpen;
      hsize_t dims[2];
      hsize_t start[2];
      hsize_t stride[2];
      hsize_t count[2];
      hsize_t block[2];

    };
  }
}

#endif
#endif
