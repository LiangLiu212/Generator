export LHAPDF_INC=`lhapdf-config --incdir`
export LHAPDF_LIB=`lhapdf-config --libdir`
export LOG4CPP_INC=`log4cpp-config --prefix`/include
export LOG4CPP_LIB=`log4cpp-config --libdir`

./configure \
  --enable-gsl \
  --enable-rwght \
  --enable-lhapdf6 \
  --disable-lhapdf5 \
  --disable-pythia8 \
  --enable-pythia6 \
  --enable-boost \
  --enable-incl \
  --with-log4cpp-inc=${LOG4CPP_INC} \
  --with-log4cpp-lib=${LOG4CPP_LIB} \
  --with-libxml2-inc=${LIBXML2}/include/libxml2/ \
  --with-libxml2-lib=${LIBXML2}/lib \
  --with-lhapdf6-lib=${LHAPDF_LIB} \
  --with-lhapdf6-inc=${LHAPDF_INC} \
  --with-incl-lib=${INCLXX_DIR}/lib \
  --with-incl-inc=${INCLXX_DIR}/include \
  --with-pythia6-lib=${PYTHIA6}/lib \
  --with-boost-lib=${BOOST}/lib \
  --with-boost-inc=${BOOST}/include
