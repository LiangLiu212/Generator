THIS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export GENIEBASE=/root/genie
export GENIE=$GENIEBASE/Generator
export XSECSPLINEDIR=$GENIEBASE/data
export GENIE_REWEIGHT=$GENIEBASE/Reweight
export PATH=$GENIE/bin:$GENIE_REWEIGHT/bin:$PATH
export PYTHIA=/opt/pythia6/v6_428/lib
export LHAPDF_INC=`lhapdf-config --incdir`
export LHAPDF_LIB=`lhapdf-config --libdir`
export LOG4CPP_INC=`log4cpp-config --prefix`/include
export LOG4CPP_LIB=`log4cpp-config --libdir`
export LD_LIBRARY_PATH=$GENIE/lib:$GENIE_REWEIGHT/lib:$LD_LIBRARY_PATH
unset GENIEBASE

./configure \
  --enable-gsl \
  --enable-rwght \
  --enable-lhapdf6 \
  --enable-boost \
  --enable-incl \
  --with-log4cpp-inc=${LOG4CPP_INC} \
  --with-log4cpp-lib=${LOG4CPP_LIB} \
  --with-libxml2-inc=${LIBXML2_INC} \
  --with-libxml2-lib=${LIBXML2_LIB} \
  --with-lhapdf6-lib=${LHAPDF_LIB} \
  --with-lhapdf6-inc=${LHAPDF_INC} \
  --with-pythia6-lib=${PYTHIA} \
  --with-boost-lib=/usr/local/lib \
  --with-boost-inc=/usr/local/include \
  --with-incl-lib=/root/inclxx/install/lib \
  --with-incl-inc=/root/inclxx/install/include 
#  --enable-geant4 \
