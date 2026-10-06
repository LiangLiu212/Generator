#!/bin/bash
########################################################################
# build_inclxx_external.sh
#
# (Re)build the INCL++ v6.34 binary external kept in this directory:
#
#   inclxx/bin      inclxx-config (used by GENIE's Make.include), thisinclxx.sh
#   inclxx/lib      the 9 shared libraries
#   inclxx/include  headers
#   inclxx/share    run-time data tables only (INCL, ABLA07, ABLAXX, GEMINI++)
#
# It is NOT a UPS product and contains NO INCL++ / de-excitation source code.
# GENIE is built against it (--enable-incl) and it is installed inside the
# genie product, e.g. genie v3_06_02_sbn4 -q e26:incl634:prof -- see README.md.
#
# You only need this script to change the INCL++ version or its build options.
# It needs read access to the INCL++ source (GENIE-interface fork,
# github.com/LiangLiu212/inclxx, not public), builds it with the e26:prof
# toolchain (gcc 12.1.0, root v6_28_12, boost v1_82_0) and REPLACES ./inclxx.
# Commit the result afterwards.
#
# Usage:  ./build_inclxx_external.sh          (enters the SL7 container itself)
# Env:    INCLXX_REPO  git repo to build from   (default: Liang's clone on the SBND gpvms)
#         INCLXX_REF   commit / branch / tag    (default: 8a6c1cf)
#         WORK         scratch build directory  (default: ./work, git-ignored)
########################################################################
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INCLXX_REPO=${INCLXX_REPO:-/exp/sbnd/app/users/liangliu/GENIE/inclpp/inclxx}
INCLXX_REF=${INCLXX_REF:-8a6c1cf}
INCLXX_VERSION=6.34
CQUAL=e26
BTYPE=prof
WORK=${WORK:-${HERE}/work}

die() { echo "FATAL: $*" 1>&2; exit 1; }

# The e26 toolchain only works on SL7: re-run inside the container.
if ! grep -q "release 7" /etc/redhat-release 2>/dev/null; then
  IMG=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:jsl
  APPT=/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer
  [ -x "$APPT" ] || APPT=$(command -v apptainer)
  exec "$APPT" exec --pid --ipc -B /etc/hosts,/tmp,/opt,/cvmfs,/exp/,/pnfs/,/nashome \
    "$IMG" /bin/bash "${HERE}/$(basename "${BASH_SOURCE[0]}")" "$@"
fi

echo "===== $(date) : INCL++ ${INCLXX_VERSION} external build START ($(cat /etc/redhat-release)) ====="

# Same toolchain as genie -q e26:prof (gcc 12.1.0, root v6_28_12, boost v1_82_0).
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup root v6_28_12 -q ${CQUAL}:p3915:${BTYPE} || die "setup root"
setup boost v1_82_0 -q ${CQUAL}:${BTYPE}       || die "setup boost"
setup cmake v3_27_4                            || die "setup cmake"
echo "gcc    : $(which gcc)  ($(gcc -dumpfullversion))"
echo "root   : ${ROOTSYS}"
echo "boost  : ${BOOST_FQ_DIR}"
echo "cmake  : $(which cmake)"

SRC=${WORK}/src; BLD=${WORK}/build; STAGE=${WORK}/stage; PREFIX=${STAGE}/inclxx
rm -rf "${WORK}"; mkdir -p "${BLD}" "${STAGE}" || die "cannot create ${WORK}"

# Pristine checkout of the pinned commit (the working clone is left untouched).
git clone -q --no-hardlinks "${INCLXX_REPO}" "${SRC}" || die "git clone ${INCLXX_REPO}"
git -C "${SRC}" checkout -q "${INCLXX_REF}"           || die "git checkout ${INCLXX_REF}"
COMMIT=$(git -C "${SRC}" rev-parse HEAD)
echo "source : ${INCLXX_REPO} @ ${COMMIT}"

# Same configuration as the validated install in GENIE/inclpp/install.
CMAKE_OPTS=(
  -DCMAKE_INSTALL_PREFIX=${PREFIX}
  -DCMAKE_BUILD_TYPE=RelWithDebInfo
  -DBUILD_SHARED_LIBRARY=ON
  -DUSE_INSTALL_PATH=ON
  -DINCL_ROOT_USE=ON
  -DINCL_DEBUG_LOG=ON
  -DINCL_SIGNAL_HANDLING=ON
  -DINCL_USE_ALLOCATION_POOL=ON
  -DINCL_PEDANTIC_COMPILER=ON
  -DINCL_DEFINE_EMPTY_G4THREADLOCAL=ON
  -DINCL_AVATAR_SEARCH=MinElement
  -DINCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON=Set
  -DINCL_DEEXCITATION_USE=ON
  -DINCL_DEEXCITATION_ABLA07=ON
  -DINCL_DEEXCITATION_ABLA07_USER_RNG=ON
  -DINCL_DEEXCITATION_ABLAXX=ON
  -DINCL_DEEXCITATION_ABLAXX_USER_RNG=ON
  -DINCL_DEEXCITATION_FERMI_BREAKUP=ON
  -DINCL_DEEXCITATION_GEMINIXX=ON
  -DINCL_DEEXCITATION_GEMINIXX_USER_RNG=ON
  -DINCL_DEEXCITATION_SMM=ON
  -DINCL_DEEXCITATION_SMM_USER_RNG=ON
  -DINCL_ASCII_USE=OFF
  -DINCL_HDF5_USE=OFF
  -DINCL_PROTOBUF_USE=OFF
  -DINCL_COUNT_RND_CALLS=OFF
  -DINCL_REGENERATE_AVATARS=OFF
  -DINCL_TREEADD=OFF
  -DINCL_TREECONVERTER=OFF
)
cd "${BLD}" || die "cd ${BLD}"
cmake "${CMAKE_OPTS[@]}" "${SRC}" || die "cmake configure failed"
make -j"$(nproc)"                 || die "make failed"
make install                      || die "make install failed"

echo "===== $(date) : pruning the install to binary + headers + data ====="
cd "${PREFIX}" || die "cd ${PREFIX}"

# GENIE only needs the libraries and inclxx-config; drop the standalone app.
rm -f bin/INCLCascade

# G4INCLVersion.hh is generated during the build, i.e. after cmake globbed the
# headers to install, so a build from a clean checkout does not install it.
cp -p "${SRC}/utils/include/G4INCLVersion.hh" include/ || die "G4INCLVersion.hh not generated"

# 'make install' copies the whole GEMINI++ upstream tree (sources, docs) to
# share/. At run time GEMINI++ reads only ${GINPUT}/tbl and ${GINPUT}/tl.
G=share/de-excitation/geminixx/upstream
[ -d ${G}/tbl ] && [ -d ${G}/tl ] || die "GEMINI++ data tables missing in ${G}"
find ${G} -mindepth 1 -maxdepth 1 ! -name tbl ! -name tl -exec rm -rf {} + || die "prune ${G}"

# GENIE (NucleusGenINCL.xml) looks for the ABLA++ tables in
# ${INCLXX_DATA_DIR}/de-excitation/ablaxx/upstream/data/G4ABLA3.0/, but
# INCL++ installs them together with its own tables in share/data.
mkdir -p share/de-excitation/ablaxx/upstream/data || die "mkdir ablaxx data"
ln -s ../../../../data share/de-excitation/ablaxx/upstream/data/G4ABLA3.0 || die "ln G4ABLA3.0"
[ -f share/de-excitation/ablaxx/upstream/data/G4ABLA3.0/frldm.dat ] || die "ABLA++ tables not reachable"
[ -d share/de-excitation/abla07/upstream/tables ] || die "ABLA07 tables missing"

cp -p "${SRC}/LICENSE.pdf" . 2>/dev/null

cat > BUILD_INFO <<EOF
INCL++ ${INCLXX_VERSION} (GENIE-interface fork) -- binary external of genie -q ${CQUAL}:incl634:${BTYPE}
source   : github.com/LiangLiu212/inclxx @ ${COMMIT}
built    : $(date -u +"%Y-%m-%d %H:%M UTC") on $(hostname) ($(cat /etc/redhat-release))
compiler : gcc $(gcc -dumpfullversion) (${CQUAL}), CMAKE_BUILD_TYPE=RelWithDebInfo (-O2 -g -DNDEBUG)
depends  : root v6_28_12 -q ${CQUAL}:p3915:${BTYPE}, boost v1_82_0 -q ${CQUAL}:${BTYPE}
recipe   : build_inclxx_external.sh (cmake options below)
contents : bin/inclxx-config bin/thisinclxx.sh lib/ include/ share/ (run-time data only).
           No source code is distributed; see LICENSE.pdf for the INCL++ licence.
cmake    : ${CMAKE_OPTS[*]/#-DCMAKE_INSTALL_PREFIX=*/}
EOF

# Nothing but headers may be left of the source code.
LEFT=$(find . -type f \( -name '*.cc' -o -name '*.cpp' -o -name '*.cxx' -o -name '*.c' -o -name '*.f' \
        -o -name '*.F' -o -name '*.f90' -o -name '*.inc' -o -name '*.in' -o -name 'CMakeLists.txt' -o -name 'Makefile' \))
[ -z "${LEFT}" ] || die "source files left in the bundle:\n${LEFT}"
[ -z "$(find share -type f \( -name '*.h' -o -name '*.hh' \))" ] || die "headers left under share/"
[ "$(ls lib/*.so | wc -l)" = 9 ] || die "expected 9 shared libraries in lib/"
[ -x bin/inclxx-config ] || die "bin/inclxx-config missing"
[ -f include/G4INCLVersion.hh ] || die "include/G4INCLVersion.hh missing"

echo "===== $(date) : replacing ${HERE}/inclxx ====="
rm -rf "${HERE}/inclxx" && cp -a "${PREFIX}" "${HERE}/inclxx" || die "cannot install ${HERE}/inclxx"
cd "${HERE}/inclxx" || die "cd ${HERE}/inclxx"
du -sh *
echo "inclxx-config --dflags: $(bin/inclxx-config --dflags)"
echo "inclxx-config --libs  : $(bin/inclxx-config --libs)"
grep G4INCL "${SRC}/utils/include/G4INCLVersion.hh"
echo "===== $(date) : INCL++ external build DONE rc=0 ====="
