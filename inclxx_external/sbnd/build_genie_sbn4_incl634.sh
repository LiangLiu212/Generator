#!/bin/bash
########################################################################
# Build  genie v3_06_02_sbn4 -q e26:incl634:prof  as a UPS product, with
# the INCL++ v6.34 binary external of this branch installed inside it.
#
#   bootstrap.sh <product-dir>
#       GENIE source (this clone, checked-out branch, committed files only)
#       + Reweight source -> <product-dir>/genie/v3_06_02_sbn4/tar
#       NOTE: it first deletes <product-dir>/genie/v3_06_02_sbn4{,.version}
#   build_genie.sh <product-dir> e26:incl634 prof tar
#       build, install, declare and make the binary tarball
#       <product-dir>/genie-3.06.02.sbn4-sl7-x86_64-e26-incl634-prof.tar.bz2
#
# Usage:  ./build_genie_sbn4_incl634.sh <product-dir> > build.log 2>&1
#         <product-dir> = a writable UPS products area (it must contain .upsfiles/)
#         Enters the SL7 container itself; about 25 min on one core.
# Env:    GENIE_SOURCE_URL / GENIE_SOURCE_REF   build another repository / branch
#         SBNDCODE_VERSION                      (default v10_14_02_05)
########################################################################
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/common.sh"
enter_sl7 "${HERE}/$(basename "${BASH_SOURCE[0]}")" "$@"

[ -d "$1/.upsfiles" ] || { echo "usage: $0 <product-dir>   (a UPS products area containing .upsfiles/)"; exit 1; }
PROD="$(cd "$1" && pwd -P)"

echo "===== $(date) : genie ${PKG_VERSION} -q ${PKG_QUAL} build START ($(cat /etc/redhat-release)) ====="
setup_sbnd_no_genie
export PRODUCTS=${PROD}:${PRODUCTS}
echo "GENIE='${GENIE}' INCLXX_DIR='${INCLXX_DIR}' (empty=good)  gcc=$(which gcc)"

echo "===== $(date) : STEP 1/2 bootstrap.sh ====="
cd "${HERE}" || exit 1
./bootstrap.sh ${PROD} || { echo "===== $(date) : BOOTSTRAP FAILED ====="; exit 3; }
echo "===== $(date) : STEP 2/2 build_genie.sh ${PROD} e26:incl634 prof tar ====="
./build_genie.sh ${PROD} e26:incl634 prof tar || { echo "===== $(date) : BUILD FAILED ====="; exit 4; }
echo "===== $(date) : BUILD DONE rc=0 ====="

FQ=$(ls -d ${PROD}/genie/${PKG_VERSION}/Linux64bit+*-e26-incl634-prof)
grep -n "GOPT_ENABLE_PYTHIA\|GOPT_ENABLE_INCL\|GOPT_WITH_INCL\|GOPT_WITH_BOOST" ${FQ}/GENIE-Generator/src/make/Make.config
ls -d ${FQ}/inclxx/{bin,lib,include,share}
ls -l ${PROD}/${PKG_TARBALL}
echo "inclxx files in the binary tarball: $(tar -tjf ${PROD}/${PKG_TARBALL} | grep -c '/inclxx/')"
