#!/bin/bash
########################################################################
# Test  genie v3_06_02_sbn4 -q e26:incl634:prof  the way it will be used from
# cvmfs: unpack the BINARY TARBALL into an empty products area, set it up from
# there and run gevgen. strace proves that INCL++ is read from inside the
# product only (not from the build area or another INCL++ installation).
#
# Usage:  ./test_genie_sbn4_incl634.sh <product-dir> [<scratch-dir>]
#         <product-dir>  where build_genie_sbn4_incl634.sh left the tarball
#         <scratch-dir>  is created and DELETED first (default <product-dir>_reloc_test, ~0.7 GB)
# Env:    XSEC_PROD  products area with the INCL26 genie_xsec splines
#         XSEC_QUAL  genie_xsec qualifier   (default INCL2609a00000:k250:e1000)
#         NEV        number of events       (default 20)
########################################################################
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/common.sh"
enter_sl7 "${HERE}/$(basename "${BASH_SOURCE[0]}")" "$@"

[ -d "$1/.upsfiles" ] || { echo "usage: $0 <product-dir> [<scratch-dir>]"; exit 1; }
PROD="$(cd "$1" && pwd -P)"
RELOC=${2:-${PROD}_reloc_test}
XSEC_PROD=${XSEC_PROD:-/exp/sbnd/app/users/liangliu/sbnd_genie/localProducts_larsoft_v10_14_02_02_e26_prof}
XSEC_QUAL=${XSEC_QUAL:-INCL2609a00000:k250:e1000}
NEV=${NEV:-20}

echo "===== $(date) : relocation test of ${PROD}/${GENIE_TARBALL} in ${RELOC} ====="
[ -f "${PROD}/${GENIE_TARBALL}" ] || fail "no ${PROD}/${GENIE_TARBALL}"
rm -rf "${RELOC}"; mkdir -p "${RELOC}/run" || fail "cannot create ${RELOC}"
cp -a ${PROD}/.upsfiles "${RELOC}/" || fail "cannot copy .upsfiles"
tar -xjf ${PROD}/${GENIE_TARBALL} -C "${RELOC}" || fail "cannot unpack the tarball"

setup_sbnd_no_genie
export PRODUCTS=${RELOC}:${PRODUCTS}
setup genie ${GENIE_VERSION} -q ${GENIE_QUAL} || fail "setup genie ${GENIE_VERSION} -q ${GENIE_QUAL}"
setup genie_xsec v3_06_00 -q ${XSEC_QUAL} -z ${XSEC_PROD} || fail "setup genie_xsec ${XSEC_QUAL} -z ${XSEC_PROD}"
ups active | grep -E "^genie|^boost|^root |^pythia|^lhapdf"
for v in GENIE GENIE_FQ_DIR INCLXX_DIR INCLXX_DATA_DIR GENIE_XSEC_TUNE; do echo "$v=${!v}"; done
case ${GENIE}:${INCLXX_DIR}:${INCLXX_DATA_DIR} in
  ${RELOC}/*:${RELOC}/*:${RELOC}/*) ;;
  *) fail "environment does not point into ${RELOC}" ;;
esac
[ "${INCLXX_DIR}" = "${GENIE_FQ_DIR}/inclxx" ] || fail "INCLXX_DIR is not \${GENIE_FQ_DIR}/inclxx"

echo "--- INCL++ libraries seen by libGPhHadTransp.so"
ldd ${GENIE_LIB}/libGPhHadTransp.so | grep -E "INCL|ABLA|DeExc|GEMINI|SMM|FERMI" | awk '{print "    "$1" => "$3}'
[ "$(ldd ${GENIE_LIB}/libGPhHadTransp.so | grep -E "INCL|ABLA|DeExc|GEMINI|SMM|FERMI" | grep -c "${INCLXX_DIR}/lib/")" = 9 ] \
  || fail "the 9 INCL++ libraries do not all resolve inside the product"
# (GENIE libraries do not link each other, so 'ldd -r' always lists undefined genie:: symbols;
#  only missing libraries and unresolved INCL++ / ABLA symbols are errors here)
[ -z "$(ldd -r ${GENIE_LIB}/libGPhHadTransp.so 2>&1 | grep -E "not found|undefined symbol: .*(G4INCL|G4Abla|ABLA)")" ] \
  || fail "missing libraries / unresolved INCL++ symbols"

echo "--- Make.config (paths must be relative to UPS variables)"
grep -n "GOPT_WITH_INCL\|GOPT_WITH_BOOST\|GOPT_ENABLE_INCL=" ${GENIE}/src/make/Make.config
grep -q 'GOPT_WITH_INCL_LIB=${GENIE_FQ_DIR}/inclxx/lib' ${GENIE}/src/make/Make.config || fail "Make.config INCL++ path is not relocatable"

cd "${RELOC}/run" || fail "cd run"
echo "--- gevgen: ${NEV} numu Ar40 events at 1 GeV, tune ${GENIE_XSEC_TUNE}"
strace -f -e trace=open,openat,execve -o strace.out \
  gevgen -n ${NEV} -p 14 -t 1000180400 -e 1 --seed 20260824 \
         --tune ${GENIE_XSEC_TUNE} --cross-sections ${GENIEXSECFILE} --event-generator-list Default \
         -o test.ghep.root > gevgen.log 2>&1
rc=$?
echo "    gevgen rc=${rc}   FATAL lines: $(grep -c FATAL gevgen.log)"
[ ${rc} = 0 ] || fail "gevgen rc=${rc} (see ${RELOC}/run/gevgen.log)"
gevdump -f test.ghep.root 2>&1 | grep -E "^ \|" > test.dump
echo "    particle lines: $(wc -l < test.dump)"
[ -s test.dump ] || fail "empty event dump"

echo "--- INCL++ / genie files opened outside the relocated product (must be none)"
BAD=$(grep -v ENOENT strace.out | grep -o '"/[^"]*"' | grep -E '/(inclxx[^/]*|inclpp|genie/v[0-9][^/]*)/' | grep -v "^\"${RELOC}/" | sort -u)
[ -z "${BAD}" ] || { echo "${BAD}" | head -20; fail "INCL++ / genie files were read from outside the relocated product"; }
echo "    INCL++ libraries and data tables read from the product:"
grep -v ENOENT strace.out | grep -o "\"${INCLXX_DIR}/[^\"]*\.\(so\|dat\|tab\|tbl\|tl\|inv\)\"" | sort -u | sed "s%\"${INCLXX_DIR}/%      inclxx/%; s%\"$%%; s%//%/%"

echo "RESULT: PASS -- relocated tarball is self-contained (INCL++ comes from \${GENIE_FQ_DIR}/inclxx)"
