#!/bin/bash
########################################################################
# LArSoft check of  genie v3_06_02_sbn4 -q e26:incl634:prof  with sbndcode:
# nugen must find all its GENIE libraries in the new product, and a GENIEGen
# job (rockbox: GENIE + Geant4) must run.
#
# lar loads Geant4, which carries its own INCL++ with the same symbol names,
# so the INCL++ libraries shipped inside genie are put in LD_PRELOAD.
#
# Usage:  ./test_lar_sbn4_incl634.sh <product-dir> [<run-dir>]
#         <product-dir>  products area holding genie v3_06_02_sbn4
#         <run-dir>      is created and DELETED first (default <product-dir>_lar_test)
# Env:    XSEC_PROD / XSEC_QUAL  INCL26 genie_xsec splines (see test_genie_sbn4_incl634.sh)
#         FCL   fcl file to run   (default: Liang's rockbox INCL26 example)
#         NEV   number of events  (default 1)
########################################################################
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/common.sh"
enter_sl7 "${HERE}/$(basename "${BASH_SOURCE[0]}")" "$@"

[ -d "$1/.upsfiles" ] || { echo "usage: $0 <product-dir> [<run-dir>]"; exit 1; }
PROD="$(cd "$1" && pwd -P)"
RUN=${2:-${PROD}_lar_test}
XSEC_PROD=${XSEC_PROD:-/exp/sbnd/app/users/liangliu/sbnd_genie/localProducts_larsoft_v10_14_02_02_e26_prof}
XSEC_QUAL=${XSEC_QUAL:-INCL2609a00000:k250:e1000}
FCL=${FCL:-/exp/sbnd/data/users/liangliu/genie_incl_package/prodgenie_incl26_rockbox_sbnd.fcl}
NEV=${NEV:-1}

echo "===== $(date) : lar test of genie ${GENIE_VERSION} -q ${GENIE_QUAL} from ${PROD} ====="
setup_sbnd_no_genie
export PRODUCTS=${PROD}:${PRODUCTS}
setup genie ${GENIE_VERSION} -q ${GENIE_QUAL} || fail "setup genie"
setup genie_xsec v3_06_00 -q ${XSEC_QUAL} -z ${XSEC_PROD} || fail "setup genie_xsec"
ups active | grep -E "^genie|^sbndcode|^nugen|^geant4 "
echo "GENIE=${GENIE}"; echo "INCLXX_DIR=${INCLXX_DIR}"; echo "GENIE_XSEC_TUNE=${GENIE_XSEC_TUNE}"

L=${NUGEN_LIB}/libnugen_EventGeneratorBase_GENIE.so
echo "--- ${L}"
[ -z "$(ldd ${L} 2>&1 | grep "not found")" ] || { ldd ${L} 2>&1 | grep "not found" | head; fail "nugen does not find all its libraries"; }
echo "    GENIE + INCL++ libraries taken from the new product: $(ldd ${L} | grep -c "${GENIE_FQ_DIR}/")"
[ -z "$(ldd ${L} | grep "/genie/" | grep -v "${GENIE_FQ_DIR}/")" ] || fail "nugen picks up libraries of another genie"
echo "    undefined symbols (ldd -r): $(ldd -r ${L} 2>&1 | grep -c "undefined symbol")"

rm -rf "${RUN}"; mkdir -p "${RUN}" && cd "${RUN}" || fail "cannot create ${RUN}"
LD_PRELOAD=$(incl_preload) lar --print-description GENIEGen > print_description.log 2>&1 || fail "lar --print-description GENIEGen"
echo "--- lar --print-description GENIEGen: OK"

echo "--- $(date) : lar -c $(basename ${FCL}) -n ${NEV}  (INCL++ libraries in LD_PRELOAD)"
LD_PRELOAD=$(incl_preload) lar -c ${FCL} -n ${NEV} > lar.log 2>&1
rc=$?
echo "    $(date) : lar rc=${rc}"
grep -E "Tune configured|HadronTransp-Model|TrigReport Events total|Art has completed" lar.log | sort | uniq -c | sort -rn | head -8
[ ${rc} = 0 ] || { tail -30 lar.log; fail "lar rc=${rc} (see ${RUN}/lar.log)"; }
ls -lh *.root
echo "RESULT: PASS -- lar GENIEGen job ran with genie ${GENIE_VERSION} -q ${GENIE_QUAL}"
