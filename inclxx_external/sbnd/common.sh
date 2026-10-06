# Shared by the build / test drivers in this directory (sourced, not run).
#
# Environment (all optional):
#   SBNDCODE_VERSION  sbndcode release whose toolchain is used   (default v10_14_02_05)
#   SL7_IMAGE         SL7 container image
SBNDCODE_VERSION=${SBNDCODE_VERSION:-v10_14_02_05}
SL7_IMAGE=${SL7_IMAGE:-/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:jsl}
GENIE_VERSION=v3_06_02_sbn4
GENIE_QUAL=e26:incl634:prof
GENIE_TARBALL=genie-3.06.02.sbn4-sl7-x86_64-e26-incl634-prof.tar.bz2

fail() { echo "RESULT: FAIL -- $*"; exit 1; }

# The e26 products only run on SL7: re-run the calling script in the container.
enter_sl7() {  # enter_sl7 <script> [args...]
  grep -q "release 7" /etc/redhat-release 2>/dev/null && return 0
  local appt=/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer
  [ -x "${appt}" ] || appt=$(command -v apptainer)
  [ -n "${appt}" ] || { echo "FATAL: not on SL7 and no apptainer found"; exit 1; }
  exec "${appt}" exec --pid --ipc -B /etc/hosts,/tmp,/opt,/cvmfs,/exp/,/pnfs/,/nashome \
    "${SL7_IMAGE}" /bin/bash "$@"
}

# sbndcode environment without its own genie / genie_xsec.
setup_sbnd_no_genie() {
  source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
  setup sbndcode ${SBNDCODE_VERSION} -q prof:e26 || { echo "FATAL: cannot set up sbndcode ${SBNDCODE_VERSION}"; exit 1; }
  unsetup genie
  unsetup genie_xsec
}

# The 9 INCL++ libraries, ':'-separated, for LD_PRELOAD in lar jobs.
incl_preload() {
  ls ${INCLXX_DIR}/lib/lib{INCL_Utils,INCL_IO,INCL_Physics,DeExcitation,ABLA07,ABLAXX,FERMI_BREAKUP,GEMINIXX,SMM}.so | tr '\n' ':' | sed 's/:$//'
}
