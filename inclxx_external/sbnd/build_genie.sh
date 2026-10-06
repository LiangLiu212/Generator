#!/bin/bash
########################################################################
# build_genie.sh
#
# Build genie and package for UPS.
#
# Adapted from build.sh.template provided by ssibuildshims v2_09_02.
########################################################################

########################################################################
# Non-customize-able boilerplate: skip to "Customize-able functions."
########################################################################

####################################
# Useful variables.
prog=${BASH_SOURCE[0]##*/}

##################
# Local (dumb) ssi_die() function for use until we source ssi_functions.
ssi_die() {
  printf "FATAL ERROR: ${*}\n" 1>&2; exit 1
}
##################

##################
# Basic usage function for use by ensure_buildshims until we have enough
# information to provide something more specific. Do NOT customize!
build_script_usage() {
    cat <<EOF
USAGE: ${prog} <ssi-msg-config-options> [--] <product-dir> <qualifier-args> <maketar>
       ${prog} -h|--help [<product-dir>]

  Build ${package} and package for UPS.

OPTIONS

  -h|--help

    Help. Supply <product-dir> for specific, detailed information on
    invoking ${prog}.

  <ssi-msg-config-options>

    Supply <product-dir> for more information.

ARGUMENTS

  <product-dir>

    The top-level UPS directory in which the source code distribution
    was unpacked. When specified with -h|--help, this enables us to
    provide more specific help about invoking ${prog}.

  <qualifier-args>
  <maketar>

    Supply <product-dir> for more information.

EOF
}
##################

##################
# Make sure we have UPS and ssibuildshims set up.
ensure_ssibuildshims() {
  if ! (( $# )) || [ "${@: -1}" = "-h" ] || [ "${@: -1}" = "--help" ]; then
    build_script_usage
    exit 1
  fi
  # Find <product-dir> from our arguments:
  local arg
  for arg in "${@}"; do
    [[ "$arg" ==  -* ]] && continue;
    product_dir="${arg}"; break
  done
  if [ -z "${product_dir}" ]; then
    # Bad usage.
    ssi_die "Bad <product-dir> ${product_dir}\n$(build_script_usage)"
    exit 2
  elif ! { [ -d "${product_dir}/.upsfiles" ] && [ -w "${product_dir}" ]; }; then
    # Something wrong with <product-dir>.
    ssi_die "<product-dir> not a valid, writable unified UPS directory: " \
      "${product_dir}\n$(build_script_usage)"
  fi
  # Promising...
  local abs_product_dir="$(cd "${product_dir}" 2>/dev/null 2>&1 && pwd -P)"
  (( $? )) && ssi_die "unable to cd to <product-dir>: ${product_dir}\n" \
    "$(build_script_usage)"
  # Find our own bootstraps and yank 'em.
  local ssibuildshims_topdir="${product_dir}/ssibuildshims/${ssibuildshims_version}"
  [ -f "${ssibuildshims_topdir}/bin/ssi_functions" ] && \
    source "${ssibuildshims_topdir}/bin/ssi_functions" >/dev/null 2>&1 && \
    ensure_ups "${product_dir}" && \
    setup ssibuildshims ${ssibuildshims_version} -z "${product_dir}" || \
    ssi_die "unable to find ssibuildshims ${ssibuildshims_version} in
${product_dir} and set it up with UPS."
  # Absolute is safer.
  absolute_path -- "${product_dir}" product_dir
}
##################

########################################################################
# Customize-able functions.
########################################################################

##################
# Ascertain whether we built and installed correctly. Add more /
# different calls to verify_in_filesystem as appropriate. Arguments
# (e.g. -q) are passed through to verify_in_filesystem.
is_built() { # *C*
  ${trace} verify_in_filesystem -d "${@}" "${pkgdir}/lib"{64,}
}
##################

##################
# If a distribution archive is wanted and the installation is successful
# according to is_built(), then make it. Any arguments (e.g. -q) are
# passed through both to is_built() and make_basic_distribution().
#
# If make_basic_distribution() finds that the product is declared
# correctly and attempts to make the distribution, then it will cause us
# to exit with the status code of that operation. If not, it returns an
# error code.
#
# See make_basic_distribution -h for useful options to add to
# mbd_opts(). *C*
maybe_make_distribution() {
  if [ "${maketar}" = "tar" ]; then
    # Customize mbd_opts. *C*
    local mbd_opts=()
    ${trace} is_built "$@" && \
      ${trace} make_basic_distribution "${@}" "${mbd_opts[@]}"
  fi
}
##################

########################################################################
# Package-specific functions.
########################################################################

# ...

########################################################################
# Main.
########################################################################

####################################
# Useful variables. *C*
trace=maybe_xtrace # Could set to xtrace instead for unconditional.
(( ssi_trace = 1 )) # Trace certain commands by default.
####################################

####################################
# Package name and version information. *C*
package=genie
origpkgver=v3_06_02_sbn4
pkgver=${origpkgver}
ssibuildshims_version=v2_16_03
pkgdotver=${origpkgver//_/.}; pkgdotver=${pkgdotver#v}
geniedot=`echo ${origpkgver} | sed -e 's/_b/b/g' | sed -e 's/_/./g'`
geniever=`echo ${origpkgver} | sed -e 's/_b/b/' | sed -e 's/^v/R-/g'`
reweightver=R-1_04_00
geniedir=GENIE-Generator
reweightdir=GENIE-Reweight
####################################

##################
# Other useful variables. *C*
pkgtarname=${package}-${pkgdotver}
pkgtarfile=${pkgtarname}.tar.gz
##################

####################################
# Handle options and arguments.

# Set up ssibuildshims and UPS / product_dir as early as possible.
ensure_ssibuildshims "${@}"

##################
# Now parse options and arguments, setting basequal, extraqual, maketar,
# etc.
#
# Use verify_quals() to parse command line options, verify that the user
# has provided qualifiers supported by this build script, and set
# expected variables such as basequal, extraqual, cqual, etc.
#
# Customize the extra usage text as appropriate. Note that DEBUG_PATCHES
# and VERBOSE_PATCHES are honored by prepare_source() and
# patch_source(). DEBUG_(CONFIG|BUILD|TEST|INSTALL) should be accounted
# for (or not) by the relevant command in this script.
#
# See verify_quals -h for options and examples.
verify_quals \
  --cqual={e20,e26,c7,c14} \
  --bq-addl={,geant4,incl634} \
  --btype={debug,prof} \
  --parse --usage-file=- -- "$@" <<EOF || exit # *C*
ENVIRONMENT

  DEBUG_(CONFIG|BUILD|TEST|INSTALL)

    Verbose output from the corresponding step.

  DEBUG_PATCHES
  VERBOSE_PATCHES

    Command tracing and increased verbosity for source unpacking and patching
    operations. DEBUG_PATCHES will cause exit after patches have been applied.
EOF
##################

##################
# Define globals: pkgdir, sourcedir, blddir, fullqual, etc. See
# define_basics -h for steering options.
${trace} define_basics NO_DEBUG_SOURCE TRACK_PROGRESS || exit # *C*
##################

####################################
# Declare ourselves privately, and then set ourselves up to allow setup
# of dependencies. *C*
#

# The INCL++ external only exists as an e26:prof binary.
case ${addl}:${cqual}:${extraqual} in
   incl634:e26:prof) ;;
   incl634:*) ssi_die "Qualifier incl634 is only supported as e26:incl634 prof." ;;
esac

# Comment or delete if not required:
${trace} fake_declare_product "${ssi_args[@]}" -x -- || exit
####################################

####################################
# Build-only dependencies. *C*
#
# pythia6 is NOT pulled in by genie.table (only pythia8 is), but v3_06_02
# still builds the pythia6 interface (GOPT_ENABLE_PYTHIA6=YES), so set it up
# here to provide $PYLIB. pythia6 is prof-only.
${trace} setup pythia v6_4_28x -q ${cqual}:prof || ssi_die "failed to setup pythia6"

# INCL++ (qualifier incl634) is NOT a UPS product: INCL++ v6.34 comes with the
# GENIE source as a binary external (libraries + headers + run-time data, no
# source) in ${geniedir}/inclxx_external/inclxx. It is moved below to
# ${pkgdir}/inclxx, so it is installed and distributed as part of this genie
# product. genie.table points INCLXX_DIR etc. at ${GENIE_FQ_DIR}/inclxx.
if [ "${addl}" = "incl634" ]; then
  # boost is required by GENIE's --enable-incl (the INCL++ interface).
  ${trace} setup boost v1_82_0 -q ${cqual}:${extraqual} || ssi_die "failed to setup boost"
fi
####################################

# Print active products.
echo
ups active
echo

# Possible short-circuit if we're already built and declared.
${trace} maybe_make_distribution --silent

ssi_info "building ${package} for ${OS}-${plat}-${qualdir} (flavor ${pkgflvr})"

if [ -z ${GENIE} ]
then
   echo "ERROR: failed to setup genie"
   exit 1
else
   echo "GENIE: ${GENIE}"
   echo "GENIE_REWEIGHT: ${GENIE_REWEIGHT}"
fi
if [ -z ${PYLIB} ]
then
   echo "ERROR: failed to setup pythia"
   exit 1
else
   echo "PYLIB: ${PYLIB}"
fi
if [ -z ${ROOTSYS} ]
then
   echo "ERROR: failed to setup root"
   exit 1
else
   echo "ROOTSYS: ${ROOTSYS}"
fi
if [ -z ${LIBXML2_INC} ]
then
   echo "ERROR: failed to setup libxml"
   exit 1
else
   echo "LIBXML2_INC: ${LIBXML2_INC}"
fi
if [ -z ${LOG4CPP_INC} ]
then
   echo "ERROR: failed to setup log4cpp"
   exit 1
else
   echo "LOG4CPP_INC: ${LOG4CPP_INC}"
fi
if [ -z ${LHAPDF_INC} ]
then
   echo "ERROR: failed to setup lhapdf"
   exit 1
else
   echo "LHAPDF_INC: ${LHAPDF_INC}"
fi
# pythia8 is NOT used by this build (--disable-pythia8 below, 2026-09-11);
# genie.table still lists it as a runtime dependency, so just report it.
echo "PYTHIA8_FQ_DIR (unused, pythia8 disabled): ${PYTHIA8_FQ_DIR:-<not set>}"

# Ensure directories: pkgdir, sourcedir, blddir.
${trace} init_pkg_dirs || exit

####################################
# Prepare source for building. *C*

${trace} cd "${pkgdir}" \
  || ssi_die "unable to cd to source directory ${pkgdir}"

# Source tarballs and patches. See prepare_source -h / patch_source -h
# for options.
#${trace} prepare_source "${ssi_args[@]}" -P "${patchdir}" \
#  "${tardir}/${pkgtarfile}" \
#  || exit

${trace} tar -xf ${tardir}/${geniedir}.tar || ssi_die "failed.to unwind ${geniedir}.tar" || exit
${trace} tar -xf ${tardir}/${reweightdir}.tar || ssi_die "failed.to unwind ${reweightdir}.tar" || exit

if [ "${addl}" = "incl634" ]; then
  ${trace} rm -rf "${pkgdir}/inclxx"
  ${trace} mv "${pkgdir}/${geniedir}/inclxx_external/inclxx" "${pkgdir}/inclxx" || \
    ssi_die "no INCL++ binary external in ${geniedir}/inclxx_external/inclxx"
  INCLXX_INSTALL=${pkgdir}/inclxx
  [ -x ${INCLXX_INSTALL}/bin/inclxx-config ] || ssi_die "no INCL++ external in ${INCLXX_INSTALL}"
  [ "$(cd "${INCLXX_DIR}" 2>/dev/null && pwd -P)" = "$(cd "${INCLXX_INSTALL}" && pwd -P)" ] || \
    ssi_die "genie.table set INCLXX_DIR=${INCLXX_DIR}, expected ${INCLXX_INSTALL}"
fi

# No patches: your changes are already committed in the forked source pulled by bootstrap.sh.
${trace} cd "${pkgdir}" \
  || ssi_die "unable to cd to source directory ${pkgdir}"

# Bail early if DEBUG_PATCHES is set.
[ -n "${DEBUG_PATCHES}" ] && exit 1
####################################

####################################
# Prepare to build. *C*
#
# Set configure / compile flags, environment variables, etc.

if [ "${extraqual}" = "debug" ]
then
    extracommand="--enable-debug --with-optimiz-level=O0"
elif [ "${extraqual}" = "prof" ]
then
    extracommand="--disable-debug --with-optimiz-level=O3"
    cxxflg="-g -DNDEBUG -fno-omit-frame-pointer"
    # compile flags for optimiz-level=O3: g++ -c  -Wall -fPIC  -O3   -Wno-strict-aliasing -ffriend-injection
fi

case ${cqual} in
  e2[06]) cc=gcc; cxx=g++; fc=gfortran;
           cxxflg="${cxxflg} -std=c++17 -Wno-deprecated-declarations -Wno-unused-variable -Wno-maybe-uninitialized";;
  c7|c14) cc=clang; cxx=clang++; fc=gfortran;
      cxxflg="${cxxflg} -std=c++17 -Wno-deprecated-declarations -Wno-register -Wno-unused-variable -Wno-shadow";;
  *) ssi_die "Qualifier $cqual not recognized."
esac

# special builds
case ${addl} in
  '')  echo "configured with neither Geant4 nor INCL++"
    extra_config="--disable-incl"
    ;;
  geant4)
    extra_config="--disable-incl
                  --enable-geant4
                  --with-geant4-inc=${GEANT4_FQ_DIR}/include
                  --with-geant4-lib=${GEANT4_FQ_DIR}/lib64"
    ;;
  incl634)
    echo "INCL++ v6.34 enabled, bundled external at ${INCLXX_INSTALL}"
    extra_config="--enable-incl
                  --with-incl-inc=${INCLXX_INSTALL}/include
                  --with-incl-lib=${INCLXX_INSTALL}/lib
                  --with-boost-inc=${BOOST_FQ_DIR}/include
                  --with-boost-lib=${BOOST_FQ_DIR}/lib"
    ;;
  *) ssi_die -I "accepted unknown addl=${addl}"
esac

##################
# Configure parallelism. *C*
ncores=$(ncores)

# For cmake --build:
#
# [[ -n "${CMAKE_BUILD_PARALLEL_LEVEL}" ]] || \
#   export CMAKE_BUILD_PARALLEL_LEVEL=${ncores}
#
# For ctest:
#
# [[ -n "${CTEST_PARALLEL_LEVEL}" ]] || \
#   export CTEST_PARALLEL_LEVEL=${ncores}
##################
####################################

####################################
# Configure, build, and install.
#
# Use ${blddir}/${pkgtarname} for in-source builds, ${blddir} otherwise.

${trace} cd "${blddir}" || \
  ssi_die "unable to cd to build directory ${blddir}"

##################
# Configure. *C*

# e.g.
#   ${trace} ${sourcedir}/${pkgtarname}/configure \
#     --prefix="${pkgdir}" ... || ssi_die "configure failed"
# OR
#   ${trace} cmake <opts> ${DEBUG_CONFIG:+--verbose} \
#     "${sourcedir}/${pkgtarname}" || \
#     ssi_die "configure failed"
##################

##################
# Build. *C*

# e.g.
#   ${trace} cmake --build "${blddir}" ${DEBUG_BUILD:+--verbose} || \
#     ssi_die "build failed"

set -x

# build genie generator
cd ${GENIE} || ssi_die "cannot find ${GENIE}"

# configure generates .../src/make/Make.config
# Flags mirror the official genie v3_06_02_sbn2 build (its installed
# src/make/Make.config): lhapdf6, pythia6; EXCEPT pythia8 is disabled here
# (2026-09-11, pythia6-only build). INCL++ is enabled by the incl634
# qualifier only (see extra_config above).
./configure  --prefix=${pkgdir} ${extracommand} \
--with-compiler=${cc} \
--enable-gsl \
--enable-rwght \
--enable-lhapdf6 \
--disable-lhapdf5 \
--enable-pythia6 \
--disable-pythia8 \
--enable-atmo \
--enable-boosted-dark-matter \
--enable-dark-neutrino \
--enable-dylibversion \
--enable-event-library \
--enable-flux-drivers \
--enable-fnal \
--enable-geom-drivers \
--enable-gfortran \
--enable-neutral-heavy-lepton \
--enable-heavy-neutral-lepton \
--enable-nnbar-oscillation \
--enable-nucleon-decay \
--disable-lowlevel-msg \
--with-libxml2-inc=${LIBXML2_INC} \
--with-libxml2-lib=${LIBXML2_FQ_DIR}/lib \
--with-log4cpp-inc=${LOG4CPP_INC} \
--with-log4cpp-lib=${LOG4CPP_FQ_DIR}/lib \
--with-lhapdf6-inc=${LHAPDF_INC} \
--with-lhapdf6-lib=${LHAPDF_FQ_DIR}/lib \
--with-pythia6-lib=${PYLIB} \
${extra_config} || ssi_die "configure failed."

# you must pass extra compiler flags directly to make
make GOPT_WITH_CXX_USERDEF_FLAGS="${cxxflg}" || ssi_die "make failed."
make install || ssi_die "make install failed."

# build genie reweight
cd ${GENIE_REWEIGHT} || ssi_die "cannot find ${GENIE_REWEIGHT}"

# you must pass extra compiler flags directly to make
make GOPT_WITH_CXX_USERDEF_FLAGS="${cxxflg}" || ssi_die "make failed."
make install || ssi_die "make install failed."

#clean now
cd ${GENIE} || ssi_die "cannot find ${GENIE}"
make clean || ssi_die "make clean failed."
cd ${GENIE_REWEIGHT} || ssi_die "cannot find ${GENIE_REWEIGHT}"
make clean || ssi_die "make clean failed."


set +x
##################

##################
# Test. *C*

# e.g.
# (( WANT_TESTS )) && \
#   { ${trace} ctest --output-on-failure || ssi_die "tests failed"; }
##################

##################
# Install. *C*

# e.g.
#   ${trace} cmake --install "${blddir}" || ssi_die "install failed"
##################

ssi_info "finished building ${package} ${pkgver}"

####################################
# Pre-declaration processing. *C*
#
# Move things, install from contrib, clean build products from in-source
# builds, etc.
####################################

##################
# Declare and set up the installed package. *C*
#
# Add -c for current declaration if required.
${trace} declare_product -x -- || exit
##################

####################################
# Post-declaration processing. *C*
#
# Operations requiring the *final* declaration and setup (as opposed to
# the pre-build fake declaration and setup). This should be rare.

echo ""
echo "FIX PATHS in ${GENIE}/src/make/Make.config"
echo ""

filename=/tmp/geniefix-$$$$.sh

if [ -e ${filename} ]; then rm -f ${filename}; fi

echo "#!/bin/bash" > ${filename}
echo " set -x " >>  ${filename}
echo "mv ${GENIE}/src/make/Make.config ${GENIE}/src/make/Make.config.bak" >> ${filename}
echo "cat ${GENIE}/src/make/Make.config.bak | sed 's%${PYLIB}%\${PYLIB}%' | sed 's%${GENIE_FQ_DIR}%\${GENIE_FQ_DIR}%'  | sed 's%${LHAPDF_FQ_DIR}%\${LHAPDF_FQ_DIR}%' | sed 's%${LOG4CPP_FQ_DIR}%\${LOG4CPP_FQ_DIR}%' | sed 's%${LIBXML2_FQ_DIR}%\${LIBXML2_FQ_DIR}%' | sed 's%${BOOST_FQ_DIR:-NO_BOOST_FQ_DIR}%\${BOOST_FQ_DIR}%'> ${GENIE}/src/make/Make.config" >> ${filename}
echo " set +x " >>  ${filename}

chmod +x ${filename}
${filename}

if [ -e ${filename} ]; then rm -f ${filename}; fi

####################################

#################
# Verify installation and make distribution archive. *C*
#
# To customize verification, see is_built(), since it is also used by
# maybe_make_distribution(). See the documentation for those functions
# (above) for customization opportunities.

# Add a message to ssi_build_failed if appropriate.
${trace} is_built && ${trace} ssi_build_completed && \
  ssi_info "${package} is installed at ${pkgdir}" && \
  ${trace} maybe_make_distribution || \
  ssi_build_failed # *C*
########################################################################
# Nothing beyond this point.
########################################################################
