#!/bin/bash
########################################################################
# bootstrap.sh
#
# Build a source code distribution for genie.
#
# Adapted from bootstrap.sh.template provided by ssibuildshims v2_09_02.
#
########################################################################

########################################################################
# Usage instructions.
#
# 1. Obtain the repository:
#      git clone ssh://p-build-framework@cdcvs.fnal.gov/cvs/projects/build-framework-genie-ssi-build
# 2. cd to the correct place and execute this bootstrap.sh to generate
#    the source distribution archive:
#      cd build-framework-genie-ssi-build && \
#        ./bootstrap.sh <ups-top-dir>
########################################################################

####################################
# Useful variables
prog=${BASH_SOURCE[0]##*/}
thisdir=$(cd "${BASH_SOURCE%/*}" && pwd -P)
trace=maybe_xtrace # Could set to xtrace instead for unconditional.

####################################
# Useful functions.

# Report usage.
usage() {
  cat <<EOF
USAGE: ${prog} ${ssi_msg_config_usage:-<ssi-msg-config-options>} [--] <product-dir>
       ${prog} -h|--help [<product-dir>]

  Install the source code distribution for genie in <product-dir> and
  create the source distribution archive.

EOF
  (( ${1:-0} < 2 )) && return;
  cat <<EOF
OPTIONS

  -h|--help

    This help.

${ssi_msg_config_options:-"  <ssi-msg-config-options>

    [Supply <product-dir> to enable more details]
"}\
ARGUMENTS

  <product-dir>

    The top-level UPS directory under which the package should be
    installed.

ENVIRONMENT

${ssi_msg_config_environment}
EOF
}

# Local (dumb) ssi_die() function for use until we source ssi_functions.
ssi_die() {
  printf "FATAL ERROR: ${*}\n" 1>&2; exit 1
}

# Ensure we have access to the ssibuildshims product we need.
ensure_ssibuildshims() {
  product_dir="${@: -1}"
  # Try as hard as we can to give useful information if necessary or
  # requested.
  if ! (( $# )) || [ "${product_dir}" = -h ] || [ "${product_dir}" = --help ]; then
    # We\'ve been asked for help: give as much as we can.
    usage 2; exit 1
  elif [[ "${product_dir}" == -* ]]; then
    # Bad usage.
    ssi_die "Bad <product-dir> ${product_dir}\n$(usage)"
  elif ! { [ -d "${product_dir}/.upsfiles" ] && [ -w "${product_dir}" ]; }; then
    # Something wrong with <product-dir>.
    ssi_die "<product-dir> not a valid, writable unified UPS directory: ${product_dir}"
  fi
  # Promising...
  local abs_product_dir="$(cd "${product_dir}" 2>/dev/null 2>&1 && pwd -P)"
  (( $? )) && ssi_die "unable to cd to <product-dir>: ${product_dir}"
  # Absolute is safer.
  product_dir="${abs_product_dir}"
  # Find our own bootstraps...
  local ssibuildshims_topdir="${product_dir}/ssibuildshims/${ssibuildshims_version}"
  if ! [ -f "${ssibuildshims_topdir}/bin/ssi_functions" ]; then
    echo "Installing ssibuildshims ${ssibuildshims_version} in ${product_dir}"
    local ssidotver=${ssibuildshims_version//_/.}; ssidotver=${ssidotver#v}
    curl --fail -o "${product_dir}/ssibuildshims-${ssidotver}-noarch.tar.bz2" \
      "https://scisoft.fnal.gov/scisoft/packages/ssibuildshims/${ssibuildshims_version}/ssibuildshims-${ssidotver}-noarch.tar.bz2" || \
      ssi_die "unable to obtain ssibuildshims ${ssibuildshims_version} from https://scisoft.fnal.gov/"
    tar -C "${product_dir}" -xf "${product_dir}/ssibuildshims-${ssidotver}-noarch.tar.bz2" || \
      ssi_die "unable to expand ssibuildshims package ssibuildshims-${ssidotver}-noarch.tar.bz2"
  fi
  # ...and yank \'em.
  if source \
    "${ssibuildshims_topdir}/bin/ssi_functions" \
    >/dev/null 2>&1; then
    ssi_info "found ssibuildshims ${ssibuildshims_version} at ${ssibuildshims_topdir}"
  else
    ssi_die "unable to verify, install and initialize ssibuildshims
${ssibuildshims_version} in product directory ${product_dir}"
  fi
}

########################################################################
# Main.
########################################################################

####################################
# Package name and version information.
package=genie
origpkgver=v3_06_02_sbn4
pkgver=${origpkgver}
ssibuildshims_version=v2_16_03
geniedir=GENIE-Generator
reweightdir=GENIE-Reweight
# GENIE source: by default the Generator clone this recipe sits in, at the
# branch that is checked out there. Only COMMITTED files are built (git
# archive), and they include the INCL++ binary external inclxx_external/inclxx.
# Override with GENIE_SOURCE_URL / GENIE_SOURCE_REF (branch or tag), e.g.
#   GENIE_SOURCE_URL=https://github.com/LiangLiu212/Generator.git
#   GENIE_SOURCE_REF=feature/genie-incl634-sbnd
gitinfo() { ( cd "${thisdir}" && git rev-parse "$@" 2>/dev/null ); }
sourceurl=${GENIE_SOURCE_URL:-file://$(gitinfo --show-toplevel)}
geniever=${GENIE_SOURCE_REF:-$(gitinfo --abbrev-ref HEAD)}
case "${sourceurl}:${geniever}" in
  file://:*|*:|*:HEAD) ssi_die "cannot tell which GENIE source to build: set GENIE_SOURCE_URL and GENIE_SOURCE_REF (a branch or tag)";;
esac
reweighturl=${REWEIGHT_SOURCE_URL:-https://github.com/LiangLiu212/Reweight.git}
reweightver=${REWEIGHT_SOURCE_REF:-master}

####################################
# Handle options and arguments.

# Obtain ssibuildshims and set product_dir as early as possible.
ensure_ssibuildshims "${@}"

# Initialization and defaults.
(( ssi_trace = 1 ))

# Parse options.
while (( $# )); do
  case $1 in
    -h|--help) usage 2; exit 1;;
    --) shift; break;;
    -*) ssi_handle_config_options "${1}" || \
      ssi_die 2 "unrecognized option $1\n$(usage)";;
    *) shift; break # Swallow <product-dir>: dealt with.
  esac
  shift
done

# Environment overrides.
ssi_honor_environment_config

# Useful shorthand.
pkgdir="${product_dir}/${package}/${pkgver}"

# Report.
ssi_info "making source code distribution for ${package} ${pkgver} in ${pkgdir}"
ssi_info "GENIE source: ${sourceurl} @ ${geniever}; Reweight: ${reweighturl} @ ${reweightver}"

########################################################################
# Package-specific operations.

echo "calling make_source_code_base"
${trace} make_source_code_base "${ssi_args[@]}" -e tar \
  ${product_dir} \
  ${package} \
  ${pkgver} \
  ${thisdir}

# now checkout and archive the source code
pkgdir=${product_dir}/${package}/${pkgver}
if [ ! -d ${pkgdir}/tar ]
then
   echo "ERROR: cannot find ${pkgdir}/tar"
   exit 1
fi

cd ${pkgdir}/tar || exit 1

# Assemble the distribution, keeping the user informed.

# GENIE source
${trace} download_file -F ${geniedir}.tar -G ${geniever} "${sourceurl}" || \
  ssi_die "unable to obtain GENIE ${geniever} source"

# Reweight
${trace} download_file  -F ${reweightdir}.tar -G ${reweightver} "${reweighturl}" || \
  ssi_die "unable to obtain Reweight ${reweightver} source"

# The INCL++ binary external must have come with the GENIE source.
tar -tf ${geniedir}.tar ${geniedir}/inclxx_external/inclxx/bin/inclxx-config >/dev/null 2>&1 || \
  ssi_die "${geniever} of ${sourceurl} has no inclxx_external/inclxx (INCL++ binary external)"

# Make the archive distribution file.
${trace} make_source_code_tarball "${ssi_args[@]}"  ${product_dir} ${package} ${pkgver}
