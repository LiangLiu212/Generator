function build_genie() {
  prep_build genie ${1} ${2}:${build_type} || return 0
  ./build_genie.sh ${product_topdir} ${2} ${build_type} ${maketar} >& "${logfile}"
}

# Local Variables:
# mode: sh
# eval: (sh-set-shell "bash")
# End:
