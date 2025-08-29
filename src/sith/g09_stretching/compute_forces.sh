#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH --exclusive
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool computes the forces in from a chk file. The output replaces the
string 'conopt' for 'forces'.

  -c  run in cascade.
  -f  <chk file> of the configuration that you want to compute the forces.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

compute_forces () {
  echo "construct Z-matrix from $1"
  opt_name=$1
  for_name=${opt_name//conopt/forces}
  for_name=${for_name%.chk}

  echo -e "here \n newzmat -ichk -ozmat -rebuildzmat -bmodel "$1" $for_name.com"

  newzmat -ichk -ozmat -rebuildzmat -bmodel "$1" $for_name.com || \
    { lnbck=$(search_last_bck ${1%.chk}) ; \
      newzmat -ichk -ozmat -rebuildzmat -bmodel "${1%.chk}-bck_$lnbck.chk" \
      $for_name.com || \
      fail "creating matrix" ;
    }
  sed -i "/#P bmk\/6-31+g opt(modredun/c %mem=60000MB\n%chk=$for_name\n%NProcShared=8\n#P bmk\/6-31+g force" $for_name.com
  verbose "executes g09 computation of forces for $for_name.com"
  g09 $for_name.com || fail "computing forces"
  formchk -3 $for_name.chk
}

# ----- definition of functions finishes --------------------------------------

# ----- general setup ---------------------------------------------------------
cascade='false'
while getopts 'f:cvh' flag; do
  case "${flag}" in
    c) cascade='true' ;;
    f) chkfile=${OPTARG} ;;
  
    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" Forces $verbose

verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

if $cascade
then
  load_modules || fail "loading modules"
fi

# ----- BODY ------------------------------------------------------------------

[[ -f $chkfile ]] || fail "$chkfile does not exist"

compute_forces "$chkfile"

finish
