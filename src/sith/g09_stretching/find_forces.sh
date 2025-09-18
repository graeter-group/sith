#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e


# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool computes the forces from a chk files that contains a given structure
and saves them in a directory called 'forces' that has to be previously
created.

  -c  run in cascade.
  -f  <file> chk file.
  -p  <pattern> pattern present in the chk files that will be replaced with the
      word 'forces'.

  -v  verbose.
  -h  prints this message.

Note
----

  Take care with the  files that already exist in the directory 'forces'. They
  may be overwritten. 
"
exit 0
}

compute_forces () {
  chk_name=$1
  for_name=${chk_name//${pattern}/forces}
  for_name=${for_name%.chk}
  verbose "construct Z-matrix for $1 into $for_name"
  echo newzmat -ichk -ozmat -rebuildzmat -bmodel "$1" $for_name.com
  newzmat -ichk -ozmat -rebuildzmat -bmodel "$1" $for_name.com || \
    {
      lnbck=$(search_last_bck ${1%.chk}) ; \
      newzmat -ichk -ozmat -rebuildzmat -bmodel "${1%.chk}-bck_$lnbck.chk" \
      $for_name.com || fail "Creating the matrix"
    }
  sed -i "1i %NProcShared=$n_processors" $for_name.com
  sed -i "1i %chk=$for_name" $for_name.com
  sed -i "s/opt(modredun,calcfc)/force/g" $for_name.com
  verbose "Executes gaussian computation of forces for $1"
  gaussian $for_name.com || fail "computing forces"
  sith extract_forces -f $for_name.log -c -v
}

# ----- definition of functions finishes --------------------------------------

# ---- set-up starts ----------------------------------------------------------
cascade='false'
directory='./'
pattern=''
verbose='false'
n_processors=1
while getopts 'cf:n:p:vh' flag; do
  case "${flag}" in
    c) cascade='true' ;;
    f) chkfile=${OPTARG} ;;
    n) n_processors=${OPTARG} ;;
    p) pattern=${OPTARG} ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" FIND_FORCES $verbose

verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

if $cascade
then
  load_modules
  n_processors=$SLURM_CPUS_ON_NODE
fi

# ---- set-up ends ------------------------------------------------------------

# ---- BODY -------------------------------------------------------------------
verbose "Create forces directory and extracting forces from $chkfile"
[[ -d forces ]] || fail "A directory called 'forces' have to exist to run
  execute the computation of forces"

compute_forces "$chkfile"
name=${chkfile//${pattern}/forces}
verbose "Moving result to forces/${name%.*}.*"
for fil in ${name%.chk}.*
do
  mv $fil "forces/${name%.*}.${fil##*.}" || fail "moving results to forces
    directory."
done

finish
