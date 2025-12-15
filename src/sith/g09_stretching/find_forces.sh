#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e


# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool computes the forces from a chk files that contains a given structure
and saves them in a directory called 'forces' that has to be previously
created.

  -c  run in cluster.
  -f  <file> chk file or com file to compute the forces from. In case of com
      file, it is assumed that it contains the right keywords to compute
      forces. In case of chk file, a com file is created with the right
      keywords to compute forces and replacing <pattern> with the word
      'forces'.
  -n  <n_processors=1> number of processors to be used in the gaussian job.
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
  newzmat -ichk -ozmat -rebuildzmat -bmodel "$1" $for_name.com > /dev/null || \
    {
      lnbck=$(search_last_bck ${1%.chk}) ; \
      newzmat -ichk -ozmat -rebuildzmat -bmodel "${1%.chk}-bck_$lnbck.chk" \
      $for_name.com > /dev/null || fail "Creating the matrix"
    }
  sed -i "1i %NProcShared=$n_processors" $for_name.com
  sed -i "1i %chk=$for_name" $for_name.com
  sed -i "s/opt(modredun,calcfc)/force/g" $for_name.com
  # remove unnecessary keywords, they are added sometimes to restart failed
  # processes of optimization:
  sed -i "s/Guess=Read Geom=Check//g" $for_name.com

  verbose "Executes gaussian computation of forces for $1"
  gaussian $for_name.com || fail "computing forces"
  sith extract_forces -f $for_name.log -c -v
}

# ----- definition of functions finishes --------------------------------------

# ---- set-up starts ----------------------------------------------------------
cluster='false'
directory='./'
pattern=''
verbose=''
n_processors=1
while getopts 'cf:n:p:vh' flag; do
  case "${flag}" in
    c) cluster='true' ;;
    f) file=${OPTARG} ;;
    n) n_processors=${OPTARG} ;;
    p) pattern=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" FindForces $verbose

verbose -t "JOB information"
verbose -t "==============="
verbose -t " \* Date:"
verbose -t $(date)
verbose -t " \* Command:"
verbose -t "$0" "$@"

if $cluster
then
  load_modules
  if [[ -z "$n_processors" ]]
  then
    if [[ ! -z "$SLURM_CPUS_ON_NODE" ]]
    then
      n_processors=$SLURM_CPUS_ON_NODE
    else
      n_processors=1
    fi
  fi
fi

# ---- set-up ends ------------------------------------------------------------

# ---- BODY -------------------------------------------------------------------

if [[ "${file##*.}" == "chk" ]]
then
  compute_forces "$file"
  name=${file//${pattern}/forces}
  verbose "Moving result to forces/${name%.*}.*"
  mkdir -p forces
  for fil in ${name%.chk}.*
  do
    mv $fil "forces/${name%.*}.${fil##*.}" || fail "moving results to forces
      directory."
  done
elif [[ "${file##*.}" == "com" ]]
then
  gaussian $file || fail "computing forces"
  for_name=${file//.com/}
  sith extract_forces -f $for_name.log -c -v
fi

finish
