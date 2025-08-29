#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
This code submit an optimization job and uses the output to compute the
forces.

  -f  name if the gaussian input file without extension (.com).
  -c  run in server.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ---- set up -----------------------------------------------------------------
c_flag=""
cascade='false'
while getopts 'f:cvh' flag; do
  case "${flag}" in
    f) file=${OPTARG} ;;
    c) cascade='true' ;;
    
    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done
source "$(sith basics -path)" opt_forces $verbose

verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

if $cascade
then
  load_modules
  c_flag="-c"
fi
# ----- BODY ------------------------------------------------------------------
verbose "submit constrained optimization $file"
grep -q "%mem" $file.com || sed -i "1a %mem=60000MB" $file.com
g09 "$file.com" "$file.log"

if $(grep -q "NtrErr Called from FileIO." "$file.log")
then
  $(sith resubmit_failed -path) \
          -e "$(sith opt_and_forces -path ) -c -v -f $file" \
          -c "$file.com" -l "$file.log" -j $SLURM_JOB_NAME \
          -f "$(head -n 1 ../frozen_dofs.dat)" -v || \
    fail "resubmitting $file after NtrErr Called from FileIO"
  fail "$file failed, it was submitted again"
fi

grep -q "Normal termination of Gaussian" "$file.log" || \
  fail "optimization did not work for $file"

if [[ "$(whoami)" == "hits_"* ]]
then
  single_part="--partition=cpu-single"
else
  single_part=""
fi

verbose "submit forces computation"
sbatch --job-name="${SLURM_JOB_NAME}_forces" $single_part \
       $(sith compute_forces -path) -f $file.chk -c -v || fail "submitting forces"

finish "optmimization"
