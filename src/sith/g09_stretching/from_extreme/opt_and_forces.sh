#!/bin/bash

#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
This code submits a gaussian job (typically an optimization) and uses the
output to compute the forces.

  -c  Use this flag to run in a cluster. When -p is not defined, and you run in
      a slurm job manager system, the number of processors is equal to the
      number of cores asked in the submission of the job.
  -f  <com file> name of the gaussian input file without extension (.com). The
      output is has the same name, but .log extension.
  -F  use this flag fo AVIOD force calculation after optimization.
  -p  <processors=1> number of processors per gaussian job. See description of
      flag -c.
  -S  <job_options=''> options for submitting a new job. This flag only makes
      sense in slurm cluster. Please, do not include a name and add the options
      as in the next example: \"--partition=cpu --nice\".

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ---- set up -----------------------------------------------------------------
cluster='false'
force_calc='true'
n_processors=''
job_options=''
verbose=''
while getopts 'cf:F:p:S:vh' flag; do
  case "${flag}" in
    c) cluster='true' ;;
    f) file=${OPTARG} ;;
    F) force_calc='false' ;;
    p) n_processors=${OPTARG} ;;
    S) job_options=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" OptAndForces $verbose

verbose -t "JOB information"
verbose -t "==============="
verbose -t " \* Date:"
verbose -t $(date)
verbose -t " \* Command:"
verbose -t "$0" "$@"

c_flag=''
if $cluster
then
  load_modules
  c_flag="-c"
  if [[ -z "$n_processors" ]]
  then
    if [[ ! -z "$SLURM_CPUS_ON_NODE" ]]
    then
      n_processors=$SLURM_CPUS_ON_NODE
    else
      n_processors=1
    fi
  fi
  job_options="$job_options -J='${SLURM_JOB_NAME}_F' -n $n_processors"
  job_options="sbatch $job_options"
fi

# ----- BODY ------------------------------------------------------------------
verbose "submit constrained optimization $file"
gaussian "$file.com" "$file.log"

grep -q "Normal termination of Gaussian" "$file.log" || \
  fail "Optimization did not work for $file"

if $force_calc
then
  verbose "Submit forces computation $file.chk"
  $job_options \
    $(sith find_forces -path) $c_flag -f $file.chk -p "conopt" \
                              $verbose || fail "submitting forces"
fi

finish
