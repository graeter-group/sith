#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e


# ----- definition of functions -----------------------------------------------
print_help() {
echo "
This tool creates all the stretched structures for a molecule and computes the
needed quantities for sith. You can use this code to submit a Job in cluster or
to execute it locally. Consider the next options:

  -b  <number of breakages=1> The simulation will run until getting this number
      of ruptures in the bonds.
  -c  run in cluster (see the documentation of the installation of sith
      -execute 'sith doc' in your terminal-)
  -i  <index1,index2> indexes of the atoms to use for increasing the distance.
      If this flag is not used and 'molecule' (flag -m) is a pdb, indexes 1 and
      2 will correspond to the CH3 atoms in ACE and NME residues defined in the
      pdb if they exist.
  -l  <xc,base="bmk,6-31+g"> evel of DFT theory.
  -m  <molecule> definition of the molecule (xyz, pdb, ...).
  -M  <method=0> Index of stretching method. To see the options, use
      'sith change_distance -h' to see the order.
  -n  <n_processors=1> number of processors per gaussian job.
  -r  restart. In this case, run from the directory of the pre-created
      stretched molecule.
  -s  <size[A]=0.2> of the step that increases the distances.
  -S  <job_options=''> options for submitting a new job. This flag only makes
      sense in slurm cluster. Please, do not include a name and add the options
      as in the next example: \"--partition=cpu --nice\".

  -v  verbose.
  -h  prints this message.
"
exit 0
}

resubmit () {
  sleep 23h 58m

  verbose "Resubmitting the job for another 24h";
  sbatch "$job_options" -J "$SLURM_JOB_NAME" "$( sith workflow -path)" \
    -b "$breakages" \
    $cluster \
    -i "$indexes" \
    -l "$level" \
    -m "$molecule" \
    -M "$method" \
    -n "$n_processors" \
    -r \
    -s "$size" \
    -S "$job_options" \
    $verbose || fail "Resubmission of the job failed"
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
breakages=1
cluster=''
indexes=''
level="bmk,6-31+g"
method=0
n_processors=''
restart=''
size=0.2
job_options=""

verbose=''
while getopts 'b:ci:l:m:M:n:rs:S:vh' flag;
do
  case "${flag}" in
    b) breakages=${OPTARG} ;;
    c) cluster='-c' ;;
    i) indexes=${OPTARG} ;;
    l) level=${OPTARG} ;;
    m) molecule=${OPTARG} ;;
    M) method=${OPTARG} ;;
    n) n_processors=${OPTARG} ;;
    r) restart='-r' ;;
    s) size=${OPTARG} ;;
    S) job_options=${OPTARG} ;;

    v)  verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" WORKFLOW $verbose

# starting information
verbose -t "JOB information"
verbose -t "==============="
verbose -t " \* Date:"
verbose -t "$(date)"
verbose -t " \* Command:"
verbose -t "$0" "$@"

resubmit &
# ---- Set up -------------------------------------------------------------------

# load modules
if [[ "$cluster" == "-c" ]]
then
  load_modules "$molecule" "$method" "$breakages" "$size"
  if [[ -z "$n_processors" ]] 
  then
    if [[ -n "$SLURM_CPUS_ON_NODE" ]]
    then
      n_processors=$SLURM_CPUS_ON_NODE
    else
      n_processors=1
    fi
  fi
fi

ase -h &> /dev/null || fail "This code needs ASE."
command -V gaussian &> /dev/null || fail "Remeber to define the function
  gaussian (check the documentation of the installation -sith doc-)."
# gmx -h &> /dev/null || fail "This code needs gmx."
sith -h &> /dev/null || fail "This code needs sith."

# ---- BODY -------------------------------------------------------------------
verbose "execute stretching"

sith stretching -b "$breakages" $cluster -e "$method" -i "'$indexes'" \
                -l "$level" \
                -m "$molecule" -p "$n_processors" $restart -s "$size" \
                $verbose || fail "Stretching of $molecule failed"

# TODO: add Classical energies

# compute forces
verbose "submitting comptutation of forces.";

create_bck forces
mkdir -p forces
mkdir -p bck
mv ./*-bck*.* bck

Fcounter=0
for i in "${molecule%.*}"-stretched*.chk
do
  verbose "Forces of $i"
  sbatch $job_options -J "F${molecule%.*}" \
    "$(sith find_forces -path)" "$cluster" -f "$i" -p stretched $verbose || \
    fail "Submitting forces calculation"
  Fcounter=$(( Fcounter + 1 ))
done

cd ../
sbatch $job_options -J "FE$molecule" \
  $(sith workflow_from_extreme -path) -a "$molecule" $cluster -m "$molecule" \
    -t "$molecule/$molecule-stretched00.pdb" $verbose -S "$job_options"


#sbatch $single_part -J ${molecule}_WAR $(sith workflow_from_extreme -path) -c -p "." -v

#verbose "running grappa and amber";
#final_force=$(sith F_max_stretch ../ $molecule)
#sbatch $single_part -J ${molecule}_ff $(sith pulling_with_ff -path) -F $final_force -f $molecule-stretched00.pdb 


finish
