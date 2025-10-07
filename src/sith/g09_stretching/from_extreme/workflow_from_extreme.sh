#!/bin/bash

#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e

# ----- definition of functions starts ----------------------------------------
print_help() {
cat << EOF
This tool creates the files to do the sith analysis by optimizing a molecule,
then takes the intermedia steps and prepares a constrained optimization of each
one of them with :bashscript:`sith.g09_stretching.from_extreme.continuous_path`
after creating intermedias states that guarantee continuity in the degrees of
freedom. This script submits in parallel each one of those optimizations with
:bashscript:`sith.g09_stretching.from_extreme.opt_and_forces` (-S here refers
to these subjobs).

  -a  <alias> new name of the xyz file and subsequent files in directory
      'from_extreme'.
  -c  Use this flag to run in a cluster. When -p is not defined, and you run in
      a slurm job manager system, the number of processors is equal to the
      number of cores asked in the submission of the job.
  -i  <index1,index2> indexes of the atoms used for constraining the distance
      of the intermedia structures when optimizing. If this flag is not used
      but a pdb file is given (-t), indexes 1 and 2 will correspond to the CH3
      atoms in ACE and NME residues defined in the pdb if they exist.
  -l  <xc,base="bmk,6-31+g"> level of DFT theory.
  -m  <molecule> directory or coordinates file of configuration to be relaxed.
      For example, \"./AAA/\" a trialanine configuration ('last' after
      organizing all xyz files alphabetically all AAA*.xyz files in ./AAA/).
  -p  <processors=1> number of processors per gaussian job. See description of
      flag -c.
  -r  Use it to restart, in which case no directory will be created and the new
      <molecule>-optext.log is assumed to exist in the working directory.
  -S  <job_options=''> options for submitting a new job. This flag only makes
      sense in slurm cluster. Please, do not include a name (-J), nor the
      number of cores (-n, use -p for this). The input should be as in the next
      example: \"--partition=cpu --nice\".
  -t  <template.pdb> template pdb file to define indexes if -i is not used.

  -v  verbose.
  -h  prints this message.

Note
----

  The outputs are stored in a directory called 'from_extreme'/

  This tool requires gaussian and sklearn.
EOF
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
cluster='false'
indexes=''
level="bmk,6-31+g"
molecule=''
n_processors=''
restart='false'
restart='false'
job_options=''

verbose=''
while getopts 'a:ci:l:m:p:rS:t:vh' flag;
do
  case "${flag}" in
    a) alias=${OPTARG} ;;
    c) cluster='true' ;;
    i) indexes=${OPTARG} ;;
    l) level=${OPTARG} ;;
    m) molecule=${OPTARG} ;;
    p) n_processors=${OPTARG} ;;
    r) restart='true' ;;
    S) job_options=${OPTARG} ;;
    t) template=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" "WFromE" $verbose

# starting information
verbose -t "JOB information"
verbose -t "==============="
verbose -t " \* Date:"
verbose -t $(date)
verbose -t " \* Command:"
verbose -t "$0" "$@"

c_flag=''
if $cluster
then
  load_modules # TODO: ADD the parameters to resubmit
  c_flag='-c'
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

xc_functional=$(echo $level | cut -d ',' -f 1)
basis_set=$(echo $level | cut -d ',' -f 2)

if [[ -z "$indexes" ]] && [[ -f "$template" ]]
then
  # reading indexes from pdb file
  index1=$( grep ACE "$template" | grep CH3 | awk '{print $2}' )
  index2=$( grep NME "$template" | grep CH3 | awk '{print $2}' )
  indexes="$index1,$index2"
  verbose -t "Indexes read from $template: $indexes"
fi

if [ -z $molecule ]
then 
  fail "This code needs a reference  molecule. Please, define it using the flag
        -m. For more info, use \"sith workflow_from_extreme -h\""
fi

# ---- BODY -------------------------------------------------------------------
# In case pep is a directory, it searches the last xyz in this dir.
if [ -d $molecule ]
then
    cd $molecule
    molecule=${molecule##*/}
    mapfile -t previous < <( find . -maxdepth 1 -type f\
                                    -name "$molecule*.xyz" \
                                    -not -name "*bck*" | sort )
    molecule=${previous[-1]}
    verbose -t "'$molecule' found to be the last structure."
fi

# ==== Optimization from extreme
xyz=${molecule##*/}
name=${xyz%.*}

[ -z "$alias" ] && alias=$name || name=$alias

verbose -t "Create $name-optext.com file."
if [[ "$restart" == "false" ]]
then
  if [ ! -f $molecule ]
  then
      fail "$molecule does not exist."
  fi
  # create from_extreme directory
  mkdir -p from_extreme
  cp $molecule from_extreme/$alias.xyz
  cd from_extreme
else
  tmp_path=$(pwd)
  [[ "${tmp_path##*/}" == "from_extreme" ]] && cd ../
  cd from_extreme || fail "directory 'from_extreme' is not [and does not exist
      in] the current directory."
  if  ! grep -q "Normal termination" "$name-optext.log" "$name-optext-bck*.log"
  then  
    sith log2xyz "$name-optext.log" > /dev/null || fail "extracting coordinates from process
      from $name-optext.log"
    create_bck $xyz
    cp $name-optext.xyz $xyz
  fi
fi

# run gaussian
verbose -t "Running optimization from $name-optex.com."
if  ! grep -q "Normal termination" "$name-optext.log" "$name-optext-bck*.log"
then
  # creates gaussian input that optimizes the structure
  sith change_distance "$name.xyz" "$name-optext" no_frozen_dofs 0 0 \
    "scale_distance" --xc "'$xc_functional'" --basis "'$basis_set'" \
    > /dev/null || fail "Preparating the input of gaussian"
  sed -i "1a %NProcShared=$n_processors" "$name-optext.com"
  sed -i "/#P/a opt(modredun,calcfc)" "$name-optext.com"
  sed -i "1a %mem=60000MB" "$name-optext.com"

  gaussian "$name-optext.com" "$name-optext.log" || \
    { if [ "$(grep -c "Atoms too close." \
           "$name-optext.com")" \
           -eq 1 ]; then fail "Atoms too close for ${nameiplusone}" ; \
      fi ; } || fail "running gaussian optimization"
  # check convergence from output
  output=$(grep -i optimized "$name-optext.log" | \
           grep -c -i Non )

  [ "$output" -ne 0 ] && fail "Optimization didn't converged"
fi

if [[ $(ls *optext*.log | wc -l ) -gt 1 ]];
then
  verbose -t "Concatenate all the generated logfiles: $name-optext*.log"
  create_bck $name-optext.log
  for bck_logfile in $(ls $name-optext*.log | sort )
  do
    cat $bck_logfile >> $name-optext.log
    rm $bck_logfile
  done
fi

# ==== Reduce number of structures with reduced changes of DOFs
verbose -t  "Create continuous structutes path with 'sith info_from_opt'."
# The output are the xyz files without peak energies, output
# <name>-conopt<n>.xyz
sith info_from_opt $name-optext.log ${name}-conopt > /dev/null || \
  fail "extracting xyz files from log file from $name-optext.log"

verbose "Starting 'sith continuous_path' after having all
  ${name}-conopt<n>.xyz files"

$(sith continuous_path -path) $c_flag -i "$indexes" -l "$level" -n "$name" \
                              -p "$n_processors" -P "conopt" -S "$job_options"
                              $verbose || \
  fail "submiting continuous path"

finish "continuous path of $name finished."
