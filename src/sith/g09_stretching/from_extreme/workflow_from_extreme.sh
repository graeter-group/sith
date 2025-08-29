#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 9
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive


# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool creates the files to do the sith analysis by optimizing a molecule
that was just about to get a first rupture, then takes the intermedia steps and
find the internal forces. Consider the next options:

  -c  run in cascade. (modules are loaded)
  -p  <peptide>. directory or xyzfile of last conf. Chains of aminoacids to
      be evaluated. For example, \"./AAA/\" would optimize the last
      stretched a trialanine peptide (where last means after organizing
      alphabetically).
  -r  Use if the job corresponds to a restart, in which case, no directory will
      be created and the new <peptide>-optext.com is assumed to exist.

  -v  verbose.
  -h  prints this message.

Note: it is assumed that the file of the last configuration is named as:
<amino acids-code>-<description><number of stretching>.xyz
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
cascade='false'
ref=''
verbose='false'
restart='false'
while getopts 'cp:rvh' flag;
do
  case "${flag}" in
    c) cascade='true' ;;
    p) ref=${OPTARG} ;;
    r) restart='true' ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" WF_FROM_EXTREME $verbose

if $cascade
then
    load_modules # ADD the parameters to resubmit
fi

if [ "${#ref}" -eq 0 ]
then 
  fail "This code needs one reference. Please, define it using the flag -p.
    For more info, use \"sith workflow_from_extreme -h\""
fi

# In case pep is a directory, it searches the last xyz in this dir.
if [ -d $ref ]
then
    cd $ref
    ref=${ref##*/}
    mapfile -t previous < <( find . -maxdepth 1 -type f -name "*.xyz" \
                                    -not -name "*bck*" | sort )
    ref=${previous[-1]}
fi

if [ ! -f $ref ]
then
    fail "$ref does not exist"
fi

# ---- BODY -------------------------------------------------------------------
# ==== Optimization from extreme
xyz=${ref##*/}
name=${xyz%-*}

verbose "The first g09 process is an optimization starting from $ref"

if [[ "$restart" == "false" ]]
then
  # create from_extreme directory
  mkdir from_extreme
  cp $xyz from_extreme
  cp *00.pdb from_extreme
  cd from_extreme
  
  # creates gaussian input that optimizes the structure
  sith change_distance "$xyz" "$name-optext" no_frozen_dofs 0 0 \
    "scale_distance" || fail "Preparating the input of gaussian"
  rm "$xyz"
  sed -i "1a %NProcShared=8" "$name-optext.com"
  sed -i "/#P/a opt(modredun,calcfc)" "$name-optext.com"
fi

# run gaussian
verbose "Running optmization of stretching ${nameiplusone}"

if ! grep -q "Normal termination" "$name-optext.log"
then
  g09 "$name-optext.com" "$name-optext.log" || \
    { if [ "$(grep -c "Atoms too close." \
           "$name-optext.com")" \
           -eq 1 ]; then fail "Atoms too close for ${nameiplusone}" ; \
      fi ; } || fail "running gaussian optimization"

  # Restart in case of i/0 problems
  if $(grep -q "NtrErr Called from FileIO." "$name-optext.log")
  then
    verbose "resubmit because of FileIO error"
    $(sith resubmit_failed -path) \
      -e "$(sith workflow_from_extreme -path) -r -p ${name}-optext.log -c " \
      -c ${name}-optext.com -l ${name}-optext.log \
      -j ${SLURM_JOB_NAME} || \
      fail "resubmitting $file after NtrErr Called from FileIO"
    fail "$file failed, it was submitted again"
  fi
  # check convergence from output
  output=$(grep -i optimized "$name-optext.log" | \
           grep -c -i Non )

  [ "$output" -ne 0 ] && fail "Optimization didn't converged"
fi

# concatenate all the generated logfiles
if [[ $(ls *optext*.log | wc -l ) -gt 1 ]];
then
  # in case of an optimization was already made in a bck process
  grep -q "Normal termination of Gaussian" *optext-bck*.log && rm *-optext.log
  create_bck $name-optext.log
  for bck_logfile in $(ls $name-optext*.log | sort )
  do
    cat $bck_logfile >> $name-optext.log
  done
fi

verbose "starting after optimization"
$(sith after_optimization -path) -l "$name-optext.log" -n $name -v

finish "$name finish"
