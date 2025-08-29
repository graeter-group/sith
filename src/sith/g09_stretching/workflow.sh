#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 9
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive


# ----- definition of functions -----------------------------------------------
print_help() {
echo "
This tool creates all the stretched structures for a peptide and computes the
needed quantities for sith. You can use this code to submit a Job in cascade or
to execute it locally. Consider the next options:

  -b  <number of breakages=1> The simulation will run until getting this number
      of ruptures in the bonds.
  -c  run in cascade. (modules are loaded)
  -d  <reference document=00-aminos.txt> file containing the existing peptides
      to avoid repetition. Mainly used for generation of random peptides.
  -e  <endo> or <exo> states for initial state of proline. Default 'random'.
  -m  <method=0> Index ofstretching method. To see the options, use
      'sith change_distance -h' to see the order.
  -n  <options> Pepgen options (use \\\" for this)
  -p  <peptide> Chains of aminoacids to be evaluated. For example, \"AAA\"
      would analyse a trialanine peptide.
  -R  random pepeptide. Give the number of amino acids with this argument.
  -r  restart. In this case, run from the directory of the pre-created
      peptide.
  -s  <size[A]=0.2> of the step that increases the distances.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

resubmit () {
  sleep 23h 58m ; \
  if [[ "$(whoami)" == "hits_"* ]]
  then
    single_part="--partition=cpu-single"
  else
    single_part=""
  fi
  sbatch $single_part -J $SLURM_JOB_NAME "$( sith workflow -path)" -p "$1" -c -r -m "$2" -b "$3" -s "$4" ; \
  echo "new JOB submitted"
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
breakages=1
cascade='false'
endoexo='random'
method=0
pep=''
restart=''
size=0.2
ref_doc='00-aminos.txt'

verbose='false'
while getopts 'd:b:ce:m:n:p:rR:s:vh' flag;
do
  case "${flag}" in
    b) breakages=${OPTARG} ;;
    c) cascade='true' ;;
    d) ref_doc=${OPTARG} ;;
    e) endoexo=${OPTARG} ;;
    m) method=${OPTARG} ;;
    n) pep_options=${OPTARG} ;;
    p) pep=${OPTARG} ;;
    r) restart='-r' ;;
    R) random=${OPTARG} ;;
    s) size=${OPTARG} ;;

    v)  verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" WORKFLOW $verbose

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

# ---- Set up -------------------------------------------------------------------
# random peptide
if [ ! "${#random}" -eq 0 ]
then
  [ -f "$ref_doc" ] || fail "Non-recognized $ref_doc, check flag -d"
  pep=$( sith gen_randpep "$random" ) || fail "Creating random peptide"
  while awk '!/^#/ {print $1}' "$ref_doc" | grep -q "$pep"
  do
    pep=$( sith gen_randpep "$random" ) || fail "Creating random peptide"
  done
  echo "$pep" "   R" >> "$ref_doc"
  warning "The code created the random peptide $pep, the workflow will run
    with this peptide even if you also passed -p argument."
fi

# debug peptides
if [ "${#pep}" -eq 0 ]
then 
  fail "This code needs one peptide. Please, define it using the flag -p or
        -R. For more info, use \"sith workflow -h\""
fi

# load modules
if $cascade
then
  load_modules "$pep" "$method" "$breakages" "$size"
fi

ase -h &> /dev/null || fail "This code needs ASE"
command -V g09 &> /dev/null || fail "This code needs gaussian"
gmx -h &> /dev/null || fail "This code needs gmx"
sith -h &> /dev/null || fail "This code needs sith"
perl -E "say '+' x 80"

# ---- BODY -------------------------------------------------------------------
# ---- firstly, backup previous directories with the same name
if [[ "$restart" != "-r" ]]
then
  # check pepgen
  pepgen -h &> /dev/null || fail "This code needs pepgen"

  # create back up
  create_bck "$pep"

  # Creation of the peptide directory and moving inside.
  mkdir "$pep"
  cd "$pep" || fail "directory $pep does not exist"
  verbose "generating peptide"
  # Creation of peptide
  # shellcheck disable=SC2086
  pepgen "$pep" tmp -r -s flat $pep_options || fail "Creating peptide $pep"
  mv tmp/pep.pdb "./$pep-stretched00.pdb"
  sith classical_minimization -f "./$pep-stretched00.pdb" \
                                 -o "./$pep-stretched00.pdb"
  verbose "define proline state"
  sith proline_mod -f "$pep-stretched00.pdb" -s "$endoexo" || \
    fail "Proline estates configuration"
  mv "$pep-stretched00modpro.pdb" "$pep-stretched00.pdb" 
  verbose "protonate/deprotonate"
  sith protonate "./$pep-stretched00.pdb" "./$pep-stretched00.pdb" || \
    fail "protonizing"
  rm -r tmp
else
  # moving to the peptide directory
  cd "$pep" || fail "directory $pep does not exist"
  warning "$pep restarted"
fi

# construct the stretched configurations
sith stretching -b "$breakages" -p "$pep" "$restart" -m "$method" \
                   -s "$size" || fail "Stretching of $pep failed"

# Compute classical energies
verbose "computing classical energies."
sith classical_energies

# compute forces
verbose "submitting comptutation of forces.";

if [[ "$(whoami)" == "hits_"* ]]
then
  single_part="--partition=cpu-single"
else
  single_part=""
fi

sbatch $single_part -J ${pep}_forces "$( sith find_forces -path )" -c  -p $pep &&
  echo "computation of forces submitted"

sbatch $single_part -J ${pep}_WAR $(sith workflow_from_extreme -path) -c -p "." -v

verbose "running grappa and amber";
final_force=$(sith F_max_stretch ../ $pep)
sbatch $single_part -J ${pep}_ff $(sith pulling_with_ff -path) -F $final_force -f $pep-stretched00.pdb 


finish "$pep finished"
