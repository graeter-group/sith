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
This tool creates a peptide as a change of amino acids. It adds ACE and NME
capping groups in the ends. The output is a pdb file with the name of the
peptide saved in a directory with the same name.

  -d  <reference document=00-aminos.txt> file containing the existing peptides
      to avoid repetition. Mainly used for generation of random peptides.
  -e  <endo> or <exo> states for initial state of proline. Default 'random'.
  -n  <options> Pepgen options (use \\\" for this)
  -p  <peptide> Chains of aminoacids to be evaluated. For example, \"AAA\"
      would analyse a trialanine peptide.
  -R  random pepeptide. Give the number of amino acids with this argument.

  -v  verbose.
  -h  prints this message.

Note
----

  This tools requires pepgen already installed.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
endoexo='random'
pep=''
ref_doc='00-aminos.txt'
random=''

verbose=''
while getopts 'd:e:n:p:R:vh' flag;
do
  case "${flag}" in
    d) ref_doc=${OPTARG} ;;
    e) endoexo=${OPTARG} ;;
    n) pep_options=${OPTARG} ;;
    p) pep=${OPTARG} ;;
    R) random=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: pkgdeveloper <function> -h" >&2 ; exit 1 ;;
  esac
done

source $(sith basics -path) CreateAminoAcid "$verbose"

# starting information
verbose -t "JOB information"
verbose -t "==============="
verbose -t " * Date:"
verbose -t $(date)
verbose -t " * Command:"
verbose -t "$0" "$@"

# --- set up ------------------------------------------------------------------
ase -h &> /dev/null || fail "This code needs ASE."
pepgen -h &> /dev/null || fail "This code needs pepgen"

# ---- BODY -------------------------------------------------------------------
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

# create back up
create_bck "$pep"

# Creation of the peptide directory and moving inside.
mkdir "$pep"
cd "$pep" || fail "directory $pep does not exist"
verbose "generating peptide"
# Creation of peptide
# shellcheck disable=SC2086
pepgen "$pep" tmp -r -s flat $pep_options|| fail "Creating peptide $pep"
mv tmp/pep.pdb "./$pep-stretched00.pdb"
# TODO: add again classical minimization
#sith classical_minimization -f "./$pep-stretched00.pdb" \
#                                -o "./$pep-stretched00.pdb"
verbose "define proline state"
sith proline_mod -f "$pep-stretched00.pdb" -s "$endoexo" > /dev/null || \
  fail "Proline estates configuration"
mv "$pep-stretched00modpro.pdb" "$pep-stretched00.pdb" 
verbose "protonate/deprotonate"
sith protonate "./$pep-stretched00.pdb" "./$pep-stretched00.pdb" > /dev/null || \
  fail "protonizing"
rm -r tmp

finish
