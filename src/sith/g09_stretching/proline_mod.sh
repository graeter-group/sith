#!/usr/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
Changes the state of the proline to endo, exo or random.

  -f  <path> pdb file.
  -o  <path> output pdb file.
  -l  <path> log file of the gromacs outputs. Default /dev/null
  -s  <state> proline state. So far, random, endo and exo are accepted.

  -v  verbose.
  -h  prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
proline_state='random'
outfile=''
pdbfile=''
outgromacs='/dev/null'
verbose='false'
while getopts 'f:o:l:s:vh' flag;
do
  case "${flag}" in
    f) pdbfile=${OPTARG} ;;
    s) proline_state=${OPTARG} ;;
    o) outfile=${OPTARG} ;;
    l) outgromacs=${OPTARG} ;;

    v)  verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" PROLINE_MODE $verbose

if [ "${#outfile}" -eq 0 ]
then
  outfile="${pdbfile%.*}modpro.pdb"
fi

source "$(sith basics -path)" PROLINE_MODE
# ---- BODY -------------------------------------------------------------------

# checking dependencies
[ "${#pdbfile}" -eq 0  ] && fail "To use proline modification, you have to
  provide the pdb file. use 'sith proline_mod -h' to see your options."

# changing proline states.
sith classical_minimization  -f "$pdbfile" -o "$outfile"  -l "$outgromacs" \
  || fail "minimization before proline definition of states"

verbose "define proline states"
sith proline_state "$outfile" "$proline_state" || fail "defining proline states"

sith classical_minimization -f "${outfile%.*}modpro.pdb" -l "$outgromacs" \
  || fail "minimization after proline definition of states"

mv "${outfile%.*}modpro.pdb" "$outfile"

finish
