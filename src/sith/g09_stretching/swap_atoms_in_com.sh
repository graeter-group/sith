#!/bin/bash

# ---- functions --------------------------------------------------------------
print_help() {
echo "
Permutes the incides of two atoms in com file where the molecule is defined as
a z-matrix.

  -a  <atom1> index of one of the atoms to be permuted. 
  -b  <atom2> index of the other atom to be permuted.
  -f  <file.com> input gaussian file with the z-matrix already defined.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

permutate_pattern() {
    pattern1="$1"
    pattern2="$2"
    fil="$3"
    # permutate index in connectivity
    sed -i "s/$pattern1/tmp_pattern1/g" $fil
    sed -i "s/$pattern2/$pattern1/g" $fil
    sed -i "s/tmp_pattern1/$pattern2/g" $fil
}

# ----- set up starts ---------------------------------------------------------
# General variables
verbose=''
while getopts 'a:b:f:vh' flag;
do
  case "${flag}" in
    a) atom1=${OPTARG} ;;
    b) atom2=${OPTARG} ;;
    f) file=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

# ---- BODY -------------------------------------------------------------------
source "$(sith basics -path)" SwapAtomsInCom $verbose

# change lines
def_atom1=$( grep ",R$atom1," $file )
def_atom2=$( grep ",R$atom2," $file )

# permutate definition
permutate_pattern "$def_atom1" "$def_atom2" $file

# permutate index in connectivity
permutate_pattern ",$atom1," ",$atom2," $file
permutate_pattern "\ $atom1\ " "\ $atom2\ " $file

# permutate index in R,A,D
permutate_pattern "R$atom1" "R$atom2" $file
permutate_pattern "A$atom1" "A$atom2" $file
permutate_pattern "D$atom1" "D$atom2" $file
