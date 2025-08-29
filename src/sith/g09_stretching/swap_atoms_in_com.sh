#!/bin/bash
print_help() {
echo "
# TODO: repair

Use this template to create your scripts with a standard structure

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
def_var="inse here your default"
cascade='false'
verbose='false'
while getopts 'a:b:f:vh' flag;
do
  case "${flag}" in
    a) atom1=${OPTARG} ;;
    b) atom2=${OPTARG} ;;
    f) file=${OPTARG} ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

permutate_pattern() {
    pattern1="$1"
    pattern2="$2"
    fil="$3"
    # permutate index in connectivity
    sed -i "s/$pattern1/tmp_pattern1/g" $fil
    sed -i "s/$pattern2/$pattern1/g" $fil
    sed -i "s/tmp_pattern1/$pattern2/g" $fil
}

# change lines
def_atom1=$( grep ",R$atom1," $file )
def_atom2=$( grep ",R$atom2," $file )

# permutate definition
permutate_pattern "$def_atom1" "$def_atom2" $file

# permutate index in connectivity
permutate_pattern ",$atom1," ",$atom2," $file
permutate_pattern "\ $atom1\ " "\ $atom2\ " $file

# permutate indez in R,A,D
permutate_pattern "R$atom1" "R$atom2" $file
permutate_pattern "A$atom1" "A$atom2" $file
permutate_pattern "D$atom1" "D$atom2" $file
