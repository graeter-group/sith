#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
Creates a com file for each xyz file with a pattern that it finds in a given
directory. The name of each com file is the same than the xyz, but with different
extension.

  -d  directory where the xyz files are.
  -n  pattern in the name of the desired xyz file.
  -p  pdb of reference.

  -h  prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- general setup ---------------------------------------------------------
directory="./"
pattern=""
verbose=''
while getopts 'd:n:p:cvh' flag; do
  case "${flag}" in
    d) directory=${OPTARG} ;;
    n) name_pattern=${OPTARG} ;;
    p) pdb_ref=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

# ---- BODY -------------------------------------------------------------------
source "$(sith basics -path)" OptFromXYZ $verbose
original_forc2xyz_dir=$(pwd)

index1=$( grep ACE $pdb_ref | grep CH3 | awk '{print $2}' )
index2=$( grep NME $pdb_ref | grep CH3 | awk '{print $2}' )

cd $directory
for file in *$name_pattern*.xyz
do
  echo ${file%.xyz}
  tail -n +3 $file > tmp.xyz
  newzmat -ixyz -ozmat -rebuildzmat tmp ${file%.xyz} > /dev/null || \
    fail "executing newzmat"
  sed -i "s/-- No Title Specified --/Computation of forces/g" ${file%.xyz}.com
  sed -i "s/\# HF\/6-31G\* Test/%chk=$name_opt\n%NProcShared=8\n#P bmk\/6-31+g opt(modredun)/g" ${file%*.xyz}.com
  echo -e "$index1 $index2 F" >>  ${file%.xyz}.com
done
rm tmp.xyz

cd $original_forc2xyz_dir
