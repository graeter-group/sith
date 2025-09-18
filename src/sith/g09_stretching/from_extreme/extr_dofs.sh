#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool extract the dofs of a set of xyz files. The ouput are files called
<xyz file with *pattern*>-dofs.dat, where pattern comes from -f flag.

  -f  <xyz files pattern>. The code will look for *pattern*.xyz

  -h  prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- general setup ---------------------------------------------------------
verbose=''
while getopts 'f:p:vh' flag; do
  case "${flag}" in
    f) xyzs=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" Extract_DOFs
# ---- BODY -------------------------------------------------------------------
# Reduce, optimize and then try to find intermedias.
for xyzfile in *"${xyzs}"*.xyz
do
  tail -n +3 $xyzfile > tmp.xyz
  newzmat -ixyz -ozmat -rebuildzmat -bmodel \
    tmp.xyz ${xyzfile%.xyz}-forces.com || fail "z-matrix"
  n=$(grep -n "Variables:" ${xyzfile%.xyz}-forces.com | awk '{print $1}')
  n=${n%:}

  # check that the new structure has the same dofs in tha z-matrix
  if [ ! -f mat_inf.dat ]
  then
    head -n $n ${xyzfile%.xyz}-forces.com > mat_inf.dat
  else
    head -n $n ${xyzfile%.xyz}-forces.com > tmp2.dat
    diff -q mat_inf.dat tmp2.dat || fail "different matrix definition in
      $xyzfile"
  fi

  # save dofs
  end=$(grep -n "^ D" ${xyzfile%.xyz}-forces.com | tail -n 1)
  end=${end%:*}
  head -n $end ${xyzfile%.xyz}-forces.com | \
    tail -n +$n > ${xyzfile%.xyz}-dofs.dat
  rm ${xyzfile%.xyz}-forces.com
done

rm tmp.xyz
rm tmp2.dat
rm mat_inf.dat

finish
