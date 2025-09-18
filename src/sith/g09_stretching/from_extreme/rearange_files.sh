#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool rearange the files called *<n>*.xyz changing the number <n>
for the correspondending in increasing order. it must be executed in the directory
where the files to organize are.

  -h  prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- general setup ---------------------------------------------------------
verbose=''
while getopts 'vh' flag; do
  case "${flag}" in
    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

# ---- BODY -------------------------------------------------------------------
source "$(sith basics -path)" RearangeFiles $verbose

n=0
for xyz_file in $(ls *.xyz | sort)
do
  num=$(echo "$xyz_file" | grep -oE '[0-9]+')
  nn=$(printf "%03d" $n )
  n=$(( n + 1 ))
  echo -n "($num/$nn)"
  if [ "$nn" -ne "$num" ]
  then
    for j in ${xyz_file%.xyz}*
    do
      rename=${j//$num/$nn}
      mv $j $rename || fail "moving file $j to $rename"
    done
  fi
done
echo
