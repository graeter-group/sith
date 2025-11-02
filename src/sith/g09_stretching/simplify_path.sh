#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Takes all log files with certain pattern, extracts the last configuration,
extract the dofs and removes the steps corresponding to oscillations or the
system stacked in a configuration.

  -p  <pattern> pattern of the log files to be considered

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
pattern=""
verbose=''
cluster='false'
while getopts 'cp:vh' flag;
do
  case "${flag}" in
    c) cluster='true' ;;
    p) pattern=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith simplify_path -h" >&2 ; exit 1 ;;
  esac
done

# shellcheck disable=SC1090
source "$(sith basics -path)" SimplifyPath "$verbose"

# starting information
verbose -t "JOB information"
verbose -t "==============="
verbose -t " Command:" "$0" "$@"

# load modules
if $cluster
then
  load_modules
fi

# ---- BODY -------------------------------------------------------------------
mkdir reduced_path
for i in *$pattern*.log
do
  cp $i reduced_path
done

cd reduced_path
for i in *.log
do
  verbose -t "$i to xyz"
  sith log2xyz $i  > /dev/null || fail "extracting xyz from $i"
done

sith extr_dofs -f "$pattern" $verbose

sith reduce_structs './' 'test' --print_reduced True > selected.dat \
  || fail "reducing structures"

mapfile -t selected < <(sith find_blocks -f selected.dat -e '\\[\\[' -o terminal\
                        | grep -oE '[0-9]+')

rm -r *

new_count=0
for i in "${selected[@]}"
do
  file2copy=$(ls ../*$i*.fchk)
  index=$(printf "%03d" "$new_count")
  verbose -t "Copying $file2copy"
  file_name=${file2copy##*/}
  cp $file2copy ${file_name//$i/$index}
  ((new_count++))
done

# ---- END --------------------------------------------------------------------
finish "message to finish"

