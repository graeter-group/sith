#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Submit opt_and_forces for each one of the xyz files given as positional
arguments.
"
exit 0
}

# ---- set up -----------------------------------------------------------------
verbose=''
while getopts 'vh' flag;
do
  case "${flag}" in
    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done
source "$(sith basics -path)" PrepareAndSubmit

# ---- BODY -------------------------------------------------------------------
for pep in $@
do
  verbose ${pep%*/}
  cd $pep || fail "$pep does not exist"
  index1=$( grep ACE "${pep%*/}-stretched00.pdb" | grep CH3 | \
            awk '{print $2}' )
  index2=$( grep NME "${pep%*/}-stretched00.pdb" | grep CH3 | \
            awk '{print $2}' )

  # check that the indexes were read properly:
  [[ "$index1" -eq 0 && "$index2" -eq 0 ]] && fail "Not recognized indexes"
  [[ "$index1" -eq 1 && "$index2" -eq 1 ]] && fail "Not recognized indexes"
  echo "$pep freezing $index1 and $index2"

  cd forces
  for file in *.xyz
  do
    sith change_distance \
            $file ${file%.xyz}-opt \
            "nofile" 0 0 "scale_distance" \
            || fail "Preparating gaussian input"
    sed -i '$d' ${file%.xyz}-opt.com
    echo "$index1 $index2 F" >> ${file%.xyz}-opt.com
    sed -i "1a %NProcShared=8" "${file%.xyz}-opt.com"
    sed -i "/#P/a opt(modredun,calcfc)" "${file%.xyz}-opt.com"

    if [[ "$(whoami)" == "hits_"* ]]
    then
      single_part="--partition=cpu-single"
    else
      single_part=""
    fi

    sbatch --job-name="${file:0:6}_opt" $single_part \
           --output="${file:0:6}_opt.o" \
           --error="${file:0:6}_opt.e" \
           $(sith opt_and_forces -path) -f ${file%.xyz}-opt -c || \
      fail "submitting opt_and_forces"
  done
  cd ../..
done

finish
