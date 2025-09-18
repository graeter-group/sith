#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
Creates the com files from the xyz structures extracted from a gaussian log file
and submit the corresponding jobs to compute the forces.

  -l  <log_file> optimization gaussian logfile.
  -n  <name> standard name. Usually pep name.

  -h  prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
verbose=''
while getopts 'l:n:vh' flag;
do
  case "${flag}" in
    l) logfile=${OPTARG} ;;
    n) name=${OPTARG} ;;
    
    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source $(sith basics -path) "AfterOpt" $verbose

load_modules
# ---- BODY -------------------------------------------------------------------

# ==== Reduce number of structures with reduced changes of DOFs
verbose "Create continuous structutes path."
# The output are the xyz files without peak energies, output name-conopt<n>.xyz
sith info_from_opt $logfile ../${name}-stretched00.pdb ${name}-conopt \
  || fail "extracting xyz files from log file from $logfile"
# Extract the dofs from the created xyzs. out; <name>-conopt-dofs.dat
sith extr_dofs -f ${name}-conopt > /dev/null || \
  fail "extracting dofs from xyzs"
# reduce irrelevant changes, store the new subset in a dir called subset
sith reduce_structs "." ${name}-conopt > /dev/null || \
  fail "reducing structures"

# ==== Create com gaussian files
# Create com file template
sith opt_from_xyzs -d . -n ${name}-conopt000 -p ../${name}-stretched00.pdb \
  > /dev/null|| fail "creating com files using forces_from_xyz"
# clean files: only leaves the template
mv ${name}-conopt000.com template.com
echo "" >> template.com
rm *conopt*

# import xyz files of the subset
mv subset/* .
rm -r subset

# Create .com files
verbose "Create com files."
str_index=0
for file in ${name}-conopt*.dat
do
  struct_name=${file%.dat}
  echo $struct_name
  sith find_blocks -f template.com -e "Variables:" -o tmp > /dev/null
  mv tmp_001.out  $struct_name.com
  echo "     Variables:" >> $struct_name.com
  cat $file >> $struct_name.com
  echo "" >> $struct_name.com
  sith find_blocks -s "\^\$" -e "\^\$" -f template.com -o tmp > /dev/null
  cat tmp_001.out >> $struct_name.com
  sed -i "/chk=/c %chk=$struct_name" $struct_name.com

  if [[ "$(whoami)" == "hits_"* ]]
  then
    single_part="--partition=cpu-single"
  else
    single_part=""
  fi

  # This one can sbatch each job
  sbatch -J ${name}${str_index}cono $single_part \
    $( sith opt_and_forces -path ) -f $struct_name -c -v
  str_index=$(( 10#$str_index + 1 ))
done

rm tmp_001.out
rm template.com
rm *.dat
rm *.pdb

finish
