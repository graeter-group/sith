#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
cat << EOF
Creates the com files from the a set of xyz structures, and submit the jobs
with :bashscript:`sith.g09_stretching.from_extreme.opt_and_forces` (-S here
refers to these subjobs)

  -c  Use this flag to run in a cluster. When -p is not defined, and you run in
      a slurm job manager system, the number of processors is equal to the
      number of cores asked in the submission of the job.
  -i  <index1,index2> indexes of the atoms used for constraining the distance
      of the intermedia structures when optimizing. If this flag is not used
      and there is a pdb file with the same <name> as defined with the flag -n,
      indexes 1 and 2 will correspond to the CH3 atoms in ACE and NME residues
      defined in the pdb if they exist.
  -l  <xc,base="bmk,6-31+g"> level of DFT theory.
  -n  <name> this script will collect all the *<name>*.xyz files sorted in
      alphabetic order as initial path. If this is a directory, this script
      will create a subdirectory there called 'conopt' and copy the
      *<name>*.xyz files to that subdirectory first.
  -p  <processors=1> number of processors per gaussian job. See description of
      flag -c.
  -S  <job_options=''> options for submitting a new job. This flag only makes
      sense in slurm cluster. Please, do not include a name (-J), nor the
      number of cores (-n, use -p for this). The input should be as in the next
      example: \"--partition=cpu --nice\".

  -h  prints this message.
EOF

exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
cluster='false'
indexes=''
level="bmk,6-31+g"
n_processors=''
job_options=''
verbose=''
while getopts 'ci:l:n:p:S:vh' flag;
do
  case "${flag}" in
    c) cluster='true' ;;
    i) indexes=${OPTARG} ;;
    l) level=${OPTARG} ;;
    n) name=${OPTARG} ;;
    p) n_processors=${OPTARG} ;;
    S) job_options=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done


source $(sith basics -path) "ContiPath" $verbose

# starting information
verbose -t "JOB information"
verbose -t "==============="
verbose -t " \* Date:"
verbose -t $(date)
verbose -t " \* Command:"
verbose -t "$0" "$@"

c_flag=''
if $cluster
then
  load_modules
  c_flag='-c'
  if [[ -z "$n_processors" ]] 
  then
    if [[ ! -z "$SLURM_CPUS_ON_NODE" ]]
    then
      n_processors=$SLURM_CPUS_ON_NODE
    else
      n_processors=1
    fi
  fi
  opt_forces_job_options="$job_options"
  job_options="$job_options -n $n_processors"
  job_options="sbatch $job_options"
fi

# C-CAP indexes in gaussian convention
if [[ -z "$indexes" ]] && [[ -f "${name%.*}.pdb" ]]
then
  # reading indexes from pdb file
  index1=$( grep ACE "$mol.pdb" | grep CH3 | awk '{print $2}' )
  index2=$( grep NME "$mol.pdb" | grep CH3 | awk '{print $2}' )
  indexes="$index1,$index2"
fi

# reading indexes from user input
index1=$( echo "$indexes" | cut -d ',' -f 1 )
index2=$( echo "$indexes" | cut -d ',' -f 2 )


# check that the indexes were read properly:
[[ -z "$index1" ]] && fail "Not recognized indexes: index 1: '$index1',
  index 2: '$index2' from indexes: '$indexes'"
[[ -z "$index2" ]] && fail "Not recognized indexes: index 1: $index1, index 2:
  $index2 from indexes: $indexes"
[[ "$index1" == "$index2" ]] && fail "Not recognized indexes: index 1:
  '$index1', index 2: '$index2' from indexes: '$indexes'"

xc_functional=$(echo $level | cut -d ',' -f 1)
basis_set=$(echo $level | cut -d ',' -f 2)

# ---- BODY -------------------------------------------------------------------
if [ -d "$name" ]
then
  cd $name
  mkdir conopt
  if [[ "${name: 0: 2}" == "./" ]]
  then
    name=${name#./}
  fi
  name=${name%/}
  cp *"${name}"*.xyz conopt
  cd conopt
fi

verbose "collecting and renaming xyz files"
j=0
for xyz_file in $(ls *$name*.xyz | sort )
do
  if [[ "$xyz_file" != "$name-conopt$(printf "%03d" $j).xyz" ]]
  then
    mv $xyz_file $name-conopt$(printf "%03d" $j).xyz
  fi
  verbose -t "$xyz_file --> $name-conopt$(printf "%03d" $j).xyz"
  j=$(( j + 1 ))
done

verbose "Extract dofs from xyzs. output: $name-conopt<n>-dofs.dat" 
sith extr_dofs -f ${name}-conopt > /dev/null || \
  fail "extracting dofs from xyzs"
# reduce irrelevant changes and add intermedias when the changes are too large,
# store the new subset in a dir called 'subset'.
verbose "Reduce structures"
sith reduce_structs "." ${name}-conopt > /dev/null || \
  fail "reducing structures"

# ==== Create com gaussian files
verbose "create template"
# Create templete first
tail -n +3 ${name}-conopt000.xyz > tmp.xyz
newzmat -ixyz -ozmat \
        -rebuildzmat tmp template > /dev/null || \
  fail "executing newzmat"
rm tmp.xyz

sed -i "s/-- No Title Specified --/Computation of forces/g" template.com
sed -i "s/\# HF\/6-31G\* Test/%chk=replaced_later/g" template.com
sed -i "/%chk/a %NProcShared=$n_processors" template.com

sed -i "/%NProcShared/a #P ${xc_functional}\/${basis_set} opt(modredun)" template.com
sed -i "1a %mem=60000MB" template.com

# import xyz files of the subset
rm ${name}-conopt*.dat
mv subset/* .
rm -r subset

# Create .com files
verbose "Create com files and submitting job"

sith find_blocks -f template.com -e "Variables:" -o tmp $verbose > /dev/null
mv tmp_000.out heading_template.out

sith find_blocks -s "\^\$" -e "'^ $'" \
                   -f template.com -o tmp  $verbose > /dev/null
mv tmp_002.out connectivity_template.out
rm tmp*.out

create_bck forces
mkdir -p forces

str_index=0
speficic_job_options=''
for file in ${name}-conopt*.dat
do
  str_index=$(( 10#$str_index + 1 ))
  struct_name=${file%.dat}
  cp heading_template.out  $struct_name.com
  sed -i "/chk=/c %chk=$struct_name" $struct_name.com
  echo "     Variables:" >> $struct_name.com
  cat $file >> $struct_name.com
  echo "" >> $struct_name.com
  echo "$index1 $index2 F" >> $struct_name.com
  echo "" >> $struct_name.com
  verbose -t "-  $struct_name $str_index"
  [ -z "$job_options" ] || \
    speficic_job_options="$job_options -J $(printf "%03d" \
                          $str_index)O$struct_name"
  $speficic_job_options \
    $( sith opt_and_forces -path ) $c_flag -f "$struct_name" \
                                   -p "$n_processors" \
                                   -S "$opt_forces_job_options" \
                                   $verbose
done

rm heading_template.out
rm connectivity_template.out
rm template.com
rm *.dat

finish
