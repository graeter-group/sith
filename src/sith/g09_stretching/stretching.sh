#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool obtains the stretched configurations of a molecule by increasing the
distance between two atoms, constraining and optimizing at every step.

  -b  <number_of_breakages=1> The simulation will run until get this number of
      ruptures.
  -c  Use this flag to run in a cluster set with slurm. In that case, the
      of processors is equal to the number of cores asked in the submission of
      the job. When this flag is present, -p is ignored.
  -e  <extend_method=0> index of stretching method. To see the options, use
      'sith change_distance -h' to see the order.
      carbons of the capping groups
  -i  <index1,index2> indexes of the atoms to use for increasing the distance.
  -l  <xc,base="bmk,6-31+g"> evel of DFT theory.
  -m  <molecule> molecule name. In this directory, a file called
      <molecule>-stretched00.pdb must exist.
  -p  <processors=1> number of processors per gaussian job.
  -r  restart stretching. In this case, this code must be executed from
      the molecule's directory.
  -s  <size[A]=0.2> Size of the step that increases the distances at each step.

  -v  verbose
  -h  prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
breakages=1
method=0
restart='false'
size=0.2
verbose=''
indexes=''
level="bmk,6-31+g"
cluster='false'
n_processors=1
retake='true'
while getopts 'b:ce:i:l:m:p:rs:vh' flag; do
  case "${flag}" in
    b) breakages=${OPTARG} ;;
    c) cluster='true' ;;
    e) extend_method=${OPTARG} ;;
    i) indexes=${OPTARG} ;;
    l) level=${OPTARG} ;;
    m) mol=${OPTARG} ;;
    p) n_processors=${OPTARG} ;;
    r) restart='true' ;;
    s) size=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" STRETCHING $verbose

if $cluster
then
  load_modules
  n_processors=$SLURM_CPUS_ON_NODE
fi

# starting information
verbose -t "JOB information"
verbose -t "==============="
verbose -t " * Date:"
verbose -t $(date)
verbose -t " * Command:"
verbose -t "$0" "$@"

# stretching method
if [[ "$extend_method" -eq 0 ]]
then
  if [ "$breakages" -eq 1 ]
  then
    method='scale_distance'
  else
    method='increase_distance_with_constraints'
  fi
fi
verbose "The stretching method will be '$method'"

# level of theory
xc_functional=$(echo $level | cut -d ',' -f 1)
basis_set=$(echo $level | cut -d ',' -f 2)

# check dependencies
ase -h &> /dev/null || fail "This code needs ASE"

mol_file=$mol
mol=${mol_file%.*}
# ----- set up finishes -------------------------------------------------------

# ---- BODY -------------------------------------------------------------------
# ----- checking restart ------------------------------------------------------
if $restart
then
  [[ -f 'frozen_dofs.dat' ]] || fail "frozen_dofs.dat doesn't exist"

  # extracting last i with xyz file already created
  mapfile -t previous < <( find . -maxdepth 1 -type f -name "$mol*.xyz" \
                                  -not -name "*bck*" | sort )
  
  [ ${#previous} -eq 0 ] && fail "Non previous xyz files were found"

  wext=${previous[-1]}
  last=${wext%.*}
  if [ "${last: -1}" == 'a' ]
  then
    last=${last::-2}
  fi
  i=$(( 10#${last:0-2} ))
  verbose "Restarting $mol, searching last optimization, $i is the last
           stretching step detected"

  # searching incomplete optimization trials
  nameiplusone=$(printf "%02d" "$(( i + 1 ))")

  # searching advances in i+1.
  # retake='true' means that the last configuration was not taken from an
  # incomplete job. In the next block, we search for advances in i+1 and if it
  # finds one, it restarts from there and sets retake='false' as a consecuence,
  # this variable is used later in the loop.
  sith log2xyz "$mol-stretched${nameiplusone}.log" 2> /dev/null && \
    create_bck "$mol-stretched${nameiplusone}."* &&
    lastone=$( search_last_bck $mol-stretched${nameiplusone} ) &&
    if [ $(( lastone )) -gt 2 ]; then fail "this optimization was" \
      "restarted more than 3 times and didn't converged."; fi &&
    echo "coping $mol-stretched${nameiplusone}-bck_$lastone.xyz" &&
    sith change_distance "$mol-stretched${nameiplusone}-bck_$lastone.xyz" \
      "$mol-stretched${nameiplusone}" frozen_dofs.dat 0 0 "$method" \
      --xc "'$xc_functional'" --basis "'$basis_set'" && \
    retake='false' && \
    warning "The stretching of molecule $mol will be restarted
      from $(( i + 1 ))"

  # if i+1 trial doesn't exist
  $retake && \
      warning "The stretching of molecule $mol will be restarted from $i"
else
  # C-CAP indexes in gaussian convention
  if [[ -z "$indexes" ]] && [[ "${mol_file##*.}" == 'pdb' ]]
  then
    # reading indexes from pdb file
    index1=$( grep ACE "$mol.pdb" | grep CH3 | awk '{print $2}' )
    index2=$( grep NME "$mol.pdb" | grep CH3 | awk '{print $2}' )
  else
    # reading indexes from user input
    index1=$( echo "$indexes" | cut -d ',' -f 1 )
    index2=$( echo "$indexes" | cut -d ',' -f 2 )
  fi

  # check that the indexes were read properly:
  [[ "$index1" -eq 0 && "$index2" -eq 0 ]] && fail "Not recognized indexes"
  [[ "$index1" -eq 1 && "$index2" -eq 1 ]] && fail "Not recognized indexes"
  [[ "$index1" == "$index2" ]] && fail "Not recognized indexes"

  if ! [[ -f 'frozen_dofs.dat' ]]
  then
    echo "$index1 $index2 F" > frozen_dofs.dat
  fi

  verbose "This code will stretch the atoms with the indexes $index1 $index2
  (1-based indexing)"

  # in case of not restarting
  i=-1
fi
# ----- checking restart finishes ---------------------------------------------

# ----- stretching starts -----------------------------------------------------
verbose "Stretching of $mol starts and will run until getting $breakages
  ruptures"

while [[ "$( wc -l < "frozen_dofs.dat" )" -le "$breakages" ]]
do
  # names by index
  namei=$(printf "%02d" "$i")
  nameiplusone=$(printf "%02d" $(( i + 1)))
  nameiplustwo=$(printf "%02d" $(( i + 2)))

  verbose "Stretched ${nameiplusone} starts"
  # Creates gaussian input file
  if [ $(( i + 1)) -eq 0 ]
  then
    # initial gaussian optimization
    verbose "The first gaussian process is an optimization"
    sith change_distance "$mol_file" \
      "$mol-stretched00" frozen_dofs.dat 0 0 "$method" \
      --xc "'$xc_functional'" --basis "'$basis_set'" || \
      fail "Creating initial gaussian input"
    sed -i  '/^TV  /d' "$mol-stretched00.com"
    sed -i "/opt/d" "$mol-stretched00.com"
  else
    if "$retake"
    then
      # if the last configuration was not taken from an incomplete job
      sith change_distance \
        "$mol-stretched$namei.xyz" "$mol-stretched${nameiplusone}" \
        frozen_dofs.dat "$size" 0 "$method" \
        --xc "'$xc_functional'" --basis "'$basis_set'"\
        || fail "Creating gaussian input of stretching $nameiplusone"
    fi
    retake='true'
    sed -i '$d' "$mol-stretched${nameiplusone}.com"
    # add constrains
    cat frozen_dofs.dat >> \
    "$mol-stretched${nameiplusone}.com"  
  fi
  sed -i "1a %NProcShared=$n_processors" "$mol-stretched${nameiplusone}.com"
  sed -i "/#P/a opt(modredun,calcfc)" "$mol-stretched${nameiplusone}.com"

  # run gaussian
  verbose "Running optmization of stretching ${nameiplusone}"
  gaussian "$mol-stretched${nameiplusone}.com" \
      "$mol-stretched${nameiplusone}.log" || \
    { if [ "$(grep -c "Atoms too close." \
            "$mol-stretched${nameiplusone}.log")" \
            -eq 1 ]; then fail "Atoms too close for ${nameiplusone}" ; \
      fi ; }

  # check convergence from output
  output=$(grep -i optimized "$mol-stretched${nameiplusone}.log" | \
           grep -c -i Non )
  if [ "$output" -ne 0 ]
  then
    # If the code enters here is because, in a first optimization, it
    # didn't converge so it has to run again to get the optimized 
    # structure. As a second chance to converge.
    verbose "Optimization did not converge with distance $(( i + 1 )) *
      $size . Then, a new trial will start now"
    sith log2xyz "$mol-stretched${nameiplusone}.log" || fail "
      Transforming log file to xyz in second trial of optimization"
    sith change_distance \
            "$mol-stretched${nameiplusone}.xyz" \
            "$mol-stretched${nameiplustwo}" frozen_dofs.dat 0 0 "$method" \
            --xc "'$xc_functional'" --basis "'$basis_set'" || fail "changing distance"
    # save the failed files in ...-stretched<number>a.*
    create_bck "$mol-stretched${nameiplusone}"*
    # then restart the optimization
    mv "$mol-stretched${nameiplustwo}.com" "$mol-stretched${nameiplusone}.com"
    sed -i "s/stretched${nameiplustwo}/stretched${nameiplusone}/g" \
      "$mol-stretched${nameiplusone}.com"
    sed -i "1a %NProcShared=$n_processors" "$mol-stretched${nameiplusone}.com"
    sed -i "/#P/a opt(modredun,calcfc)" "$mol-stretched${nameiplusone}.com"
    sed -i '$d' "$mol-stretched${nameiplusone}.com"
    cat frozen_dofs.dat >> \
      "$mol-stretched${nameiplusone}.com"
    # run optimization
    verbose "Re-running optimization"
    gaussian "$mol-stretched${nameiplusone}.com" \
        "$mol-stretched${nameiplusone}.log"
  fi

  # check the output again
  output=$(grep -i optimized "$mol-stretched${nameiplusone}.log" | \
           grep -c -i Non )
  [ "$output" -ne 0 ] && failed "Optimization when the stretched distance was
      $(( i + 1 ))*0.2 didn't converge. No more stretching will be applied"

  # creating fchk file
  formchk -3 "$mol-stretched${nameiplusone}.chk" || fail "Creating fchk file"

  # ==== Testing DOFs
  verbose "Testing dofs"
  sith log2xyz "$mol-stretched${nameiplusone}.log" || fail "Transforming
    log file to xyz"
  
  if [ "$i" -eq -1 ]
  then
    extrad=".."
  else
    # check the bonds in 'i' that are not in 'i+1'. If any, they are added to
    # frozen_dofs.dat
    extrad=$( sith diff_bonds "$mol-stretched${namei}.xyz" \
              "$mol-stretched${nameiplusone}.xyz" )
  fi

  if [ ${#extrad} -ne 2 ]
  then
    # create bck in rupture directory in case it already exists
    verbose "rupture found in $extrad, this bond will be frozen"
    mkdir -p rupture ; cd rupture
    create_bck "$mol-stretched${nameiplusone}"
    cd .. 
    mv "$mol-stretched${nameiplusone}"* rupture/
  else
     verbose "Non-rupture detected in stretched ${nameiplusone}"
  fi

  verbose "Stretched ${nameiplusone} finished"
  # next i
  i=$(( i + 1 ))
done

# ----- stretching finishes ---------------------------------------------------

finish "$mol finished"
