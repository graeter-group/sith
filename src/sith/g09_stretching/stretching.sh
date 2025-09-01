#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool obtains the stretched configurations of a peptide by increasing the
distance between carbons of the capping groups, constraining and optimizing
using BMK exchange-correlation.

  -b  <number of breakages=1> The simulation will run until get this number of
      ruptures.
  -p  <peptide> One letter code of the amino acids forming the peptides. In
      this directory, a file called <peptide>-stretched00.pdb has to exist.
  -m  <method=0> index of stretching method. To see the options, use
      'sith change_distance -h' to see the order.
  -r  restart stretching. In this case, this conde must be executed from
      the peptide's directory.
  -s  <size[A]=0.2> of the step that increases the distances.

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
n_processors=8
restart='false'
size=0.2
verbose='false'
level='bmk,6-31+g'
while getopts 'b:i:p:l:m:rs:vh' flag; do
  case "${flag}" in
    b) breakages=${OPTARG} ;;
    c) change_method=${OPTARG} ;;
    l) level=${OPTARG} ;;
    r) restart='true' ;;
    s) size=${OPTARG} ;;
    
    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" STRETCHING $verbose

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

# stretching method
if [[ "$change_method" -eq 0 ]]
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
command -V g09 &> /dev/null || fail "This code needs gaussian"

# C-CAP indexes in g09 convention
index1=$( grep ACE "$pep-stretched00.pdb" | grep CH3 | awk '{print $2}' )
index2=$( grep NME "$pep-stretched00.pdb" | grep CH3 | awk '{print $2}' )
# check that the indexes were read properly:
[[ "$index1" -eq 0 && "$index2" -eq 0 ]] && fail "Not recognized indexes"
[[ "$index1" -eq 1 && "$index2" -eq 1 ]] && fail "Not recognized indexes"
if ! [[ -f 'frozen_dofs.dat' ]]
then
  echo "$index1 $index2 F" > frozen_dofs.dat
fi
verbose "This code will stretch the atoms with the indexes $index1 $index2
  (g09 convention)"

# ----- set up finishes -------------------------------------------------------

# ---- BODY -------------------------------------------------------------------
# ----- checking restart starts -----------------------------------------------
if $restart
then
  # extracting last i with xyz file already created
  mapfile -t previous < <( find . -maxdepth 1 -type f -name "*$pep*.xyz" \
                          -not -name "*bck*" | sort )
  wext=${previous[-1]}
  last=${wext%.*}
  if [ "${last: -1}" == 'a' ]
  then
    last=${last::-2}
  fi
  i=$(( 10#${last:0-2} ))
  verbose "Restarting $pep, searching last optimization, $i is the last
           stretching detected"

  # searching incomplete optimization trials
  nameiplusone=$(printf "%02d" "$(( i + 1 ))")

  # searching advances in i+1.
  # retake='true' means that the last configuration was not taken from an
  # incomplete job. In the next block, we search for advances in i+1 and if it
  # finds one, it restarts from there and sets retake='false' as a consecuence,
  # this variable is used later in the loop.
  retake='true'
    if [ $(( lastone )) -gt 2 ]; then fail "this optimization was" \
      "restarted more than 3 times and didn't converged."; fi &&
    echo "coping $pep-stretched${nameiplusone}-bck_$lastone.xyz" &&
    sith change_distance "$pep-stretched${nameiplusone}-bck_$lastone.xyz" \
      "$pep-stretched${nameiplusone}" frozen_dofs.dat 0 0 \
      --xc $xc_functional --basis $basis_set && \
    retake='false' && \
    warning "The stretching of peptide $pep will be restarted
      from $(( i + 1 ))"

  # if i+1 trial doesn't exist
  $retake && \
      warning "The stretching of peptide $pep will be restarted from $i"
else
  # in case of not restarting
  i=-1
fi
# ----- checking restart finishes ---------------------------------------------

# ----- stretching starts -----------------------------------------------------
verbose "Stretching of $pep starts and will run until getting $breakages
  ruptures"

while [[ "$( wc -l < "frozen_dofs.dat" )" -le "$breakages" ]]
do
  # names by index
  namei=$(printf "%02d" "$i")
  nameiplusone=$(printf "%02d" $(( i + 1)))
  nameiplustwo=$(printf "%02d" $(( i + 2)))

  verbose "Stretched ${nameiplusone} starts"
  # Creates g09 input file
  if [ $(( i + 1)) -eq 0 ]
  then
    # initial g09 optimization
    verbose "The first g09 process is an optimization"
      --xc $xc_functional --basis $basis_set || \
      fail "Preparating the input of gaussian"
    sed -i  '/^TV  /d' "$pep-stretched00.com"
    sed -i "/opt/d" "$pep-stretched00.com"
  else
    if "$retake"
    then
      sith change_distance \
        "$pep-stretched$namei.xyz" "$pep-stretched${nameiplusone}" \
        frozen_dofs.dat "$size" 0 "$method" \
        --xc $xc_functional --basis $basis_set\
        || fail "Preparating g09 input"
    fi
    retake='true'
    sed -i '$d' "$pep-stretched${nameiplusone}.com"
    # add constrains
    cat frozen_dofs.dat >> \
    "$pep-stretched${nameiplusone}.com"  
  fi
  sed -i "1a %NProcShared=$n_processors" "$pep-stretched${nameiplusone}.com"
  sed -i "/#P/a opt(modredun,calcfc)" "$pep-stretched${nameiplusone}.com"

  # run gaussian
  verbose "Running optmization of stretching ${nameiplusone}"
  g09 "$pep-stretched${nameiplusone}.com" \
      "$pep-stretched${nameiplusone}.log" || \
    { if [ "$(grep -c "Atoms too close." \
            "$pep-stretched${nameiplusone}.log")" \
            -eq 1 ]; then fail "Atoms too close for ${nameiplusone}" ; \
      fi ; }

  # check convergence from output
  output=$(grep -i optimized "$pep-stretched${nameiplusone}.log" | \
           grep -c -i Non )
  if [ "$output" -ne 0 ]
  then
    # If the code enters here is because, in a first optimization, it
    # didn't converge so it has to run again to get the optimized 
    # structure. As a second chance to converge.
    verbose "Optimization did not converge with distance $(( i + 1 )) *
      $size . Then, a new trial will start now"
    sith log2xyz "$pep-stretched${nameiplusone}.log" || fail "
      Transforming log file to xyz in second trial of optimization"
    sith change_distance \
            "$pep-stretched${nameiplusone}.xyz" \
            "$pep-stretched${nameiplustwo}" frozen_dofs.dat 0 0 \
            --xc $xc_functional --basis $basis_set || fail "changing distance"
    # save the failed files in ...-stretched<number>a.*
    create_bck "$pep-stretched${nameiplusone}"*
    # then restart the optimization
    mv "$pep-stretched${nameiplustwo}.com" "$pep-stretched${nameiplusone}.com"
    sed -i "s/stretched${nameiplustwo}/stretched${nameiplusone}/g" \
      "$pep-stretched${nameiplusone}.com"
    sed -i "1a %NProcShared=$n_processors" "$pep-stretched${nameiplusone}.com"
    sed -i "/#P/a opt(modredun,calcfc)" "$pep-stretched${nameiplusone}.com"
    sed -i '$d' "$pep-stretched${nameiplusone}.com"
    cat frozen_dofs.dat >> \
      "$pep-stretched${nameiplusone}.com"
    # run optimization
    verbose "Re-running optimization"
    g09 "$pep-stretched${nameiplusone}.com" \
        "$pep-stretched${nameiplusone}.log"
  fi

  # check the output again
  output=$(grep -i optimized "$pep-stretched${nameiplusone}.log" | \
           grep -c -i Non )
  [ "$output" -ne 0 ] && failed "Optimization when the stretched distance was
      $(( i + 1 ))*0.2 didn't converge. No more stretching will be applied"

  # Testing DOFs
  verbose "Testing dofs"
  sith log2xyz "$pep-stretched${nameiplusone}.log" || fail "Transforming
    log file to xyz"
  
  if [ "$i" -eq -1 ]
  then
    extrad=".."
  else
    # Add extra values to frozen
    extrad=$( sith diff_bonds "$pep-stretched${namei}.xyz" \
              "$pep-stretched${nameiplusone}.xyz" )
  fi

  if [ ${#extrad} -ne 2 ]
  then
    # if a rupture is detected, this dof is detected and a new state 
    verbose "rupture found in $extrad, this bond will be frozen"
    if ! [[ -d rupture ]]
    then
      mkdir rupture
    fi
    cd rupture || fail "rupture directory not found"
    create_bck "$pep-stretched${nameiplusone}"
    cd .. || fail "moving to back directory"
    mv "$pep-stretched${nameiplusone}"* rupture/
  else
     verbose "Non-rupture detected in stretched ${nameiplusone}"
  fi

  verbose "Stretched ${nameiplusone} finished"

  # creating fchk file
  formchk -3 "$pep-stretched${nameiplusone}.chk" || fail "Creating fchk file"
  # next i
  i=$(( i + 1 ))
done

# ----- stretching finishes ---------------------------------------------------

finish "$pep finished"
