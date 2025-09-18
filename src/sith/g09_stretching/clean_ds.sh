#!/usr/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This code checks the aminoacids in the 'data runnig' directory and evaluates
the state of the running (running, completed, failed). In case to be completed,
it is moved to the dataset directory.

  -r  running directory. default ./
  -s  data set directory. default ../random3
  -u  user. defatult sucerqdl

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
running_directory="./"
dataset="../random3"
user="sucerqdl"
verbose=''
while getopts 'r:s:u:vh' flag;
do
  case "${flag}" in
    r) running_directory=${OPTARG} ;;
    s) dataset=${OPTARG} ;;
    u) user=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(sith basics -path)" CHECK_DS $verbose 

verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

# ---- BODY -------------------------------------------------------------------

# requierements
if [ ! -d "$dataset" ]
then
  fail "The directory specified as data set does not exist. check
    'sith clean_ds -h' for more info"
fi

if [ ! -d "$running_directory" ]
then
  fail "The directory specified as running directory does not exist. check
    'sith clean_ds -h' for more info"
fi

mapfile -t job_ids < <( squeue -u "$user" | tail -n +2 |\
                        awk '{print $1}' ) || \
                        fail "You have to run this code in a cluster using
                          slurm and specify an existing user with -u option."

# ==== creating output
echo "# amino state (R:running, C: completed, E: error)" > 00-aminos.txt

# ---- collect running peptides
running_pep=( )
for id in "${job_ids[@]}"
do
  if  ! find . -maxdepth 1 -type f -name "*$id.o" > /dev/null 2>&1
  then
    fail "non peptide output found. This code is expeting to find a file
      with the id of each job in the running directory"
  else
    random="$(grep -A 1 Command "peptides_analysis-$id.o" | tail -n 1 | \
              grep -c "\-R")"
    if [ "$random" -eq 0 ]
    then
      com="$(grep -A 1 Command "peptides_analysis-$id.o" | tail -n 1 )"
      pep=$(echo "${com##*p }" | cut -d ' ' -f 1)
    else
      string=$(grep "WARNING The code created the random peptide " \
                    "peptides_analysis-$id.o" )
      pep=$(echo "${string#*peptide }" | cut -d ',' -f 1)
    fi

    if [ "${#pep}" -eq 0 ]
    then
      fail "$id  is running but the output (peptides_analysis-<id>.o)
        does not exist or does not show the peptide yet."
    fi
    running_pep+=( "$pep" )
    echo "$pep" "   R" >> 00-aminos.txt
    # move other files
    mapfile -t other_files < <(find . -maxdepth 1 -type f -name "*.o" \
                                      -exec grep -l "$pep" {} \;)
    for fil in "${other_files[@]}"
    do
      if [ "$fil" !=  "./peptides_analysis-$id.o" ]
      then
        verbose "moving $fil into $pep"
        mv "$fil" "$pep"
        mv "${fil%*.o}.e" "$pep"
      fi
    done
  fi
done

# add  peptides already in dataset
mapfile -t completed_pep < <(ls "$dataset")
for pep in "${completed_pep[@]}"
do
  echo "$pep" "   C" >> 00-aminos.txt
done

# ---- collect non-running peptides
verbose "Non running peptides"
mapfile -t peps_here < <( find . -maxdepth 1 -type d | tail -n +2)

for pep in "${peps_here[@]}"
do
  # shellcheck disable=SC2199
  if [[ ! "${running_pep[@]}" =~ ${pep:2} ]]
  then
    # move files into the proper directory
    mapfile -t out_files < <(find . -maxdepth 1 -type f -name "*.o" \
                                    -exec grep -l "${pep:2}" {} \;)
    for fil in "${out_files[@]}"
    do
      mv "$fil" "$pep"
      mv "${fil%*.o}.e" "$pep"
    done
    if [ -d "$pep/rupture" ]
    then
      warning "moving $pep to the dataset directory"
      mv "$pep" "$dataset" || fail "trying to move ${pep:2} to $dataset"
      echo "${pep:2}" "   C" >> 00-aminos.txt
    else
      warning "something happened with the stretching of peptide
        ${pep:2}"
      echo "${pep:2}" "   E" >> 00-aminos.txt
    fi
  fi
done

# Now separates the non-running by completed and error. Understanding error as
# anything that aviods to complete

verbose "summary"
echo "completed Jobs: $(grep -c "  C" 00-aminos.txt )"
echo "running Jobs: $(grep -c "  R" 00-aminos.txt )"
echo "failed Jobs: $(grep -c "  E" 00-aminos.txt )"

finish
