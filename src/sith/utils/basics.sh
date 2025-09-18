#!/bin/bash

print_help() {
echo -e "
This script contains basics tools for other bash scripts. With the code in
here, you can use the next functions

- <first_argument> Add a label in every output of each bash script by giving
  it as a first argument.
- <second_argument> Add anything in order to have verbose printing in the
  output. otherwise, the your script will ignore all the verbose.
- adjust <text>: add the label and print what <text> in an adjusted column of
  80 characters.
- verbose <text>: besides of the label, it adds the keyword VERBOSE in to the
  begining of <text> and print the adjusted text.
- warning <text>: besides of the label, it adds the keyword WARNING in to the
  begining of <text> and print the adjusted text.
- finish <text>: besides of the label, it use verbose to print <text> of the
  word "finish" if <text> is not given. It also stops the script with 'exit 0'
- fail <text>: besides of the label, it adds the keyword VERBOSE in to the
  begining of <text> and print the adjusted text. It also print the message
  in the std_error and stops the script with 'exit 1'
- create_bck [<name1> <name2> ...]: function that moves an existing file or
  directory to <basic_functions_name>-bck_[n][.ext] where n is the number of
  the backup with 3 digits (leading zeros if necessary) and ext is
  automatically extracted from the original file.
- search_last_bck <name>: finds the last [n] created by the function
  create_bck.
- load_modules: loads the modules I need in my codes when running in a cluster.
- wait_until_next_file_exist <file_name>: literaly waits untils <file_name>
  appears after executing ls. This is important for some file systems that
  take a little while to recognize files created in a script. 
"

exit 0
}

while getopts 'h' flag;
do
  case "${flag}" in
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done


# Definition functions and variables that are used along the whole package.

# ------ variables ------------------------------------------------------------
array_bfnames=( "$1" "${array_bfnames[@]}" )
basic_functions_name=${array_bfnames[0]}

if [[ "$2" == "true" ]]
then
  eval "BASICVERBOSE_${basic_functions_name[0]}=true"
else
  eval "BASICVERBOSE_${basic_functions_name[0]}=false"
fi

# ------ functions ------------------------------------------------------------
# Function that adjustes the text to 80 characters
adjust () {
  text="++++ ${basic_functions_name[0]}: $*"
  addchar=$(( 80 - ${#text} % 80 ))
  text="$text $( perl -E "say '+' x $addchar" )"
  nlines=$(( ${#text} / 80 ))
  for (( w=0; w<=nlines-1; w++ ))
  do
    echo "${text:$(( w * 79 )):79}"
  done
}

adjust_text () {
  text="$*"
  nlines=$(( ${#text} / 80 ))
  for (( w=0; w<=nlines; w++ ))
  do
    echo "${text:$(( w * 79 )):79}"
  done
}

# prints some text adjusted to 80 characters per line, filling empty spaces
# with +
verbose () {
  if [[ "$(eval "echo \$BASICVERBOSE_${basic_functions_name[0]}")" == "true" ]]
  then
    if [[ "$1" == "-t" ]]
    then
      shift
      # shellcheck disable=SC2068
      adjust_text $@
    else
      # shellcheck disable=SC2068
      adjust "VERBOSE" $@ "$( date )"
    fi
  fi
}

warning () {
  # shellcheck disable=SC2068
  adjust "WARNING" $@ "$( date )" >&2
}

finish () {
  if [ "$#" -ne 0 ]
  then
    # shellcheck disable=SC2068
    verbose $@ "$( date )"
  else
    verbose finish "$( date )"
  fi
  echo
  array_bfnames=( "${array_bfnames[@]:1}" )
  basic_functions_name=${array_bfnames[0]}
  exit 0
}

# Function that returns the error message and stops the run if something fails.
fail () {
  # shellcheck disable=SC2068
  adjust "ERROR" $@ "$( date )"
  # shellcheck disable=SC2068
  finish "ERROR" $@ "$( date )" >&2
  exit 1
}

# function that moves an existing file or directory to 
# <basic_functions_name>-bck_n[.ext] where n is the number of the backup and
# ext is automatically extracted from the original file
create_bck () {
  for to_bck in "$@"
  do
    # in case creating backup directory
    bck=$to_bck-bck_001
    if [ -d "$to_bck" ]
    then
      bck_i=$(printf "%03d" 2)
      while [ -d "$bck" ]
      do
        bck=$to_bck-bck_$bck_i
        bck_i=$(printf "%03d" $(( 10#$bck_i + 1 )) )
      done
      warning "$to_bck directory already exist. This directory will be
        backed up in $bck"
      mv "$to_bck" "$bck"
    fi

    # in case creating backup file
    new_fil=${to_bck%.*} # file name
    ext=${to_bck##*.}    # file extension
    bck=$new_fil-bck_001.$ext
    if [ -f "$to_bck" ]
    then
      bck_i=$(printf "%03d" 2)
      while [ -f "$bck" ]
      do
        bck=$new_fil-bck_$bck_i.$ext
        bck_i=$(printf "%03d" $(( 10#$bck_i + 1 )) )
      done
      warning "$to_bck file already exist. This directory will be
        backed up in $bck"
      mv "$to_bck" "$bck"
    fi
  done
}

search_last_bck() {
  name_file=$1
  mapfile -t all_bcks < <( ls -1 "$1"-bck_???.* | sort )
  last_woext=${all_bcks[-1]%.*}
  # prints the number of the last config
  echo ${last_woext:0-3}
  name_file=$1
  mapfile -t all_bcks < <( ls -1 "$1"-bck_???.* | sort )
  last_woext=${all_bcks[-1]%.*}
  # prints the number of the last config
  echo ${last_woext:0-3}

}

load_modules() {
  # This function load the modules softwares needed for the different
  # calculations.
  args=$@
  if [ ${#args} -ne 0 ]
  then
    resubmit $args &
  fi
  echo " * This JOB will be run in the Node:"
  echo "$SLURM_JOB_NODELIST"

  source $HOME/sw/load_modules.sh
}

wait_until_next_file_exist() {
  while ! ls | grep -q $1
  do
    continue
  done
}  

verbose "STARTS"
