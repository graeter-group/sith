#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e


# ----- definition of functions -----------------------------------------------
print_help() {
echo "
creates the com file to restart a gaussian09 job from an existing
.chk file assuming that the com file of the process exists. I basically
uses the same com file but replacing the geometry by the one in the .chk file.

  -c  run un a cluster.
  -d  <directory=./> directory where to find the .chk and .com files.
  -n  <name> name (without extension) of the .chk and .com files.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
dir='./'
name=''
cluster='false'
verbose=''
while getopts 'cd:n:vh' flag;
do
  case "${flag}" in
    c) cluster='true' ;;
    d) dir=${OPTARG} ;;
    n) name=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith restart_g09 -h" >&2 ; exit 1 ;;
  esac
done

# shellcheck disable=SC1090
source "$(sith basics -path)" Restartg09 "$verbose"

# starting information
verbose -t "JOB information"
verbose -t "==============="
verbose -t " Command:" "$0" "$@"

# load modules
if $cluster
then
  echo load modules for this script
fi

# ---- BODY -------------------------------------------------------------------
verbose "important separation"

cp "$file.com" "tmp-$file.com"
create_bck "$file.com"
n = search_last_bck

mv "tmp-$file.com" "$file.com"
if ! grep -q "Geom=Check" "$file.com"
then
  sed -i "/#/a Guess=Read Geom=Check/g" "$file.com"
  n_i=$(grep -n "^$" "$file.com" | cut -d : -f 1 | sed -n '2p')
  n_j=$(grep -n "^$" "$file.com" | cut -d : -f 1 | sed -n '3p')
  sed -i "$(( n_i + 2 )),$(( n_j - 1 ))d" "$file.com"
fi

# ---- END --------------------------------------------------------------------
finish "message to finish"
