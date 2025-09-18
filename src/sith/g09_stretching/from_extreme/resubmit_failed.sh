#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Takes the gaussian log files given by -l flag, and if the log file does not
report a proper termination, a new job is resubmitted (creating a backup first)
using the command given by the flag -e.

  -d  <path='./'> directory where the log and com files are.
  -e  <exec='$( sith single_g09 -path ) -c -f '> execution command to be
      resubmitted. This command is very important to be given inside of \" such
      that it is understood as only one value.
  -f  <frozen=''> line to freeze dofs. eg: '2 5 F'.
  -c  <comfile> input file used to run the previous trial.
  -l  <logfile> log file used to run the previous trial.
  -j  <jobname=comfile without extension> name of the new Job to be resubmited.

  -v  verbose
  -h  prints this message.

Note: it assumes that ../frozen_dofs.dat exist.

Note: In principle, this can be done easier by the chk, but then the input is
not a Z-matrix anymore, then the output would not contain the info of the
forces.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
frozen_file='../frozen_dofs.dat'
verbose=''
directory='./'
toexecute=""
jobname=""
while getopts 'd:e:c:f:l:j:vh' flag;
do
  case "${flag}" in
    d) directory=${OPTARG} ;;
    e) toexecute=${OPTARG} ;;
    f) frozen_file=${OPTARG} ;;
    c) comfile=${OPTARG} ;;
    l) logfile=${OPTARG} ;;
    j) jobname=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: sith <function> -h" >&2 ; exit 1 ;;
  esac
done

if [[ $toexecute == "" ]]
then
  toexecute="$( sith single_g09 -path ) -c -f ${comfile%.com}"
fi

if [[ $jobname == "" ]]
then
  jobname=${comfile%.com}
fi

# ---- BODY -------------------------------------------------------------------
source "$(sith basics -path)" ResubmitFailed $verbose
origin_resub=$( pwd )
cd $directory

# save heading in comfile

# create backup of prev_trials
cp $logfile tmp_log
cp $comfile tmp_com
# Create tmp_first-block.com containing the execution block
sith find_blocks -f $comfile -e "^$"
mv output_001.out tmp_first-block.com
create_bck ${logfile%.log}.*
create_bck ${comfile%.log}.*
mv tmp_log $logfile
mv tmp_com $comfile
grep -iq "%mem" $comfile || sed -i "1i \%mem=60000MB" $comfile


# extract xyz
sith log2xyz "$logfile" || fail "Extracting xyz from logfile"
file=${logfile%.log}.xyz

# TODO: change this to just consider the frozen dofs when they are defined with the flag -f
sith shake_except $file $frozen_file

# create comfile
sith change_distance \
  $file ${file%.xyz} \
  "nofile" 0 0 "scale_distance" \
  || fail "Preparating gaussian input"
sith find_blocks -f ${file%.xyz}.com -s "^$"
mv tmp_first-block.com $comfile
echo "" >> $comfile
cat output_001.out >> $comfile
sed -i '$d' $comfile
head -n 1 $frozen_file >> $comfile
grep -iq "%mem" $comfile || sed -i "1i \%mem=60000MB" $comfile

# remove unnecessary
rm output_001.out
rm $logfile

if [[ "$(whoami)" == "hits_"* ]]
then
  single_part="--partition=cpu-single"
else
  single_part=""
fi

sbatch --job-name=$jobname $single_part \
  $toexecute || \
  fail "submitting $toexecute "

cd $origin_resub

finish
