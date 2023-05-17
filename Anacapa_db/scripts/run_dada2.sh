#!/bin/bash

set -x

### this script is run as follows
# ~/Anacapa_db/scripts/run_dada2.sh -o <out_dir> -d <database_directory> -m <metabarcode name> -b <min ASV read count acceptable>
OUT=""
DB=""
MB=""
TYP=""
MIN_ASV=""

while getopts "o:d:m:t:e:b:" opt; do
    case $opt in
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        m) MB="$OPTARG"  # metabarcode name
        ;;
        t) TYP="$OPTARG"  # type of reads
        ;;
        b) MIN_ASV="$OPTARG"
        ;;
    esac
done
####################################script & software

##load module
source $DB/scripts/anacapa_vars.sh

mkdir -p ${OUT}/Run_info/dada2_out
mkdir -p ${OUT}/${MB}/${MB}dada2_out
echo ""
echo "Running dada2 on ${TYP} reads"
echo "${MIN_ASV}"

##load module
# ${MODULE_SOURCE} # use if you need to load modules from an HPC

#### critical or the dependency 'RcppParallel' will not install
# ${R} &> ${OUT}/Run_info/dada2_out/dada2_out_${TYP}
# ${GCC} &>> ${OUT}/Run_info/dada2_out/dada2_out_${TYP}

Rscript --vanilla ${DB}/scripts/dada2_unified_script.R ${MB} ${OUT} ${TYP} ${MIN_ASV} &>> ${OUT}/Run_info/dada2_out/dada2_out_${TYP}

echo "moving on"
