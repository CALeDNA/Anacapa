#! /bin/bash

### this script is run as follows
# sh ~/Anacapa_db/scripts/anacapa_QC_dada2.sh -i <input_dir> -o <out_dir> -d <database_directory> -f <fasta file of forward primers> -r <fasta file of reverse primers> -u <HPC_user_name>, -g (add flag (-g), no text required, if fastq files are not compressed), -c change cut adapt error 3' and 5' trimming .3 default need value, -p change cut adapt error primer sorting / trimming .3 default need value, -q minimum quality score, -x additional base pairs trimmed from 5' end forward read, -y additional base pairs trimmed from 5' end reverse read, -b number of times an ASV must be found to retain after dada2, -e file path to minimum overlap for user determined F and R primers, -k path to custom HPC job submission header
HELP=""
IN=""
OUT=""
DB=""
UN=""
FP=""
RP=""
GUNZIPED=""
CTADE=""
PCTADE=""
QUALS=""
MILEN=""
FETRIM=""
RETRIM=""
MINTIMES_ASV=""
HPC_HEADER=""

while getopts "h?:i:o:d:u:f:r:g?:c:p:q:m:x:y:b:k:" opt; do
    case $opt in
        h) HELP="TRUE"
        ;;
        i) IN="$OPTARG" # path to raw .fastq.gz files
        ;;
        o) OUT="$OPTARG" # path to desired Anacapa output
        ;;
        d) DB="$OPTARG"  # path to Anacapa_db
        ;;
        u) UN="$OPTARG"  # need username for submitting sequencing job
        ;;
        f) FP="$OPTARG"  # need forward reads for cutadapt
        ;;
        r) RP="$OPTARG"  # need reverse reads for cutadapt
        ;;
        g) GUNZIPED="TRUE" #reads not compressed
        ;;
        c) CTADE="$OPTARG"  # cutadapt error for 3' adapter and 5' primer adapter trimming
        ;;
        p) PCTADE="$OPTARG"  # cutadapt error 3' primer sorting and trimming
        ;;
        q) QUALS="$OPTARG"  # Minimum Quality score
        ;;
        m) MILEN="$OPTARG"  # Minimum Length after trimming
        ;;
        x) FETRIM="$OPTARG"  # Additional 5' trimming forward read
        ;;
        y) RETRIM="$OPTARG"  # Additional 5' trimming reverse read
        ;;
        b) MINTIMES_ASV="$OPTARG" # number of occurances required to keep ASV
        ;;
        k) HPC_HEADER="$OPTARG" # path to custom HPC job submission header
        ;;
    esac
done

# if $DB has a / at the end remove it!
case "$DB" in
*/)
    DB=${DB%/}
    ;;
*)
    echo ""
    ;;
esac

if [ "${HELP}" = "TRUE" ]
then
  printf "<<< Anacapa: Sequence QC and ASV Parsing >>>\n\nThe purpose of these script is to process raw fastq or fastq.gz files from an Illumina HiSeq or MiSeq.  It removes 3' and 5' sequencing artifacts and 5' metabarcode primers (cutadapt), removes low quality base pairs and short reads (fastX-toolkit), sorts reads by 3' metabarcode primers prior to trimming (cutadapt), and uses dada2 to denoise, dereplicate, merge and remove chimeric reads\n\n	For successful implementation \n		1. Make sure you have all of the dependencies and correct paths in the anacapa_config.sh file\n		2. Add the Metabarcode locus specific CRUX reference libraries to the Anacapa_db folder\n		3. All parameters can be modified using the arguments below.  Alternatively, all parameters can be altered in the anacapa_vars.sh folder\n\nArguments:\n- Required for either mode:\n	-i	path to .fastq.gz files, if files are already compressed use flag -g (see below)\n	-o	path to output directory\n	-d	path to Anacapa_db\n	-t	Illumina Platform: HiSeq (2 x 150) or MiSeq ( >= 2 x 250)\n    \n- Optional:\n 	-u	If running on an HPC (e.g. UCLA's Hoffman2 cluster), this is your username: e.g. eecurd\n 	-f	path to file with forward primers in fasta format \n    		e.g.	>16s\n    			GTGYCAGCMGCCGCGGTAA\n			>18S\n			GTACACACCGCCCGTC\n	-r	path to file with forward primers in fasta format \n    		e.g. 	>16s\n    			GGACTACNVGGGTWTCTAAT\n    			>18S\n			TGATCCTTCTGCAGGTTCACCTAC\n	-g	If .fastq read are not compressed: -g (no argument need)\n	-c	To modify the allowed cutadapt error for 3' adapter and 5' primer adapter trimming: 0.0 to 1.0 (default 0.3)\n	-p	To modify the allowed cutadapt error 3' primer sorting and trimming: 0.0 to 1.0 (default 0.3)\n	-q	To modify the minimum quality score allowed: 0 - 40 (default 35)\n	-m	To modify the minimum length after quality trimming: 0 - 300 (default 100)\n	-x	To modify the additional 5' trimming of forward reads: 0 - 300 (default HiSeq 10, default MiSeq 20)\n	-y	To modify the additional 5' trimming of reverse reads: 0 - 300 (default HiSeq 25, default MiSeq 50)\n	-b	To modify the number of occurrences required to keep an ASV: 0 - any integer (default 0)\n	-e	File path to a list of minimum length(s) reqired for paired F and R reads to overlap \n		(length of the locus - primer length + 20 bp). The user should take into account variability in amplicon \n		region (e.g.The amplicon size for 18S 1389f-1510r is ~260 +/- 50 bp) and make appropriate allowances.\n		e.g.	LENGTH_16S="235"\n			LENGTH_18S="200"\n	-k	Path to file with alternate HPC job submission parameters:  \n		default file = ~/Anacapa_db/scripts/Hoffman2_HPC_header.sh\n		modifiable template file = ~/Anacapa_db/scripts/anacapa_qsub_templates.sh\n\n\n-Other:\n	-h	Shows program usage then quits\n\n\n\n"
  exit
else
  echo ""
fi

####################################script & software
# This pipeline was developed and written by Emily Curd (eecurd@g.ucla.edu), Jesse Gomer (jessegomer@gmail.com), Gaurav Kandlikar (gkandlikar@ucla.edu), and Baochen Shi (biosbc@gmail.com), and with contributions from Zack Gold (zack.j.gold@gmail.com), Rachel Turba (rturba@ucla.edu) and Rachel Meyer (rsmeyer@ucla.edu).
# Last Updated 11-18-2017
#
# The purpose of these script is to process raw fastq.gz files from an Illumina sequencing and generate summarized taxonomic assignment tables for multiple metabarcoding targets.
#
# This script is currently designed to run on UCLA's Hoffman2 cluster.  Please adjust the code to work with your computing resources. (e.g. module / path names to programs, submitting jobs for processing if you have a cluster, etc)
#
# This script runs in two phases, the first is a QC and dada2 seqeunce dereplication, denoising, mergeing (if reads are paired) and chimera detection.  The second phase runs bowtie2 and our blowtie2 modified blca run_scripts.
#
######################################

###Check that User has the correct set of primer files
if [[ ! -e ${FP} && ! -e ${RP} ]];
then
  echo "Using Default Primers"
elif [[ -e ${FP} && -e ${RP} ]];
then
  echo "Using User Defined Primers"
elif [[ -e ${FP} && ! -e ${RP} ]];
then
  echo "Using User Defined Primers"
  echo "Missing File For Reverse Primer"
  echo ""
  exit
elif [[ ! -e ${FP} && -e ${RP} ]];
then
  echo "Using User Defined Primers"
  echo "Missing File For Forward Primer!"
  echo ""
  exit
fi
######

#Check that user has all of the default flags set
echo $IN
echo $OUT
echo $DB
if [[ -e ${IN} && ! -z ${OUT} && -e ${DB} ]];
then
  echo "Required Arguments Given"
  echo ""
else
  echo "Required Arguments Missing:"
  echo "check that you included arguments or correct paths for -i -o -d -a"
  echo ""
  exit
fi

# location of the config and var files
source $DB/scripts/anacapa_vars.sh  # edit to change variables and parameters
source $DB/scripts/anacapa_config.sh # edit for proper configuration

###

################################
# Preprocessing .fastq files
################################
suffix1=R1_001.fastq.gz
suffix2=R2_001.fastq.gz
###################################
mkdir -p ${OUT}
mkdir -p ${OUT}/Run_info
mkdir -p ${OUT}/Run_info/run_logs
mkdir -p ${OUT}/QC
mkdir -p ${OUT}/QC/fastq

echo " "
date
echo " "
echo "Preprocessing: 1) Generate an md5sum file"  # user can check for file corruption
md5sum ${IN}/*fastq.gz > ${OUT}/Run_info/raw_fastq.md5sum
date
echo "Preprocessing: 2) Change file suffixes"
for str in `ls ${IN}/*_${suffix1}`
do
  str1=${str%*_${suffix1}}
  i=${str1#${IN}/}
  mod=${i//_/-}
  cp ${IN}/${i}_${suffix1} ${OUT}/QC/fastq/${mod}_1.fastq.gz
  cp ${IN}/${i}_${suffix2} ${OUT}/QC/fastq/${mod}_2.fastq.gz
done
date

###

################################
# QC the preprocessed .fastq files
#############################

echo "QC: 1) Run cutadapt to remove 5' sequencing adapters and 3' primers + sequencing adapters, sort for length, and quality."

# Generate cut adapt primer files -> merge reverse complemented primers with adapters for cutting 3'end sequencing past the end of the metabarcode region, and add cutadapt specific characters to primers and primer/adapter combos so that the appropriate ends of reads are trimmed
mkdir -p ${OUT}/Run_info/cutadapt_primers_and_adapters

echo " "
echo "Generating Primer and Primer + Adapter files for cutadapt steps."
cp ${DB}/adapters_and_PrimAdapt_rc/*_truseq_*_adapter.txt ${OUT}/Run_info/cutadapt_primers_and_adapters  # make a copy of the appropriate adapter file in your ourput directory
cp ${DB}/adapters_and_PrimAdapt_rc/*_nextera_*_adapter.txt ${OUT}/Run_info/cutadapt_primers_and_adapters  # make a copy of the appropriate adapter file in your ourput directory

cat ${OUT}/Run_info/cutadapt_primers_and_adapters/g_nextera_Forward_adapter.txt >> ${OUT}/Run_info/cutadapt_primers_and_adapters/g_Forward_adapter.txt
sed -i "s/F_adapt/F_adapt_nextera/" ${OUT}/Run_info/cutadapt_primers_and_adapters/g_Forward_adapter.txt
echo "" >> ${OUT}/Run_info/cutadapt_primers_and_adapters/g_Forward_adapter.txt
cat ${OUT}/Run_info/cutadapt_primers_and_adapters/g_truseq_Forward_adapter.txt >> ${OUT}/Run_info/cutadapt_primers_and_adapters/g_Forward_adapter.txt

cat ${OUT}/Run_info/cutadapt_primers_and_adapters/G_nextera_Reverse_adapter.txt >> ${OUT}/Run_info/cutadapt_primers_and_adapters/G_Reverse_adapter.txt


python ${DB}/scripts/anacapa_format_primers_cutadapt.py ${FP:=$FP_PATH} ${RP:=$RP_PATH} ${OUT}/Run_info/cutadapt_primers_and_adapters  # given your adapter and primer sets, make cutadapt readable fasta files for trimming adapter / primer reads

# now use the formated cutadapt primer file to trim fastq reads
mkdir -p ${OUT}/QC/cutadapt_fastq
mkdir -p ${OUT}/QC/cutadapt_fastq/untrimmed
mkdir -p ${OUT}/QC/cutadapt_fastq/primer_sort
mkdir -p ${OUT}/Run_info/cutadapt_out
###
for str in `ls ${OUT}/QC/fastq/*_1.fastq.gz`
do
  # first chop of the 5' adapter and 3' adapter and primer combo (reverse complemented)
  str1=${str%_*}
  j=${str1#${OUT}/QC/fastq/}
  echo " "
  echo ${j} "..."
  # this step removes all primers and adapters with the exception of the 5' forward and reverse primers.  These are needed in a later step to sort reads by primer set.  Leaving 3' primers and 5' or 3' adapters acn affect read merging and taxonomic assignment.
  # this cutadapt command allows a certain amount of error/missmatch (-e) between the query (seqeuncing read) and the primer and adapter.  It searches for and trims off all of the 5' forward adapter (-g) and the 3' reverse complement reverse primer / reverse complement reverse adapter (-a) or the 5' reverse adapter (-G) and the 3' reverse complement forward primer / reverse complement forward adapter (-A).  It processes read pairs, and results in two files one for each read pair.

  echo ${CUTADAPT} -e ${CTADE:=$ERROR_QC1} -f ${FILE_TYPE_QC1} -g ${F_ADAPT} -a ${Rrc_PRIM_ADAPT} -G ${R_ADAPT} -A ${Frc_PRIM_ADAPT} -o ${OUT}/QC/cutadapt_fastq/untrimmed/${j}_Paired_1.fastq -p ${OUT}/QC/cutadapt_fastq/untrimmed/${j}_Paired_2.fastq ${str1}_1.fastq.gz ${str1}_2.fastq.gz

  ${CUTADAPT} -e ${CTADE:=$ERROR_QC1} -f ${FILE_TYPE_QC1} -g ${F_ADAPT} -a ${Rrc_PRIM_ADAPT} -G ${R_ADAPT} -A ${Frc_PRIM_ADAPT} -o ${OUT}/QC/cutadapt_fastq/untrimmed/${j}_Paired_1.fastq -p ${OUT}/QC/cutadapt_fastq/untrimmed/${j}_Paired_2.fastq ${str1}_1.fastq.gz ${str1}_2.fastq.gz >> ${OUT}/Run_info/cutadapt_out/cutadapt-report.txt
  rm ${str1}_1.fastq.gz # remove intermediate files
  rm ${str1}_2.fastq.gz # remove intermediate files
  # stringent quality filter to get rid of the junky reads. It mostly chops the lowquality reads off of the ends. See the documentation for details. The default average quality score for retained bases is 35 and the minimum length is 100.  Any reads that do not meet that criteria are removed
  fastq_quality_trimmer -t ${QUALS:=$MIN_QUAL} -l ${MILEN:=$MIN_LEN}  -i ${OUT}/QC/cutadapt_fastq/untrimmed/${j}_Paired_1.fastq -o ${OUT}/QC/cutadapt_fastq/${j}_qcPaired_1.fastq -Q33 #trim pair one
  rm ${OUT}/QC/cutadapt_fastq/untrimmed/${j}_Paired_1.fastq # remove intermediate files
  fastq_quality_trimmer -t ${QUALS:=$MIN_QUAL} -l ${MILEN:=$MIN_LEN}  -i ${OUT}/QC/cutadapt_fastq/untrimmed/${j}_Paired_2.fastq -o ${OUT}/QC/cutadapt_fastq/${j}_qcPaired_2.fastq -Q33 #trim pair 2
  rm ${OUT}/QC/cutadapt_fastq/untrimmed/${j}_Paired_2.fastq # remove intermediate files
  # sort by metabarcode but run additional trimming.  It makes a differnce in merging reads in dada2.  Trimming varies based on seqeuncing platform.
  echo "forward..."
  # use cut adapt to search 5' end of forward reads for forward primers.  These are then sorted by primer name.  We do an additional trimming step analagous to the trimming step in the dada2 tutorial.  Because these a forward reads an tend to be higher quality we only trim  20 bp from the end by default for the MiSeq (longer Reads). Users can modify all parameters in the vars file.
  ${CUTADAPT} -e ${PCTADE:=$ERROR_PS} -f ${FILE_TYPE_PS} -g ${F_PRIM}  -u -${FETRIM:=$MS_F_TRIM} -o ${OUT}/QC/cutadapt_fastq/primer_sort/{name}_${j}_Paired_1.fastq  ${OUT}/QC/cutadapt_fastq/${j}_qcPaired_1.fastq >> ${OUT}/Run_info/cutadapt_out/cutadapt-report.txt
  echo "check"
  echo "reverse..."
  # use cut adapt to search 5' end of reverse reads for reverse primers.  These are then sorted by primer name.  We do an additional trimming step analagous to the trimming step in the dada2 tutorial.  Because these a reverse reads an tend to be lower quality we only trim  50 bp from the end by default for the MiSeq (longer Reads). Users can modify all parameters in the vars file.
  ${CUTADAPT} -e ${PCTADE:=$ERROR_PS} -f ${FILE_TYPE_PS} -g ${R_PRIM}  -u -${RETRIM:=$MS_R_TRIM} -o ${OUT}/QC/cutadapt_fastq/primer_sort/{name}_${j}_Paired_2.fastq   ${OUT}/QC/cutadapt_fastq/${j}_qcPaired_2.fastq >> ${OUT}/Run_info/cutadapt_out/cutadapt-report.txt
  echo "check"
  rm ${OUT}/QC/cutadapt_fastq/${j}_qcPaired_1.fastq # remove intermediate files
  rm ${OUT}/QC/cutadapt_fastq/${j}_qcPaired_2.fastq # remove intermediate files
done
date
###

###############################
# Make sure unassembled reads are still paired
###############################
mkdir -p ${OUT}/Run_info
mkdir -p ${OUT}/Run_info/run_scripts

metabarcodes="$( ls -l | grep -o '>.*' ${FP:=$FP_PATH}  | cut -c 2- | tr '\n' ' ' )"

echo ${metabarcodes}

echo " "
echo "Checking that Paired reads are still paired:"
for j in ${metabarcodes}
do
  echo " "
  echo ${j} "..."
  # make directories for the soon to be sorted reads.  Reads will be sorted by primer and then by paired or unpaired read status.
  shopt -s nullglob
  files=(${OUT}/QC/cutadapt_fastq/primer_sort/${j}*)
  # now check the size of the array
  if (( ${#files[@]} == 0 )); then
      echo "No reads matched the ${j} 3' primer sequences"
  else
    mkdir -p ${OUT}/${j}
    mkdir -p ${OUT}/${j}/${j}_sort_by_read_type
    mkdir -p ${OUT}/${j}/${j}_sort_by_read_type/paired/
    mkdir -p ${OUT}/${j}/${j}_sort_by_read_type/unpaired_F/
    mkdir -p ${OUT}/${j}/${j}_sort_by_read_type/unpaired_R/
    
    echo ${OUT}/QC/cutadapt_fastq/primer_sort/${j}_*_Paired_1.fastq
    for st in `ls ${OUT}/QC/cutadapt_fastq/primer_sort/${j}_*_Paired_1.fastq`
	  do
      st2=${st%*_Paired_1.fastq}
      k=${st2#${OUT}/QC/cutadapt_fastq/primer_sort/}
      # For each sample and each metabarcode, this python script checks to see if the forward and reverse files have read pairs, or singleton F or R reads.  Reads are then sorted into the directories generated above.
      python ${DB}/scripts/check_paired.py ${OUT}/QC/cutadapt_fastq/primer_sort/${k}_Paired_1.fastq ${OUT}/QC/cutadapt_fastq/primer_sort/${k}_Paired_2.fastq ${OUT}/${j}/${j}_sort_by_read_type/paired ${OUT}/${j}/${j}_sort_by_read_type/unpaired_F/ ${OUT}/${j}/${j}_sort_by_read_type/unpaired_R/
      echo ${k} "...check!"
     done
     date
   fi
done

# rm -r ${OUT}/QC/fastq # remove intermediate directories
# rm -r ${OUT}/QC/cutadapt_fastq # remove intermediate directories
# rm -r ${OUT}/QC

###############################
# Make sure unassembled reads are still paired and submit dada2 jobs
###############################

echo " "
echo "Process metabarcode reads for with dada2"
for j in `ls ${OUT}/`
do
 if [[ "${j}" != "QC" && "${j}" != "Run_info" ]]; # ignore non-metabarcode folders...
 then
    #make folders for the metabarcode specific output of dada2 and bowtie2
  echo ""
  echo "${j}"
	mkdir -p ${OUT}/${j}/${j}dada2_out

  echo "Running Dada2"
  printf "#!/bin/bash\n ${RUNNER} ${DB}/scripts/run_dada2.sh -o ${OUT} -d ${DB} -m ${j} -t paired -b ${MINTIMES_ASV:=$MIN_ASV_ABUNDANCE} \n" > ${OUT}/Run_info/run_scripts/${j}_dada2_paired_job.sh
  printf "#!/bin/bash\n ${RUNNER} ${DB}/scripts/run_dada2.sh -o ${OUT} -d ${DB} -m ${j} -t forward -b ${MINTIMES_ASV:=$MIN_ASV_ABUNDANCE} \n" > ${OUT}/Run_info/run_scripts/${j}_dada2_F_job.sh
  printf "#!/bin/bash\n ${RUNNER} ${DB}/scripts/run_dada2.sh -o ${OUT} -d ${DB} -m ${j} -t reverse -b ${MINTIMES_ASV:=$MIN_ASV_ABUNDANCE} \n" > ${OUT}/Run_info/run_scripts/${j}_dada2_R_job.sh
  ${RUNNER} ${OUT}/Run_info/run_scripts/${j}_dada2_paired_job.sh
  date
  ${RUNNER} ${OUT}/Run_info/run_scripts/${j}_dada2_F_job.sh
  date
  ${RUNNER} ${OUT}/Run_info/run_scripts/${j}_dada2_R_job.sh
  date
 fi
done

echo "All done"