This repo contains a modified version of Anacapa toolkit (specifically just QC) used by the CALeDNA project

```
conda create -n crux27 -c bioconda python=2.7 r-essentials r-base cutadapt fastx_toolkit
```

### Preparing the anacapa_config.sh file

Before running the __Anacapa__ toolkit you need to double check the anacapa_config.sh file and update the appropriate paths. For running local mode on a personal machine or virutual box set CUTADAPT ="cutadapt",and  MUSCLE="muscle"; replace all other values to "". Double check that all dependencies work in the terminal. This is the key for success.


### Running _anacapa_QC_dada2.sh_
```
/bin/bash ~/Anacapa_db/anacapa_QC_dada2.sh -h

<<< Anacapa: Sequence QC and ASV Parsing >>>

The purpose of these script is to process raw fastq or fastq.gz files from an Illumina HiSeq or MiSeq.  It removes 3' and 5' sequencing artifacts and 5' metabarcode primers (cutadapt), removes low quality base pairs and short reads (fastX-toolkit), sorts reads by 3' metabarcode primers prior to trimming (cutadapt), and uses dada2 to denoise, dereplicate, merge and remove chimeric reads

	For successful implementation
		1. Make sure you have all of the dependencies and correct paths in the anacapa_config.sh file
		2. Add the Metabarcode locus specific CRUX reference libraries to the Anacapa_db folder
		3. All parameters can be modified using the arguments below.  Alternatively, all parameters can be altered in the anacapa_vars.sh folder

Arguments:
- Required for either mode:
	-i	path to .fastq.gz files, if files are already compressed use flag -g (see below)
	-o	path to output directory
	-d	path to Anacapa_db

 - Optional:
 	-f	path to file with forward primers in fasta format
    		e.g.	 >16s
    			     GTGYCAGCMGCCGCGGTAA
			         >18S
			         GTACACACCGCCCGTC
	-r	path to file with reverse primers in fasta format
    		e.g. 	 >16s
    			     GGACTACNVGGGTWTCTAAT
    			     >18S
			         TGATCCTTCTGCAGGTTCACCTAC
	-g	If .fastq read are uncompressed: -g (no argument need)
	-c	To modify the allowed cutadapt error for 3' adapter and 5' primer adapter trimming: 0.0 to 1.0 (default 0.3)
	-p	To modify the allowed cutadapt error 3' primer sorting and trimming: 0.0 to 1.0 (default 0.3)
	-q	To modify the minimum quality score allowed: 0 - 40 (default 35)
	-m	To modify the minimum length after quality trimming: 0 - 300 (default 100)
	-x	To modify the additional 5' trimming of forward reads: 0 - 300 (default HiSeq 10, default MiSeq 20)
	-y	To modify the additional 5' trimming of reverse reads: 0 - 300 (default HiSeq 25, default MiSeq 50)
	-b	To modify the number of occurrences required to keep an ASV: 0 - any integer (default 0)
	-j	Multithreading (True/False) in dada2. Multithreading is turned on by default. If user wished to process a 		  single sample turn multithreading to FALSE.  
 	-k	Path to file with alternate HPC job submission parameters:  
		default file = ~/Anacapa_db/scripts/Hoffman2_HPC_header.sh
		modifiable template file = ~/Anacapa_db/scripts/anacapa_qsub_templates.sh


 - Other:
	-h	Shows program usage then quits


```
__NOTE__: Script does not check that the user provided optional argument are within the correct values range. __Please use care when entering arguments.__

#### The output of the anacapa_QC_dada2.sh is as follows:
  * Output directory based on user designated name
    * Subdirectories for each metabarcode target
	   * Subdirectory for clean sorted data called "<metabarcode>\_sort_by_read_type"
    * Subdirectory for run information called "Run_info"
      * Subdirectory for cutadapt logs called "cutadapt_out"
      * Subdirectory for cutadapt primers and adapters called "cutadapt_primers_and_adapters"
      * the .md5sum file for the raw sequences
      * Subdirectory for general Anacapa run logs called  "run_logs"
      * Subdirectory Anacapa run scripts called "runscripts"
