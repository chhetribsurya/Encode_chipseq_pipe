#!/bin/bash

##############################################################
### Run Motif Analysis                                     ###
### Uses Meme scripts and data base                        ###
### Dependent on SPP peak calling finishing sucessfully   ###
### Output is a directory with eps images and html display ###
##############################################################
#SOURCE_DIR="/opt/HAIB/myerslab/etc"
SOURCE_DIR="/gpfs/gpfs2/software/HAIB/myerslab/etc"

### Mandatory sourcing of bashrc for necessary environment variables. ###
if [ -e $SOURCE_DIR/bashrc ]; then 
    . $SOURCE_DIR/bashrc
else echo "[fatal] - Could not find myerslab bashrc file. Exiting"; exit 1; fi

### Mandatory sourcing of functions to get helper functions (like call_cmd). ###
if [ -e $SOURCE_DIR/functions ]; then
    . $SOURCE_DIR/functions
else echo "[fatal] - Could not find functions file. Exiting"; exit 1; fi

### Verify we are not running on the head node. ###
if [ -z "$LSB_JOBID" ]; then log_msg fatal "Please run on a compute node. Exiting"; exit 1; fi

#Set the variables:
IDR_PEAK_FILE=$1
R1=$2
R2=$3

LIBRARY=IDR_${R1}_${R2}

### Check if the IDR passed peak file exists, and is more than 100 peaks; but ideally, at least 500 peaks would be needed to call the motifs from those sequences

if [[ -f $IDR_PEAK_FILE ]];then
	wc_num=$(zcat $IDR_PEAK_FILE | wc -l)
	peak_num=$(echo $wc_num | cut -d " " -f 1)
	if [[ $peak_num -lt 100 ]]; then
		echo -e "\nToo few IDR passed peaks, so the model cannot be built, skipping motif finding\n"
		PEAK_CALLING_FAILED="TRUE"
	fi
fi

TEMP_IDR_DIR=$BASE_DIR/IDR_${R1}_${R2}/meme_chip/temp
MOTIF_OUTPUT=$BASE_DIR/IDR_${R1}_${R2}/meme_chip

if [[ ! -d $TEMP_IDR_DIR ]]; then
	mkdir -p $TEMP_IDR_DIR
fi

if [[ ! -d $MOTIF_OUTPUT ]]; then
	mkdir -p $MOTIF_OUTPUT
fi


if [ -z $PEAK_CALLING_FAILED ]; then

    SORTED_FILE=$TEMP_IDR_DIR/${LIBRARY}.sorted_file
    CENTERED_FILE=$TEMP_IDR_DIR/${LIBRARY}_summit_500bp.bed
    CENTERED_FASTA_FILE=$TEMP_IDR_DIR/${LIBRARY}_summit_500bp_file.fasta
    BG_FILE=$TEMP_IDR_DIR/${LIBRARY}_background.bed
    BG_FASTA_FILE=$TEMP_IDR_DIR/${LIBRARY}_background_file.fasta
    BG_MARKOV_MODEL_FILE=$TEMP_IDR_DIR/${LIBRARY}_background_file.model
    SUMMARY_DESC_FILE=$MOTIF_OUTPUT/summary_file
    touch $SUMMARY_DESC_FILE

    ### sorting of the peaks based on signal float value:
    #zcat $IDR_PEAK_FILE | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5,$7,$10}' | sort -k 5,5gr > $SORTED_FILE
    echo -e "\nsorting of IDR peak file completed...\n"   

    ### generating bed file of 500bp regions centered on peak-summits: 
    #awk 'BEGIN{OFS="\t"} {chromStart=$2; summit=$6; midPos=chromStart+summit; print $1, midPos-250, midPos+250;}' $SORTED_FILE > $CENTERED_FILE
    echo -e "\nCentering of IDR peak file 250 upstream and 250 downstream completed...\n"   
    
    ### Generate 2X null sequences or random sequences with matched GC content, repeat fraction with user input sequence length:
    #python $NULL_GENERATE_SCRIPT $NULL_PARAMETERS -o $BG_FILE $CENTERED_FILE hg19 $NULL_HG19_INDICES
    echo -e  "\nGeneration of randomic genomic regions matching GC% and length completed...\n"   

    ### fetching the DNA sequences of peak summit regions using hg19 male fasta reference seq:
    #$BEDTOOLS_PATH/fastaFromBed -fi $GENOME -bed $CENTERED_FILE -fo $CENTERED_FASTA_FILE
    echo -e "\nFasta extraction of summit centered file completed...\n"   

    ### fetching the DNA sequences of randomly generated null sequences with matched GC content and repeat fraction using hg19 ref seq:
    #$BEDTOOLS_PATH/fastaFromBed -fi $GENOME -bed $BG_FILE -fo $BG_FASTA_FILE
    echo -e "\nFasta extraction of background file completed...\n"   

    ### create fasta markov model as a background file:
    $MEME_SUITE_PATH/fasta-get-markov -m 1 $BG_FASTA_FILE $BG_MARKOV_MODEL_FILE
    echo -e "\nMarkov model generation for background file completed...\n"   

    ## Run meme-chip with parallel cores and nonrandom features since the input fasta is sorted with decreasing confidence of the peaks or signal float value:
    echo -e "\nStarting the meme-chip operation using the centered fasta file and bg markov model file generated in prior steps...\n"
    perl $MEME_SUITE_PATH/meme-chip -oc $MOTIF_OUTPUT -index-name meme_combined.html -db $MOTIF_DB_PATH/JASPAR/JASPAR_CORE_2016_vertebrates.meme -db $MOTIF_DB_PATH/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme -dna -bfile $BG_MARKOV_MODEL_FILE -norand -nmeme 500 -fdesc $SUMMARY_DESC_FILE -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 5 -meme-maxsize 100000 -dreme-e 0.05 -centrimo-local -centrimo-score 5 -centrimo-ethresh 10 ${CENTERED_FASTA_FILE} -spamo-skip
    echo -e "\nMeme-chip analysis completed....\n"

    for NUM in 1 2 3; do
        if [ -e $MOTIF_OUTPUT/meme_out/logo${NUM}.png ]; then
            CMD="/usr/bin/convert $MOTIF_OUTPUT/meme_out/logo${NUM}.png $MOTIF_OUTPUT/meme_out/logo${NUM}.jpg"
            OUTPUT_FILES="$MOTIF_OUTPUT/meme_out/logo${NUM}.jpg"
            run_cmd "$CMD" "$OUTPUT_FILES"
        fi
    done
fi

echo -e "\nMotif analysis completed\n"
