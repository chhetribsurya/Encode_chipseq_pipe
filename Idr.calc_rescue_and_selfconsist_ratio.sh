#!/usr/bin/bash

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
if [ -z "$LSB_JOBID" ]; then echo "Please run on a compute node. Exiting"; exit 1; fi

if [ -z "$TEMP_DIR" ]; then TEMP_DIR=$(get_temp_dir); fi
if [ ! -d "$TEMP_DIR" ]; then
    mkdir -p $TEMP_DIR
    if [ $? -ne "0" ]; then echo "Could not create output dir: $TEMP_DIR. Exiting"; exit 1; fi
fi

### This script will be passed with 4 narrowpeak zipped files - to calculate rescue and self consistency ratio:
#Nt = Best no. of peaks passing IDR threshold by comparing true replicates (# conservative set)
#Np = No. of peaks passing IDR threshold by comparing pooled pseudo-replicates
#N1 and N2 = No. of peaks passing IDR threshold by comparing self-pseudoReplicates for Rep1 and Rep2 respectively (# gives self consistent IDR peaks)

#USAGE: ./Idr.calc_rescue_and_selfconsist_ratio.sh <True_IDR_peaks> <Pooled_IDR_peaks> <Rep1_self-pseudorep_IDR_peaks> <Rep2_self-pseudorep_IDR_peaks> <output_file>

if [[ "$#" -ne 5 ]]; then
    echo -e "\nIncorrect no. of parameters, exitting....\nUSAGE: ./Idr.calc_rescue_and_selfconsist_ratio.sh <True_IDR_peaks> <Pooled_IDR_peaks> <Rep1_self-pseudorep_IDR_peaks> <Rep2_self-pseudorep_IDR_peaks> <output_file> \n"
    exit 1
else
    echo -e "\nRunScript: ${0}  ...,  
                Output file: ${5}  ...\n"
fi

export TRUE_IDR_PEAK_FILE=$1
export POOLED_IDR_PEAK_FILE=$2
export R1_SELF_PSEUDOREP_IDR_PEAK_FILE=$3
export R2_SELF_PSEUDOREP_IDR_PEAK_FILE=$4
export OUTPUT_FILE=$5

if [[ -f $OUTPUT_FILE ]];then
    rm $OUTPUT_FILE # to avoid appending on same file
fi

Nt_TRUE_IDR_PEAKS=$(zcat $TRUE_IDR_PEAK_FILE | wc -l)
Np_POOLED_IDR_PEAKS=$(zcat $POOLED_IDR_PEAK_FILE | wc -l)
N1_R1_SELF_PSEUDOREP_IDR_PEAKS=$(zcat $R1_SELF_PSEUDOREP_IDR_PEAK_FILE | wc -l)
N2_R2_SELF_PSEUDOREP_IDR_PEAKS=$(zcat $R2_SELF_PSEUDOREP_IDR_PEAK_FILE | wc -l)
echo "Number of True IDR peaks: $Nt_TRUE_IDR_PEAKS" >> $OUTPUT_FILE
echo "Number of Pooled IDR peaks: $Np_POOLED_IDR_PEAKS" >> $OUTPUT_FILE
echo "Number of Rep1 pooled IDR peaks: $N1_R1_SELF_PSEUDOREP_IDR_PEAKS" >> $OUTPUT_FILE
echo "Number of Rep2 pooled IDR peaks: $N2_R2_SELF_PSEUDOREP_IDR_PEAKS" >> $OUTPUT_FILE

Rescue_num="$Nt_TRUE_IDR_PEAKS"
Rescue_num+=" $Np_POOLED_IDR_PEAKS"
Rescue_max=$(echo $Rescue_num | xargs -n1 | sort -n | tail -1)
Rescue_min=$(echo $Rescue_num | xargs -n1 | sort -n | head -1)
#Rescue_Ratio=$(echo "scale=2; $Rescue_max/$Rescue_min" | bc)
Rescue_Ratio=$(awk "BEGIN {printf \"%.2f\",${Rescue_max}/${Rescue_min}}")
echo -e  "\n\nRESCUE RATIO = $Rescue_Ratio" >> $OUTPUT_FILE


SelfC_num="$N1_R1_SELF_PSEUDOREP_IDR_PEAKS"
SelfC_num+=" $N2_R2_SELF_PSEUDOREP_IDR_PEAKS"
SelfC_max=$(echo $SelfC_num | xargs -n1 | sort -n | tail -1)
SelfC_min=$(echo $SelfC_num | xargs -n1 | sort -n | head -1)
#Self_consistency_Ratio=$(echo "scale=2; $SelfC_max/$SelfC_min" | bc)
Self_consistency_Ratio=$(awk "BEGIN {printf \"%.2f\",${SelfC_max}/${SelfC_min}}")
echo -e  "SELF-CONSISTENCY RATIO = $Self_consistency_Ratio" >> $OUTPUT_FILE

if [[ ( $Rescue_Ratio > 2 ) && ( $Self_consistency_Ratio > 2 ) ]];then
    echo -e "\n\nConclusion:\n REPRODUCIBILITY FAIL (Both RR and SR > 2)" >> $OUTPUT_FILE

elif [[ ( $Rescue_Ratio > 2 ) || ( $Self_consistency_Ratio > 2 ) ]];then
    echo -e "\n\nConclusion:\n REPRODUCIBILITY BORDERLINE (Either RR or SR > 2)" >> $OUTPUT_FILE
else
    echo -e "\n\nConclusion:\n IDR QC PASSED (Both RR and SR < 2)" >> $OUTPUT_FILE
fi
# Making other copy of output file to keep all under same hood:
cp $OUTPUT_FILE $RATIO_OUTPUT_DIR

log_msg info "IDR Rescue and Self consistency ratio calculation complete."
