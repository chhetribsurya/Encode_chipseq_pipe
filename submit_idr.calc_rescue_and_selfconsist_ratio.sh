#!/bin/bash

#source /gpfs/gpfs1/home/schhetri/python/anaconda_python_version3.sh
#source activate python3

#BASE_DIR=/gpfs/gpfs1/home/snewberry/interrupt/encode/spp_and_idr/idr_testing

#BSUB_OPTIONS="-n 8 -R span[hosts=1]"
#if [ $# -ne 3 ]; then echo "Usage $0: <Rep1 SL#> <Rep2 SL#> <IDR directory name>"; exit 1; fi

#R1=$1
#R2=$2
#IDR_DIR=$3

#So, to overcome the issue of exporting the bash array variable, so directly: (Read lines in file into an bash array)
#readarray LIB_REP1 < ./rep1_lib.txt
#readarray LIB_REP2 < ./rep2_lib.txt
#readarray LIB_CONTROL1 < ./control1_lib.txt
#readarray LIB_CONTROL2 < ./control2_lib.txt
#
#LIB_REP1=( $(echo ${LIB_REP1[@]} | tr " " "\n" | tr "\n" " "))
#LIB_REP2=( $(echo ${LIB_REP2[@]} | tr " " "\n" | tr "\n" " "))
#LIB_CONTROL1=( $(echo ${LIB_CONTROL1[@]} | tr " " "\n" | tr "\n" " "))
#LIB_CONTROL2=( $(echo ${LIB_CONTROL2[@]} | tr " " "\n" | tr "\n" " "))

export RUN_PATH=`pwd`
export RATIO_OUTPUT_DIR=$BASE_DIR/Ratio

if [[ ! -d $RATIO_OUTPUT_DIR ]];then
    mkdir $RATIO_OUTPUT_DIR
fi

#So, to overcome the issue of exporting the bash array variable, so directly: (Read lines in file into an bash array)
readarray LIB_REP1 < $REP1_FILE_NAME
readarray LIB_REP2 < $REP2_FILE_NAME
readarray LIB_CONTROL1 < $CONTROL1_FILE_NAME
readarray LIB_CONTROL2 < $CONTROL2_FILE_NAME
readarray ALL_LIB < $ALL_TF_FILE_NAME

#Reading or geting unique values from an bash array, and then note it's reconverted into bash array by placing array bracket ( ):
LIB_REP1=( $(echo ${LIB_REP1[@]} | tr " " "\n" | tr "\n" " "))
LIB_REP2=( $(echo ${LIB_REP2[@]} | tr " " "\n" | tr "\n" " "))
LIB_CONTROL1=( $(echo ${LIB_CONTROL1[@]} | tr " " "\n" | tr "\n" " "))
LIB_CONTROL2=( $(echo ${LIB_CONTROL2[@]} | tr " " "\n" | tr "\n" " "))
UNIQ_ALL_LIB=( $(echo ${ALL_LIB[@]} | tr " " "\n" | sort -u | tr "\n" " ") )


for i in "${!LIB_REP1[@]}"; do

    R1=${LIB_REP1[$i]}
    R2=${LIB_REP2[$i]}
	IDR_DIR=$BASE_DIR/IDR_${R1}_${R2}
	
    # True Rep IDR
	Nt_A=$IDR_DIR/$R1.filt.nodup.srt.PE2SE.narrowPeak.gz 
    Nt_B=$IDR_DIR/$R2.filt.nodup.srt.PE2SE.narrowPeak.gz 
    #Nt_C=$IDR_DIR/${R1}_${R2}.Rep0.narrowPeak.gz 

	# Rep 1 self-pseudoreps (Rep1 psuedorep1 vs Rep1 pseudorep2)
	N1_A=$IDR_DIR/$R1.filt.nodup.PE2SE.pr1.narrowPeak.gz 
    N1_B=$IDR_DIR/$R1.filt.nodup.PE2SE.pr2.narrowPeak.gz 
    #N1_C=$IDR_DIR/$R1.filt.nodup.srt.PE2SE.narrowPeak.gz 

	# Rep 2 self-pseudoreps (Rep2 pseudorep1 vs Rep2 pseudorep2) 
	N2_A=$IDR_DIR/$R2.filt.nodup.PE2SE.pr1.narrowPeak.gz 
    N2_B=$IDR_DIR/$R2.filt.nodup.PE2SE.pr2.narrowPeak.gz 
    #N2_C=$IDR_DIR/$R2.filt.nodup.srt.PE2SE.narrowPeak.gz 

	# Pooled psuedoreps ( Pooled pseudorep1 (PR1 of Rep1 + PR1 of Rep2) vs Pooled true tag align file(Rep1_TA_file + Rep2_TA_file)
	Np_A=$IDR_DIR/${R1}_${R2}.Rep0.pr1.narrowPeak.gz 
    Np_B=$IDR_DIR/${R1}_${R2}.Rep0.pr2.narrowPeak.gz 
    #Np_C=$IDR_DIR/${R1}_${R2}.Rep0.narrowPeak.gz signal.value $IDR_DIR

    Nt_REP1_VS_REP2=$(basename $Nt_A .narrowPeak.gz)_VS_$(basename $Nt_B .narrowPeak.gz)
    Np_REP1_VS_REP2=$(basename $Np_A .narrowPeak.gz)_VS_$(basename $Np_B .narrowPeak.gz)
    N1_REP1_VS_REP2=$(basename $N1_A .narrowPeak.gz)_VS_$(basename $N1_B .narrowPeak.gz)
    N2_REP1_VS_REP2=$(basename $N2_A .narrowPeak.gz)_VS_$(basename $N2_B .narrowPeak.gz)
    
    Nt_IDR_THRESHOLD_FILE=$IDR_DIR/${Nt_REP1_VS_REP2}.IDR0.02.narrowPeak.gz
    Np_IDR_THRESHOLD_FILE=$IDR_DIR/${Np_REP1_VS_REP2}.IDR0.02.narrowPeak.gz
    N1_IDR_THRESHOLD_FILE=$IDR_DIR/${N1_REP1_VS_REP2}.IDR0.02.narrowPeak.gz
    N2_IDR_THRESHOLD_FILE=$IDR_DIR/${N2_REP1_VS_REP2}.IDR0.02.narrowPeak.gz

    echo -e "\nJust a Reminder ....."
	echo -e "dir for logfile : $IDR_DIR/idr.${REP1}_${REP2}.truerep.out" 
	echo -e "check all idr peaks agnostic of cutoff: $IDR_DIR/${REP1}.filt.nodup.srt.PE2SE_VS_${REP2}.filt.nodup.srt.PE2SE.IDR0.narrowPeak.gz"
	echo -e "\ncheck idr passed peaks with cutoff 0.02 (Final peaks for analysis) : $IDR_DIR/${REP1}.filt.nodup.srt.PE2SE_VS_${REP2}.filt.nodup.srt.PE2SE.IDR0.02.filt.narrowPeak.gz\n"
	
    # Calculate rescue ration and self consistency ratio:
	$RUN_PATH/Idr.calc_rescue_and_selfconsist_ratio.sh $Nt_IDR_THRESHOLD_FILE $Np_IDR_THRESHOLD_FILE $N1_IDR_THRESHOLD_FILE $N2_IDR_THRESHOLD_FILE $IDR_DIR/${R1}_${R2}_rescue_and_self_consistency_ratio.txt

done

