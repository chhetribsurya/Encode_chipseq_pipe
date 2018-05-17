#BASE_DIR=/gpfs/gpfs1/home/snewberry/interrupt/encode/spp_and_idr/idr_testing
#OUTPUT_DIR=$BASE_DIR/Analysis/$LIBRARY
#TEMP_DIR=$BASE_DIR/TEMP
#REFERENCE_DIR=$BASE_DIR/../reference

SOURCE_DIR="/gpfs/gpfs2/software/HAIB/myerslab/etc"
#SOURCE_DIR="/opt/HAIB/myerslab/etc"

### Mandatory sourcing of bashrc for necessary environment variables. ###
if [ -e $SOURCE_DIR/bashrc ]; then 
    . $SOURCE_DIR/bashrc
else echo "[fatal] - Could not find myerslab bashrc file. Exiting"; exit 1; fi

### Mandatory sourcing of functions to get helper functions (like call_cmd). ###
if [ -e $SOURCE_DIR/functions ]; then
    . $SOURCE_DIR/functions
else echo "[fatal] - Could not find functions file. Exiting"; exit 1; fi

#export R_LIBS_SITE=/gpfs/gpfs1/software/R-site-packages
export R_LIBS_SITE=/gpfs/gpfs2/software/c7R-libs

### Verify we are not running on the head node. ###
if [ -z "$LSB_JOBID" ]; then echo "Please run on a compute node. Exiting"; exit 1; fi

if [ -z "$LIBRARY" ];            then empty_param_quit "LIBRARY"; fi
if [ -z "$GENOME" ];             then empty_param_quit "GENOME"; fi

### Verify library directory exists. ###
if [ -z "$LIBRARY_DIR" ]; then LIBRARY_DIR=$(get_library_dir $LIBRARY); fi
if [ ! -d "$LIBRARY_DIR" ]; then echo "Library directory does not exist: $LIBRARY_DIR. Exiting"; exit 1; fi

if [ -z "$OUTPUT_DIR" ]; then OUTPUT_DIR=$(get_output_dir $LIBRARY);  fi
if [ ! -d "$OUTPUT_DIR" ]; then 
    mkdir -p $OUTPUT_DIR
    if [ $? -ne "0" ]; then echo "Could not create output dir: $OUTPUT_DIR. Exiting"; exit 1; fi
fi

### Set up temp directory variable ###
if [ -z "$TEMP_DIR" ]; then TEMP_DIR=$(get_temp_dir); fi
if [ ! -d "$TEMP_DIR" ]; then
    mkdir -p $TEMP_DIR
    if [ $? -ne "0" ]; then echo "Could not create output dir: $TEMP_DIR. Exiting"; exit 1; fi
fi


export LOGFILE_NAME=$OUTPUT_DIR/$LIBRARY.log

#NTHREADS=8
NTHREADS=4

log_msg info "Beginning Fastq to TagAlign"
log_msg info "    LIBRARY:              $LIBRARY"
log_msg info "    GENOME:               $GENOME"
log_msg info "    LIBRARY_DIR:          $LIBRARY_DIR"
log_msg info "    OUTPUT_DIR:           $OUTPUT_DIR"

export TEMP_DIR=$OUTPUT_DIR/temp_dir
if [ ! -d "$TEMP_DIR" ]; then
    mkdir -p $TEMP_DIR
    if [ $? -ne "0" ]; then echo "Could not create output dir: $TEMP_DIR. Exiting"; exit 1; fi
fi

log_msg info "Consolidating fastq files..."
# Set up output files
#FASTQ_FILE_1=$TEMP_DIR/$LIBRARY.fastq.gz
#SAI_FILE_1=$TEMP_DIR/$LIBRARY.sai
#RAW_BAM_PREFIX=$LIBRARY.raw.srt
#RAW_BAM_FILE=$OUTPUT_DIR/$RAW_BAM_PREFIX.bam
#RAW_BAM_FILE_MAPSTATS=$OUTPUT_DIR/$RAW_BAM_PREFIX.flagstat.qc

# Set up output files
FASTQ_FILE_1=$TEMP_DIR/${LIBRARY}_1.fastq.gz
FASTQ_FILE_2=$TEMP_DIR/${LIBRARY}_2.fastq.gz
SAI_FILE_1=$TEMP_DIR/${LIBRARY}_1.sai
SAI_FILE_2=$TEMP_DIR/${LIBRARY}_2.sai
RAW_SAM_FILE=$TEMP_DIR/${LIBRARY}.raw.sam.gz
BAM_FILE=$TEMP_DIR/${LIBRARY}.raw.bam
RAW_BAM_PREFIX=${LIBRARY}.raw.srt
RAW_BAM_FILE=$OUTPUT_DIR/${RAW_BAM_PREFIX}.bam # to be stored
RAW_BAM_FILE_MAPSTATS=$OUTPUT_DIR/${RAW_BAM_PREFIX}.flagstat.qc # QC file

# Remove read pairs with bad CIGAR strings and sort by position 
BADCIGAR_FILE="${TEMP_DIR}/badReads${LIBRARY}.tmp"  #RAW_BAM_FILE_MAPSTATS="${RAW_BAM_PREFIX}.flagstat.qc" # QC File

CMD="cat $LIBRARY_DIR/*_1.fastq.gz > $FASTQ_FILE_1"
run_cmd "$CMD" "$FASTQ_FILE_1"

CMD="cat $LIBRARY_DIR/*_2.fastq.gz > $FASTQ_FILE_2"
run_cmd "$CMD" "$FASTQ_FILE_2"

BWA_SOFTWARE=$(get_software_dir bwa-0.7.12/bwa)
log_msg info "Aligning with bwa aln..."
$BWA_SOFTWARE aln -q 5 -l 32 -k 2 -t ${NTHREADS} $REFERENCE_DIR/$GENOME.fa ${FASTQ_FILE_1} > ${SAI_FILE_1}
$BWA_SOFTWARE aln -q 5 -l 32 -k 2 -t ${NTHREADS} $REFERENCE_DIR/$GENOME.fa ${FASTQ_FILE_2} > ${SAI_FILE_2}
$BWA_SOFTWARE sampe $REFERENCE_DIR/$GENOME.fa ${SAI_FILE_1} ${SAI_FILE_2} ${FASTQ_FILE_1} ${FASTQ_FILE_2} | gzip -c > ${RAW_SAM_FILE}
rm ${SAI_FILE_1} ${SAI_FILE_2}

# Find bad CIGAR read names                             
zcat ${RAW_SAM_FILE} | awk 'BEGIN {FS="\t" ; OFS="\t"} ! /^@/ && $6!="*" { cigar=$6; gsub("[0-9]+D","",cigar); n = split(cigar,vals,"[A-Z]"); s = 0; for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10) ; if (s!=seqlen) print $1"\t"; }' | sort | uniq > ${BADCIGAR_FILE}

#Remove bad CIGAR read pairs
if [[ $(cat ${BADCIGAR_FILE} | wc -l) -gt 0 ]]
then
    zcat ${RAW_SAM_FILE} | grep -v -F -f ${BADCIGAR_FILE} | $SAMTOOLS_PATH/samtools view -Su - | $SAMTOOLS_PATH/samtools sort - $OUPUT_DIR/${RAW_BAM_PREFIX}
else
    $SAMTOOLS_PATH/samtools view -Su ${RAW_SAM_FILE} | $SAMTOOLS_PATH/samtools sort - $OUTPUT_DIR/${RAW_BAM_PREFIX}
fi
rm ${BADCIGAR_FILE} ${RAW_SAM_FILE}

# Stats for final ra bam file
$SAMTOOLS_PATH/samtools flagstat ${RAW_BAM_FILE} > ${RAW_BAM_FILE_MAPSTATS}

# ===================================
# Remove  unmapped, mate unmapped
# not primary alignment, reads failing platform
# Remove low MAPQ reads
# Only keep properly paired reads
# Obtain name sorted BAM file
# ====================================

# Set up filter bam output files
log_msg info "Post-alignment filtering step..."
FILT_BAM_PREFIX=${LIBRARY}.filt.srt
FILT_BAM_FILE=$OUTPUT_DIR/${FILT_BAM_PREFIX}.bam
TMP_FILT_BAM_PREFIX=${FILT_BAM_PREFIX}.nmsrt
TMP_FILT_BAM_FILE=$TEMP_DIR/${TMP_FILT_BAM_PREFIX}.bam
MAPQ_THRESH=30

$SAMTOOLS_PATH/samtools view -F 1804 -f 2 -q ${MAPQ_THRESH} -u ${RAW_BAM_FILE} | $SAMTOOLS_PATH/samtools sort -n - $TEMP_DIR/${TMP_FILT_BAM_PREFIX} # Will produce name sorted BAM 

# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM
samtools fixmate -r ${TMP_FILT_BAM_FILE} $TEMP_DIR/${LIBRARY}.fixmate.tmp
$SAMTOOLS_PATH/samtools view -F 1804 -f 2 -u $TEMP_DIR/${LIBRARY}.fixmate.tmp | $SAMTOOLS_PATH/samtools sort - $OUTPUT_DIR/${FILT_BAM_PREFIX}
rm $TEMP_DIR/${LIBRARY}.fixmate.tmp
rm ${TMP_FILT_BAM_FILE}

# =============
# Mark duplicates
# =============
log_msg info "Marking and deduplication step..."
TMP_FILT_BAM_FILE=$OUTPUT_DIR/${FILT_BAM_PREFIX}.dupmark.bam
DUP_FILE_QC=$OUTPUT_DIR/${FILT_BAM_PREFIX}.dup.qc
PICARD_JAR=$(get_software_dir picard-tools-1.88/MarkDuplicates.jar)
CMD="java -Xmx4G -jar ${PICARD_JAR} INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
run_cmd "$CMD" "$TMP_FILT_BAM_FILE $DUP_FILE_QC"
mv $TMP_FILT_BAM_FILE $FILT_BAM_FILE

# ============================
# Remove duplicates
# Index final position sorted BAM
# Create final name sorted BAM
# ============================

# Set up final bam output files
FINAL_BAM_PREFIX=${LIBRARY}.filt.nodup.srt
FINAL_BAM_FILE=$OUTPUT_DIR/${FINAL_BAM_PREFIX}.bam #to be stored
FINAL_BAM_INDEX_FILE=$OUTPUT_DIR/${FINAL_BAM_PREFIX}.bai
FINAL_BAM_FILE_MAPSTATS=$OUTPUT_DIR/${FINAL_BAM_PREFIX}.flagstat.qc #QC file
FINAL_NMSRT_BAM_PREFIX=${LIBRARY}.filt.nodup.nmsrt
FINAL_NMSRT_BAM_FILE=$OUTPUT_DIR/${FINAL_NMSRT_BAM_PREFIX}.bam # To be stored

$SAMTOOLS_PATH/samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE} # Final bam file
$SAMTOOLS_PATH/samtools sort -n ${FINAL_BAM_FILE} $OUTPUT_DIR/${FINAL_NMSRT_BAM_PREFIX} 
$SAMTOOLS_PATH/samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE} # Index Final BAM file
$SAMTOOLS_PATH/samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}


# =============================
# Compute library complexity
# =============================
# Sort by name
# convert to bedPE and obtain fragment coordinates
# sort by position and strand
# Obtain unique count statistics

# Set up library complexity output file
PBC_FILE_QC=$OUTPUT_DIR/${FINAL_BAM_PREFIX}.pbc.qc
BEDTOOLS=$(get_software_dir bedtools2-2.20.0/bin/bedtools)

$SAMTOOLS_PATH/samtools sort -n ${FILT_BAM_FILE} $OUTPUT_DIR/${LIBRARY}.srt.tmp
$BEDTOOLS bamtobed -bedpe -i $OUTPUT_DIR/${LIBRARY}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}
rm $OUTPUT_DIR/${LIBRARY}.srt.tmp.bam
rm $FILT_BAM_FILE

# ===================
# Create tagAlign and BEDPE file
# ===================
# Set up BEDPE and tag align output files
log_msg info "Create BEDPE files..."
FINAL_BEDPE_FILE=$OUTPUT_DIR/${FINAL_NMSRT_BAM_PREFIX}.bedpe.gz
FINAL_TA_FILE=$OUTPUT_DIR/${FINAL_BAM_PREFIX}.PE2SE.tagAlign.gz
$BEDTOOLS bamtobed -bedpe -mate1 -i ${FINAL_NMSRT_BAM_FILE} | gzip -c > ${FINAL_BEDPE_FILE}

log_msg info "Create tag align files..."
zcat ${FINAL_BEDPE_FILE} | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -c > ${FINAL_TA_FILE}

# =================================
# Subsample tagAlign file
# Restrict to one read end per pair for CC analysis
# ================================
NREADS=15000000
SUBSAMPLED_TA_FILE="$OUTPUT_DIR/${LIBRARY}.filt.nodup.sample.$((NREADS / 1000000)).MATE1.tagAlign.gz"
zcat ${FINAL_BEDPE_FILE} | grep -v “chrM” | shuf -n ${NREADS} --random-source=${FINAL_TA_FILE}  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"N","1000",$9}' | gzip -c > ${SUBSAMPLED_TA_FILE}


log_msg info "Calcuate cross-correlation QC scores..."
# Set up CC QC output files
CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"
# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag
SPP_NODUPS=/gpfs/gpfs1/home/schhetri/spp_install/phantompeakqualtools/run_spp_nodups.R
Rscript $SPP_NODUPS -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE}
sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > $TEMP_DIR/temp
mv $TEMP_DIR/temp ${CC_SCORES_FILE}

### There seems to be bug on this line of code, since it outputs the same $CC_SCORES_FILE if re-run, and that would lead to repetition of the
### same line, and thus when you cut that line with $(cut -f 3 $CC_SCORES_FILE) then you are passing the 2 fragment length. This would kill
### the program which can only acccomodate only single fragment length argument. So, maybe you could just integrate $(cut -f 3 $CC_SCORES_FILE | tr "\n" "\t" | cut -f 1)
### in submit_peak_calls.sh script ; this would avoid the duplication of the same fragment. 
#SPP_NODUPS=/gpfs/gpfs1/home/snewberry/interrupt/encode/spp_and_idr/spp_install/phantompeakqualtools/run_spp_nodups.R
#SPP_NODUPS=/gpfs/gpfs1/home/schhetri/spp_install/phantompeakqualtools/run_spp_nodups.R
#CMD="Rscript $SPP_NODUPS -c=$SUBSAMPLED_TA_FILE -p=$NTHREADS -rf -filtchr=chrM -savp=$CC_PLOT_FILE -out=$CC_SCORES_FILE"
#run_cmd "$CMD" "$CC_SCORES_FILE $CC_PLOT_FILE"
#sed -ri 's/,[^\t]+//g' $CC_SCORES_FILE

# ========================
# Create pseudoReplicates
# =======================
# Set up PR output files
log_msg info "Creating self-psuedoreplicates..."
PR_PREFIX="${LIBRARY}.filt.nodup"
PR1_TA_FILE="$OUTPUT_DIR/${PR_PREFIX}.PE2SE.pr1.tagAlign.gz"
PR2_TA_FILE="$OUTPUT_DIR/${PR_PREFIX}.PE2SE.pr2.tagAlign.gz"
joined="$TEMP_DIR/temp.bedpe"

# Make temporary fake BEDPE file from FINAL_TA_FILE
zcat ${FINAL_TA_FILE} | sed 'N;s/\n/\t/' > $joined

# Get total number of read pairs
nlines=$( cat ${joined} | wc -l )
nlines=$(( (nlines + 1) / 2 ))

# Shuffle and split BEDPE file into 2 equal parts
cat ${joined} | shuf --random-source=${FINAL_TA_FILE}  | split -d -l ${nlines} - $TEMP_DIR/${PR_PREFIX} # Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01
# Convert fake BEDPE into standard tagAlign file
awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "$TEMP_DIR/${PR_PREFIX}00" | gzip -c > ${PR1_TA_FILE}
rm "$TEMP_DIR/${PR_PREFIX}00"
awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "$TEMP_DIR/${PR_PREFIX}01" | gzip -c > ${PR2_TA_FILE}
rm "$TEMP_DIR/${PR_PREFIX}01"
#rm -f ${joined}

log_msg info "Tag align script complete. Output is in $OUTPUT_DIR"

echo -e "\nstarting the script for the pooling of tag aligns:"
### Calling the script for pooling of the true tag align and control tag align:
#bsub -We 24:00 -n 1 -R span[hosts=1] -J "Call Pooling data" -o $LOG_DIR/pool_data.out $RUN_PATH/submit_PoolDataSets.main.sh
