#!/usr/bin/env bash
#################################################################################
# Copyright (c) 2016-, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################
# author: Bo Han (bhan@pacb.com)

#########
# Basic #
#########
declare -r MODULE_NAME=Prokaryotic
declare -r MODULE_VERSION=0.0.2.161006

#########
# Const #
#########
declare -ri DEFAULT_NUM_THREADS=8

########################
# Variable declaration #
########################
declare -x  ConfigCsvFile
declare -x  OutputDir
declare -xi Threads=8
declare -x  JobName
declare -x  EmailAdd
declare     FTPAddress
declare     FTPUsername
declare     FTPPassword

#########
# Usage #
#########
function usage {
cat << EOF
${PACKAGE_NAME}::${MODULE_NAME}
======================================================
Run CCS + IsoSeq Classify + IsoAux + Cluster + Reports
======================================================

For Prokaryotic species (aka, no intron). 
The pipeline expects two files in the ${ANNOTATION_DIR} folder for a genome called "ecoli" for example. 
1. ecoli.fa: the genome fasta file
2. ecoli.genes.bed: the bed file with gene coordinates
If you do not have permission to write in ${ANNOTATION_DIR}, prepare those files in a different directory
and point the pipeline to that directory using option -A

OPTIONS:
        -h      Show usage
        -v      Print version

${REQUIRED}[ required ]
        -c      Config CSV file
${OPTIONAL}[ optional ]
        -o      Output folder, will create if not exist. Default: ${PWD}
        -A      Annotation directory, overwriting the default directory specified in config.sh file
        -t      Number of CPU to use in each job. Default: ${DEFAULT_NUM_THREADS}
        -J      Name of this Job. Default: datetime
        -E      Email address to be notified when jobs is done. Default: no notification
FTP settings
If all three of them are provided, the pipeline will upload bam/bigWig/bigBed files to this FTP,
and generate UCSC genome browser tracks.
        -F      FTP address
        -U      FTP username
        -P      FTP password
${ADVANCED}[ advanced ]
        -D      use debug mode (bash -x). Default: off
EOF
echo -e "${FONT_COLOR_RESET}"
}

##########
# Config #
##########
declare -a REQUIRED_PROGRAMS=('smrtshell' 'readlink' 'find' 'perl' 'Rscript' \
                            'bwa' 'samtools' 'picard' \
                            'trim_isoseq_polyA' 'mrna_size_from_gff' 'colmerge' \
                            'bedGraphToBigWig' 'bedToBigBed' 'bedToGenePred' 'faSize' 'genePredToGtf' \
                            'computeMatrix' 'computeMatrix' \
                            )

#############################
# ARGS reading and checking #
#############################
while getopts "hvc:o:t:J:E:DA:F:U:P:" OPTION; do
    case $OPTION in
        h)  usage && exit 0 ;;
        v)  echo ${PACKAGE_NAME}::${MODULE_NAME} v${MODULE_VERSION} && exit 0;;
        c)  ConfigCsvFile=$(readlink -f ${OPTARG});;
        A)  ANNOTATION_DIR=$(readlink -f ${OPTARG});;
        o)  OutputDir=${OPTARG};;
        t)  Threads=${OPTARG};; 
        J)  JobName="${OPTARG}";;
        E)  EmailAdd=${OPTARG};;
        F)  FTPAddress=${OPTARG};;
        U)  FTPUsername=${OPTARG};;
        P)  FTPPassword=${OPTARG};;
        D)  set -x;;
        *)  usage && exit 1;;
    esac
done
[[ -z ${ConfigCsvFile} ]] && echo2 "You have to provide a sample csv file with -c option" error
[[ -z ${JobName} ]] && declare -x JobName="$(date)"

[[ -z $OutputDir ]] && OutputDir=${RANDOM}.out
declare OutDirFull=$(readlink -f ${OutputDir})
mkdir -p $OutputDir || echo2 "Don't have permission to create folder $OutputDir" error
cd $OutputDir || echo2 "Don't have permission to access folder $OutputDir" error
( touch a && rm -f a ) || echo2 "Don't have permission to write in folder $OutputDir" error

mkdir -p annotation log jobs jobout fasta bam bed fofn table pdf html bigWig track

for program in "${REQUIRED_PROGRAMS[@]}"; do binCheck $program; done

#########
# Begin #
#########
echo2 "Parse sample csv file"
declare JOBNAME=$(md5sum ${ConfigCsvFile} | cut -d' ' -f1)
python ${MYBIN}/parse_csv.py ${ConfigCsvFile} > jobs/${JOBNAME}.sh \
    || echo2 "Error parsing sample csv file" error
source jobs/${JOBNAME}.sh

echo2 "Submit jobs for CCS and Classify"
declare bwa_job_ids=""
declare picard_job_ids=""
declare -a flncfiles=()
declare -a flncsizefiles=()
declare -a coveragefiles=()
declare -a genomebamfiles=()
declare -a bigWigForwardFiles=()
declare -a bigWigReverseFiles=()
declare -a transcriptomerefs=()
for i in $(seq 0 $((SampleSize-1))); do
    declare samplename=${SampleNames[$i]}
    echo2 "Process ${samplename}"
    declare cellpath=${SampleCellPaths[$i]}
    declare cellid=$(echo ${cellpath} | md5sum | cut -d' ' -f1)
    declare barcodefile=${BarcodeFiles[$i]}
    declare barcode_id=${BarcodeNumbers[$i]}
    declare genome=${Genomes[$i]}
    
    declare genomefa=${ANNOTATION_DIR}/${genome}.fa
    [[ ! -f ${genomefa} ]] && echo2 "Cannot find genome fasta file ${genomefa}, please move it there or generate a symbol link" error
    
    declare genebed=${ANNOTATION_DIR}/${genome}.genes.bed
    [[ ! -f ${genebed} ]] && echo2 "Cannot find transcriptome file ${genebed}, please move it there or generate a symbol link" error
    transcriptomerefs+=(${genebed})
    
    declare refflatfile=${ANNOTATION_DIR}/${genome}.refFlat.txt
    if [[ ! -f ${refflatfile} ]]; then
        refflatfile=annotation/${genome}.refFlat.txt
        if [[ -f ${refflatfile} ]]; then
            echo2 "Cannot find refFlat file ${refflatfile}, generating one from the bed" warning
            awk 'BEGIN{FS=OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\t%s,\t%s,\n", $4,$4,$1,$6,$2,$3,$2,$3,$2,$3}' ${genebed} > ${refflatfile}
        fi
    fi

    if [[ -f ${ANNOTATION_DIR}/${genome}.sizes ]]; then 
        declare genomesize=${ANNOTATION_DIR}/${genome}.sizes;
        declare genome_size_jid=""
    elif [[ -f ${ANNOTATION_DIR}/${genome}.fa.sizes ]]; then
        declare genomesize=${ANNOTATION_DIR}/${genome}.fa.sizes;
        declare genome_size_jid=""
    else
        echo2 "Missing genome size file, generate one"
        cat > jobs/${genome}.size.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
faSize -detailed -tab ${genomefa} > annotation/${genome}.sizes.tsv; 
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
        declare genome_size_jid=$(${SUBMIT_CMD} -o log -e log -N job_${genome}.size < jobs/${genome}.size.sh | cut -f3 -d' ')","
        declare genomesize=annotation/${genome}.sizes.tsv
    fi

    # build bwa index
    declare bwaindex=${BWA_INDEX_DIR}/${genome}
    if ! bwaIndexCheck ${BWA_INDEX_DIR}/$genome; then
        echo2 "Cannot find bwa index ${bwaindex}" warning
        if ! dirWritable ${BWA_INDEX_DIR}/; then
            BWA_INDEX_DIR=annotation
            bwaindex=annotation/$genome
        fi 
        if [[ ! -f jobs/${genome}.bwabuild.sh ]]; then
            echo2 "Submiting a job to generate one at $bwaindex" warning
            cat > jobs/${genome}.bwabuild.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bwa index -p ${bwaindex} ${genomefa}
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
            declare bwa_build_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${genome}.bwabuild < jobs/${genome}.bwabuild.sh | cut -f3 -d' ')
        else # already jobs/${genome}.bwabuild.sh file
            echo2 "Looks like a job has been submit to generate it" warning
        fi
    else # bwaIndexCheck ${BWA_INDEX_DIR}/$genome
        echo2 "Use existing bwa index ${bwaindex}"
        declare bwa_build_jobid=""
    fi

    declare genome_preparation_jobid=${genome_size_jid}${bwa_build_jobid}
    # barcoded?
    if [[ ${barcode_id} == 'NA' ]]; then
        echo2 "${samplename} has no barcode"
        declare ccsname=${samplename}
        declare pbtrascript_option=""
        declare -i total_barcode_number=1
        declare split_barcode_cmd="echo"
    else # barcoded 
        echo2 "${samplename} is barcoded"
        declare ccsname=$(perl -pe 's/[^A-Za-z0-9]/_/g' <<< $(basename $(dirname ${cellpath})))$(perl -pe 's/[^A-Za-z0-9]/_/g' <<< $(basename ${cellpath}))
        declare pbtrascript_option="--primer=${barcodefile}"
        declare -i total_barcode_number=$(grep '>' ${barcodefile} | wc -l)
        let total_barcode_number/=2 # barcode files are in pair
        declare split_barcode_cmd="python ${MYBIN}/split_flnc_barcodes.py fasta/${ccsname}.isoseq_flnc.fasta ${total_barcode_number}; python ${MYBIN}/split_flnc_barcodes.py fasta/${ccsname}.isoseq_nfl.fasta ${total_barcode_number}"
    fi

    # 1. CCS
    # 1.1 find bax.h5 file and generate fofn files for CCS
    echo2 "Generate fofn for ${ccsname}"
    if [[ ! -f log/${cellid} ]]; then
        # this cell has not been submited to Process
        find $(readlink -f "${cellpath}") -name "*bax.h5" | xargs -I {} readlink -f {} > fofn/${ccsname}.bax.fofn
        echo2 "Submit CCS and Classify job for ${ccsname}"
        # 1.2 genearte scripts to run CCS
        cat > jobs/${ccsname}.ccs.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# prepare
mkdir -p pdf/${ccsname}

# run CCS
if [[ ! -f fasta/${ccsname}.CCS.fa ]]; then
$SMRT_HOME/smrtcmds/bin/smrtshell -c "ConsensusTools.sh CircularConsensus \
    --minFullPasses 0 \
    --minPredictedAccuracy 75 \
    --numThreads ${Threads} \
    --fofn fofn/${ccsname}.bax.fofn \
    -o jobout/${ccsname}.CCS " \
&& cat jobout/${ccsname}.CCS/*ccs.fasta > fasta/${ccsname}.CCS.fa
fi

# run CCS report
if [[ ! -f pdf/${ccsname}.roi_readlength_hist.png \
   || ! -f pdf/${ccsname}.roi_npasses_hist.png \
   || ! -f pdf/${ccsname}.roi_accuracy_hist.png ]]; then
    rm -f fofn/${ccsname}.reads_of_insert.fofn
    for f in \$(find jobout/${ccsname}.CCS -name "*ccs.h5"); do 
        readlink -f \${f} >> fofn/${ccsname}.reads_of_insert.fofn
    done
    $SMRT_HOME/smrtcmds/bin/smrtshell -c "reads_of_insert_report.py \
        --debug \
        --output-dir pdf/${ccsname} \
        fofn/${ccsname}.reads_of_insert.fofn \
        log/${ccsname}.reads_of_insert_report.json" \
 && mv pdf/${ccsname}/roi_accuracy_hist.png    pdf/${ccsname}.roi_accuracy_hist.png  \
 && mv pdf/${ccsname}/roi_npasses_hist.png     pdf/${ccsname}.roi_npasses_hist.png \
 && mv pdf/${ccsname}/roi_readlength_hist.png  pdf/${ccsname}.roi_readlength_hist.png
fi

# run classify
if [[ ! -f fasta/${ccsname}.isoseq_flnc.fasta ]]; then
$SMRT_HOME/smrtcmds/bin/smrtshell -c "pbtranscript.py classify \
    --cpus ${Threads} \
    --primer_search_window 200 \
    --min_dist_from_end 200 \
    --min_seq_len 300 \
    ${pbtrascript_option} \
    --flnc fasta/${ccsname}.isoseq_flnc.fasta \
    --nfl fasta/${ccsname}.isoseq_nfl.fasta \
    -d jobout/${ccsname}.classifyOut \
    fasta/${ccsname}.CCS.fa \
    fasta/${ccsname}.isoseq_draft.fasta" 
fi

# split barcode
if [[ ! -f log/${ccsname}.barcode_splited.Done ]]; then
    ${split_barcode_cmd} && touch log/${ccsname}.barcode_splited.Done
fi

# run classify report
if [[ ! -f pdf/${ccsname}.fulllength_nonchimeric_readlength_hist.png ]]; then
    $SMRT_HOME/smrtcmds/bin/smrtshell -c "isoseq_classify_report.py \
        fasta/${ccsname}.isoseq_flnc.fasta \
        fasta/${ccsname}.isoseq_draft.classify_summary.txt \
        table/${ccsname}.isoseq_classify.json \
        -o pdf/${ccsname}" \
    && mv pdf/${ccsname}/fulllength_nonchimeric_readlength_hist.png pdf/${ccsname}.fulllength_nonchimeric_readlength_hist.png
fi
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
        # 1.4 submit CCS script
        declare -i ccsclassify_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${ccsname}.CCS-classify < jobs/${ccsname}.ccs.sh | cut -f3 -d' ')
        echo ${ccsclassify_jobid} > log/${cellid}
    else 
        # the CCS and classify of this cell has been submiited
        echo2 "The same cell has been submitted, most likely because it is barcoded" warning
        declare -i ccsclassify_jobid=$(cat log/${cellid})
    fi
    # rename the file accoring to the barcode
    if [[ ${barcode_id} != 'NA' ]]; then
        cat > jobs/${samplename}.rename.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if [[ ! -f fasta/${samplename}.isoseq_flnc.fasta && -f fasta/${ccsname}.isoseq_flnc.fasta.barcode${barcode_id} ]]; then
    mv fasta/${ccsname}.isoseq_flnc.fasta.barcode${barcode_id} fasta/${samplename}.isoseq_flnc.fasta
fi

if [[ ! -f fasta/${samplename}.isoseq_nfl.fasta && -f fasta/${ccsname}.isoseq_nfl.fasta.barcode${barcode_id} ]]; then
    mv fasta/${ccsname}.isoseq_nfl.fasta.barcode${barcode_id}  fasta/${samplename}.isoseq_nfl.fasta
fi
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF

        declare -i rename_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${samplename}.rename -hold_jid ${ccsclassify_jobid} < jobs/${samplename}.rename.sh | cut -f3 -d' ')
        declare bwa_ready=${genome_preparation_jobid},${rename_jobid}
    else
        declare bwa_ready=${genome_preparation_jobid},${ccsclassify_jobid}
    fi

    # 2. Run bwa alignemnt
    echo2 "Submit polyA trimming and bwa job"
    cat > jobs/${samplename}.bwa.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if [[ ! -f fasta/${samplename}.isoseq_flnc.trima.fa ]]; then
trim_isoseq_polyA -t ${Threads} \
    -i fasta/${samplename}.isoseq_flnc.fasta \
    > fasta/${samplename}.isoseq_flnc.trima.fa \
    2> log/${samplename}.isoseq_flnc.trima.log 
fi

if [[ ! -f table/${samplename}.isoseq_flnc.trima.sizes ]]; then
faSize -detailed -tab \
    fasta/${samplename}.isoseq_flnc.trima.fa \
    > table/${samplename}.isoseq_flnc.trima.sizes
fi

if [[ ! -f bam/${samplename}.isoseq_flnc.trima.${genome}.sorted.bam ]]; then
bash ${MYBIN}/bwa.sh \
    fasta/${samplename}.isoseq_flnc.trima.fa \
    ${bwaindex} \
    ${Threads}
fi

if [[ ! -f bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bed \
     && -f bam/${samplename}.isoseq_flnc.trima.${genome}.sorted.bam ]]; then
    bedtools bamtobed \
        -bed12 \
        -split \
        -i bam/${samplename}.isoseq_flnc.trima.${genome}.sorted.bam \
        > bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bed
fi

# if [[ ! -f gff/${samplename}.isoseq_flnc.trima.${genome}.gtf \
#     &&  -f bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bed ]]; then
#     bedToGenePred bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bed /dev/stdout \
#     | genePredToGtf file /dev/stdin gff/${samplename}.isoseq_flnc.trima.${genome}.gtf
# fi

if [[ ! -f bigWig/${samplename}.isoseq_flnc.trima.${genome}.Forward.bw \
   || ! -f bigWig/${samplename}.isoseq_flnc.trima.${genome}.Reverse.bw ]]; then
    bedtools genomecov -split -bg \
        -strand + \
        -g ${genomesize} \
        -i bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bed \
    > bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bedGraph \
    && bedSort \
        bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bedGraph \
        bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bedGraph \
    && bedGraphToBigWig \
        bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bedGraph \
        ${genomesize} \
        bigWig/${samplename}.isoseq_flnc.trima.${genome}.Forward.bw
    # negate the values for Crick strand
    bedtools genomecov -split -bg \
        -strand - \
        -g ${genomesize} \
        -i bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bed \
    | awk 'BEGIN{FS=OFS="\t"}{\$4=-\$4; print \$0}' \
        > bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bedGraph \
    && bedSort \
        bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bedGraph \
        bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bedGraph \
    && bedGraphToBigWig \
        bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bedGraph \
        ${genomesize} \
        bigWig/${samplename}.isoseq_flnc.trima.${genome}.Reverse.bw \
    && rm bed/${samplename}.isoseq_flnc.trima.${genome}.sorted.bedGraph
fi
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
    declare -i bwa_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${samplename}.bwa -hold_jid ${bwa_ready} < jobs/${samplename}.bwa.sh | cut -f3 -d' ')
    bwa_job_ids=${bwa_jobid},${bwa_job_ids}

# draw coverage plot 
    cat > jobs/${samplename}.CollectRnaSeqMetrics.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if [[ ! -f table/${samplename}.isoseq_flnc.trima.coverage ]]; then
    picard CollectRnaSeqMetrics \
        REF_FLAT=${refflatfile} \
        INPUT=bam/${samplename}.isoseq_flnc.trima.${genome}.sorted.bam \
        O=table/${samplename}.isoseq_flnc.trima.coverage \
        STRAND=FIRST_READ_TRANSCRIPTION_STRAND
fi

awk 'BEGIN{FS=OFS="\t"}{if(\$1>0 && \$1<101) print \$1,\$2}' table/${samplename}.isoseq_flnc.trima.coverage \
    > table/${samplename}.isoseq_flnc.trima.coverage.tsv
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
    declare -i piccard_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${samplename}.picard -hold_jid ${bwa_jobid} < jobs/${samplename}.CollectRnaSeqMetrics.sh | cut -f3 -d' ')
    picard_job_ids=${piccard_jobid},${picard_job_ids}
# clean up
    flncfiles+=("fasta/${samplename}.isoseq_flnc.trima.fa")
    flncsizefiles+=("table/${samplename}.isoseq_flnc.trima.sizes")
    coveragefiles+=(" table/${samplename}.isoseq_flnc.trima.coverage.tsv")
    genomebamfiles+=("bam/${samplename}.isoseq_flnc.trima.${genome}.sorted.bam")
    bigWigForwardFiles+=("bigWig/${samplename}.isoseq_flnc.trima.${genome}.Forward.bw")
    bigWigReverseFiles+=("bigWig/${samplename}.isoseq_flnc.trima.${genome}.Reverse.bw")

done # end of for i in $(seq 0 $((SampleSize-1)))

echo2 "Submit job to draw length distribution on the flnc files"
# depends on
#   ${bwa_job_ids} which generate trimmed flnc files
cat > jobs/size.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if [[ ! -f table/draw_size_dis.Done ]]; then
    bash ${MYBIN}/draw_size_dis.sh ${flncsizefiles[@]} \
    && touch table/draw_size_dis.Done
fi
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
declare -i size_dis_jobid=$(${SUBMIT_CMD} -o log -e log -N job_size -hold_jid ${bwa_job_ids} < jobs/size.sh | cut -f3 -d' ')

echo2 "Generate TSS and TES plot"
# depends on 
#   ${bwa_job_ids} which generate ${bigWigForwardFiles[@]} and ${bigWigReverseFiles[@]}
declare Annotation=${genebed} # TODO: currently only one annotation is supported 
declare AnnotationBed=annotation/$(basename ${genebed})
declare AnnotationWatsonBed=${AnnotationBed%bed}watson.bed
declare AnnotationCrickBed=${AnnotationBed%bed}crick.bed
rm -f fofn/bam.fofn
for b in ${genomebamfiles[@]}; do 
    echo ${b} >> fofn/bam.fofn
done
cat > jobs/TSSTES.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# generate gennotation bed file from gtf

if [[ ! -f ${AnnotationBed%bed}watson.bed ]]; then
    awk '\$6=="+"' ${Annotation} > ${AnnotationWatsonBed}
fi

if [[ ! -f ${AnnotationBed%bed}crick.bed ]]; then
    awk '\$6=="-"' ${Annotation} > ${AnnotationCrickBed}
fi

if [[ ! -f pdf/TSS.crick.png \
   || ! -f pdf/TSS.watson.png \
   || ! -f pdf/TES.crick.png \
   || ! -f pdf/TES.watson.png ]]; then
    for refp in TSS TES; do
        computeMatrix    reference-point -q --referencePoint \${refp} --sortRegions descend -S ${bigWigForwardFiles[@]} -R ${AnnotationWatsonBed} -a 500 -b 500 --maxThreshold 20 --missingDataAsZero -out table/\${refp}.+.matrix.gz -p ${Threads} \
        && computeMatrix reference-point -q --referencePoint \${refp} --sortRegions ascend  -S ${bigWigReverseFiles[@]} -R ${AnnotationCrickBed}  -a 500 -b 500 --maxThreshold 20 --missingDataAsZero -out table/\${refp}.-.matrix.gz -p ${Threads} \
        && plotHeatmap -m table/\${refp}.+.matrix.gz -out pdf/\${refp}.watson.png --samplesLabel ${SampleNames[@]} --regionsLabel $(basename ${Annotation}) --refPointLabel \${refp} --plotTitle "distance to Watson \${refp}" --colorMap RdBu \
        && plotHeatmap -m table/\${refp}.-.matrix.gz -out pdf/\${refp}.crick.png  --samplesLabel ${SampleNames[@]} --regionsLabel $(basename ${Annotation}) --refPointLabel \${refp} --plotTitle "distance to Crick \${refp}"  --colorMap RdBu \
        || echo "failed to generate TSS and TES plots"
    done
fi

# this is very time consuming
# geneBody_coverage.py -r ${AnnotationBed} -i fofn/bam.fofn -f pdf -o pdf/RNA_coverage || echo "failed to run geneBody_coverage.py"
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF

declare -i tss_tes_report=$(${SUBMIT_CMD} -o log -e log -N job_TSSTES -hold_jid ${bwa_job_ids} < jobs/TSSTES.sh | cut -f3 -d' ')


# coverage plot
cat > jobs/coverage.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bash ${MYBIN}/draw_coverage.sh ${coveragefiles[@]}
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
declare -i coverage_report=$(${SUBMIT_CMD} -o log -e log -N job_coverage -hold_jid ${picard_job_ids} < jobs/coverage.sh | cut -f3 -d' ')

# generate final report
echo2 "Generate final report"
# depends on 
#   ${size_dis_jobid} for length distribution
#   ${gffcompare_jobid} for quantification 
#   ${tss_tes_report} for png files
# generate ${final_html_report}
cat > jobs/html_report.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bash ${MYBIN}/generate_Rmd_for_prok.sh
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
declare -i final_html_report=$(${SUBMIT_CMD} -o log -e log -N job_generate_html -hold_jid ${size_dis_jobid},${tss_tes_report},${coverage_report} < jobs/html_report.sh | cut -f3 -d' ')

# getting ready to finish
declare last_job=${final_html_report}
# upload to ftp
if [[ ! -z ${FTPAddress} \
   && ! -z ${FTPUsername} \
   && ! -z ${FTPPassword} ]] ; then
   cat > jobs/ftp_upload_and_track.sh 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
( cd bam;    bash ${MYBIN}/upload_to_ftp_and_generate_track.sh $FTPAddress $FTPUsername $FTPPassword $JOBNAME *bam 1> ../track/UCSC.genome_browser.tracks)
( cd bed;    bash ${MYBIN}/upload_to_ftp_and_generate_track.sh $FTPAddress $FTPUsername $FTPPassword $JOBNAME *bb  1> ../track/UCSC.genome_browser.tracks)
( cd bigWig; bash ${MYBIN}/upload_to_ftp_and_generate_track.sh $FTPAddress $FTPUsername $FTPPassword $JOBNAME *bw  1> ../track/UCSC.genome_browser.tracks)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
    declare -i ftpjob=$(${SUBMIT_CMD} -o log -e log -N job_gffcompare -hold_jid ${last_job} < jobs/gffcompare.sh | cut -f3 -d' ')
    last_job=${ftpjob}
fi     

# send notification
if [[ ! -z ${EmailAdd} ]]; then
    echo "bash ${MYBIN}/send_mail.sh \"${JobName}\" ${OutDirFull} ${EmailAdd}" \
    | declare -i notify_jobid=$(${SUBMIT_CMD} -o log -e log -N notify_done -hold_jid ${last_job} | cut -f3 -d' ')
fi