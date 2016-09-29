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
declare -r MODULE_NAME=All
declare -r MODULE_VERSION=0.1.5.160912

#########
# Const #
#########
declare -ri DEFAULT_NUM_THREADS=8

#########
# Usage #
#########
function usage {
cat << EOF
${PACKAGE_NAME}::${MODULE_NAME}
======================================================
Run CCS + IsoSeq Classify + IsoAux + Cluster + Reports
======================================================

OPTIONS:
        -h      Show usage
        -v      Print version

${REQUIRED}[ required ]
        -c      Config CSV file
${OPTIONAL}[ optional ]
        -o      Output folder, will create if not exist. Default: ${PWD}
        -t      Number of CPU to use in each job. Default: ${DEFAULT_NUM_THREADS}
        -J      Name of this Job
        -E      Email address to be notified when jobs is done
${ADVANCED}[ advanced ]
        -D      use debug mode (bash -x). Default: off
EOF
echo -e "${FONT_COLOR_RESET}"
}

##########
# Config #
##########
declare -a REQUIRED_PROGRAMS=('smrtshell' 'readlink' 'find' 'gmap' 'gmap_build' 'samtools' 'picard' \
                            'perl' 'Rscript' 'faSize' 'trim_isoseq_polyA' 'faSize' 'bedToBigBedl' 'colmerge' \
                            'bedToGenePred' 'genePredToGtf' 'gffcompare' 'cuffmerge' 'mrna_size_from_gff' \
                            'bam2wig.py' 'computeMatrix' 'computeMatrix' \
                            )

#############################
# ARGS reading and checking #
#############################
while getopts "hvc:o:t:J:E:D" OPTION; do
    case $OPTION in
        h)  usage && exit 0 ;;
        v)  echo ${PACKAGE_NAME}::${MODULE_NAME} v${MODULE_VERSION} && exit 0;;
        c)  declare -x ConfigCsvFile=$(readlink -f ${OPTARG});;
        o)  declare OutputDir=${OPTARG};;
        t)  declare -i Threads=${OPTARG};; 
        J)  declare -x JobName="${OPTARG}";;
        E)  declare EmailAdd=${OPTARG};;
        D)  set -x;;
        *)  usage && exit 1;;
    esac
done
[[ -z ${ConfigCsvFile} ]] && echo2 "You have to provide a sample csv file with -c option" error
[[ -z $Threads ]] && declare -i Threads=${DEFAULT_NUM_THREADS}
[[ -z ${JobName} ]] && declare -x JobName="$(date)"

[[ -z $OutputDir ]] && OutputDir=${RANDOM}.out
declare OutDirFull=$(readlink -f ${OutputDir})
mkdir -p $OutputDir || echo2 "Don't have permission to create folder $OutputDir" error
cd $OutputDir || echo2 "Don't have permission to access folder $OutputDir" error
( touch a && rm -f a ) || echo2 "Don't have permission to write in folder $OutputDir" error

mkdir -p annotation log jobs jobout fasta bam gff bed fofn table pdf html bigWig

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
declare gmap_job_ids=""
declare picard_job_ids=""
declare cluster_jobids=""
declare -a flncfiles=()
declare -a flncsizefiles=()
declare -a coveragefiles=()
declare -a genomebamfiles=()
declare -a genomegtffiles=()
declare -a bigWigForwardFiles=()
declare -a bigWigReverseFiles=()
declare -a transcriptomerefs=()
declare -a clusterMappedGtfs=()
# declare -a primerinfofiles=()
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
    
    declare genegff=${ANNOTATION_DIR}/${genome}.genes.gtf
    [[ ! -f ${genegff} ]] && echo2 "Cannot find transcriptome file ${genegff}, please move it there or generate a symbol link" error
    transcriptomerefs+=(${genegff})
    
    declare refflatfile=${ANNOTATION_DIR}/${genome}.refFlat.txt
    [[ ! -f ${refflatfile} ]] && echo2 "Cannot find refFlat file ${refflatfile}, please move it there or generate a symbol link" error

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

    # build gmap index if not exist
    declare gmapindex=${GMAP_INDEX_DIR}/${genome}
    if ! gmapIndexCheck ${GMAP_INDEX_DIR} $genome; then
        if [[ ! -f jobs/${genome}.gmapbuild.sh ]]; then
            echo2 "cannot find gmap index ${gmapindex}, submiting a job to generate it" warning
            cat > jobs/${genome}.gmapbuild.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mkdir -p ${GMAP_INDEX_DIR}/${genome} \
&& gmap_build -d ${genome} -D ${GMAP_INDEX_DIR} ${genomefa}
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
            declare gmap_build_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${genome}.gmapbuild < jobs/${genome}.gmapbuild.sh | cut -f3 -d' ')
        else # already jobs/${genome}.gmapbuild.sh file
            echo2 "cannot find gmap index ${gmapindex}, but looks like a job has been submit to generate it" warning
        fi
    else
        echo2 "Use existing gmap index ${gmapindex}"
        declare gmap_build_jobid=""
    fi

    declare genome_preparation_jobid=${genome_size_jid}${gmap_build_jobid}
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
        declare gmap_ready=${genome_preparation_jobid},${rename_jobid}
    else
        declare gmap_ready=${genome_preparation_jobid},${ccsclassify_jobid}
    fi

    # 2. Run gmap alignemnt
    echo2 "Submit polyA trimming and gmap job"
    cat > jobs/${samplename}.gmap.sh << EOF
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
bash ${MYBIN}/gmap.sh \
    fasta/${samplename}.isoseq_flnc.trima.fa \
    ${GMAP_INDEX_DIR} \
    ${genome} \
    ${Threads}
fi

if [[ ! -f gff/${samplename}.isoseq_flnc.trima.${genome}.gtf ]]; then
bedtools bamtobed -bed12 -split -i bam/${samplename}.isoseq_flnc.trima.${genome}.sorted.bam \
    | bedToGenePred /dev/stdin /dev/stdout \
    | genePredToGtf file /dev/stdin gff/${samplename}.isoseq_flnc.trima.${genome}.gtf
fi

if [[ ! -f bigWig/${samplename}.isoseq_flnc.trima.${genome}.Forward.bw \
   || ! -f bigWig/${samplename}.isoseq_flnc.trima.${genome}.Reverse.bw ]]; then
    bam2wig.py \
        -s ${genomesize} \
        -i bam/${samplename}.isoseq_flnc.trima.${genome}.sorted.bam \
        -o bigWig/${samplename}.isoseq_flnc.trima.${genome} \
        -u -d '++,--' \
 && rm \
    bigWig/${samplename}.isoseq_flnc.trima.${genome}.Forward.wig \
    bigWig/${samplename}.isoseq_flnc.trima.${genome}.Reverse.wig 
fi
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
    declare -i gmap_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${samplename}.gmap -hold_jid ${gmap_ready} < jobs/${samplename}.gmap.sh | cut -f3 -d' ')
    gmap_job_ids=${gmap_jobid},${gmap_job_ids}

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
    declare -i piccard_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${samplename}.picard -hold_jid ${gmap_jobid} < jobs/${samplename}.CollectRnaSeqMetrics.sh | cut -f3 -d' ')
    picard_job_ids=${piccard_jobid},${picard_job_ids}
# clean up
    flncfiles+=("fasta/${samplename}.isoseq_flnc.trima.fa")
    flncsizefiles+=("table/${samplename}.isoseq_flnc.trima.sizes")
    coveragefiles+=(" table/${samplename}.isoseq_flnc.trima.coverage.tsv")
    genomebamfiles+=("bam/${samplename}.isoseq_flnc.trima.${genome}.sorted.bam")
    genomegtffiles+=("gff/${samplename}.isoseq_flnc.trima.${genome}.gtf")
    bigWigForwardFiles+=("bigWig/${samplename}.isoseq_flnc.trima.${genome}.Forward.bw")
    bigWigReverseFiles+=("bigWig/${samplename}.isoseq_flnc.trima.${genome}.Reverse.bw")

    # 4. run Cluster
    #  
    cat > jobs/${samplename}.Cluster.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if [[ ! -s jobout/${samplename}.clusterOut/all_quivered_hq.100_30_0.99.fasta ]]; then
    $SMRT_HOME/smrtcmds/bin/smrtshell -c \
        "pbtranscript.py cluster \
            fasta/${samplename}.isoseq_flnc.trima.fa \
            fasta/${samplename}.final.consensus.fa \
            --nfl_fa fasta/${samplename}.isoseq_nfl.fasta \
            -d jobout/${samplename}.clusterOut \
            --ccs_fofn fofn/${ccsname}.reads_of_insert.fofn \
            --bas_fofn fofn/${ccsname}.bax.fofn \
            --cDNA_size under1k \
            --quiver \
            --blasr_nproc ${Threads} \
            --quiver_nproc ${Threads} "
fi

if [[ ! -f fasta/${samplename}.all_quivered_hq.100_30_0.99.fasta \
     && -s jobout/${samplename}.clusterOut/all_quivered_hq.100_30_0.99.fasta ]]; then
    ln -s \
        $PWD/jobout/${samplename}.clusterOut/all_quivered_hq.100_30_0.99.fasta \
        fasta/${samplename}.all_quivered_hq.100_30_0.99.fasta
fi

if [[ ! -f bam/${samplename}.all_quivered_hq.100_30_0.99.${genome}.sorted.bam \
     && -f fasta/${samplename}.all_quivered_hq.100_30_0.99.fasta ]]; then
    bash ${MYBIN}/gmap.sh \
        fasta/${samplename}.all_quivered_hq.100_30_0.99.fasta \
        ${GMAP_INDEX_DIR} \
        ${genome} \
        ${Threads}
fi

# convert bam to gtf
if [[ ! -f gff/${samplename}.all_quivered_hq.100_30_0.99.${genome}.gtf \
     && -f bam/${samplename}.all_quivered_hq.100_30_0.99.${genome}.sorted.bam ]]; then
    bedtools bamtobed -bed12 -split -i bam/${samplename}.all_quivered_hq.100_30_0.99.${genome}.sorted.bam \
        | bedToGenePred /dev/stdin /dev/stdout \
        | genePredToGtf file /dev/stdin gff/${samplename}.all_quivered_hq.100_30_0.99.${genome}.gtf
fi
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
    declare -i cluster_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${samplename}.cluster -hold_jid ${ccsclassify_jobid} < jobs/${samplename}.Cluster.sh | cut -f3 -d' ')
    cluster_jobids=${cluster_jobid},${cluster_jobids}
    clusterMappedGtfs+=( gff/${samplename}.all_quivered_hq.100_30_0.99.${genome}.gtf )

done # end of for i in $(seq 0 $((SampleSize-1)))

echo2 "Submit job to draw length distribution on the flnc files"
# depends on
#   ${gmap_job_ids} which generate trimmed flnc files
cat > jobs/size.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if [[ ! -f table/draw_size_dis.Done ]]; then
    bash ${MYBIN}/draw_size_dis.sh ${flncsizefiles[@]} \
    && touch table/draw_size_dis.Done
fi
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
declare -i size_dis_jobid=$(${SUBMIT_CMD} -o log -e log -N job_size -hold_jid ${gmap_job_ids} < jobs/size.sh | cut -f3 -d' ')

echo2 "Submit job to run gffcompare"
# depends on 
#   $gmap_job_ids for bam/gff files
# generate $gffcompare_jobid
cat > jobs/gffcompare.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bash ${MYBIN}/count_gene_from_gff.sh jobs/${JOBNAME}.sh ${genegff} ${genomegtffiles[@]}
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
declare -i gffcompare_jobid=$(${SUBMIT_CMD} -o log -e log -N job_quantification -hold_jid ${gmap_job_ids} < jobs/gffcompare.sh | cut -f3 -d' ')

echo2 "Generate TSS and TES plot"
# depends on 
#   ${gmap_job_ids} which generate ${bigWigForwardFiles[@]} and ${bigWigReverseFiles[@]}
declare Annotation=${genegff} # TODO: currently only one annotation is supported 
declare AnnotationBed=annotation/$(basename ${Annotation%g[tf]f}bed)
declare AnnotationWatsonBed=${AnnotationBed%bed}watson.bed
declare AnnotationCrickBed=${AnnotationBed%bed}crick.bed
rm -f fofn/bam.fofn
for b in ${genomebamfiles[@]}; do 
    echo ${b} >> fofn/bam.fofn
done
cat > jobs/TSSTES.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# generate gennotation bed file from gtf
if [[ ! -f ${AnnotationBed} ]]; then
    gtfToGenePred ${Annotation} /dev/stdout \
  | genePredToBed /dev/stdin ${AnnotationBed}
fi
if [[ ! -f ${AnnotationBed%bed}watson.bed ]]; then
    awk '\$6=="+"' ${AnnotationBed} > ${AnnotationWatsonBed}
fi

if [[ ! -f ${AnnotationBed%bed}crick.bed ]]; then
    awk '\$6=="-"' ${AnnotationBed} > ${AnnotationCrickBed}
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

declare -i tss_tes_report=$(${SUBMIT_CMD} -o log -e log -N job_TSSTES -hold_jid ${gmap_job_ids} < jobs/TSSTES.sh | cut -f3 -d' ')


# coverage plot
cat > jobs/coverage.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bash ${MYBIN}/draw_coverage.sh ${coveragefiles[@]}
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
declare -i coverage_report=$(${SUBMIT_CMD} -o log -e log -N job_coverage -hold_jid ${picard_job_ids} < jobs/coverage.sh | cut -f3 -d' ')

# post clustering, compare annotation with gffcompare
echo -e ${clusterMappedGtfs[@]} | tr ' ' '\n' > fofn/cluster_mapped_gtf.fofn
cat > jobs/gffcompare.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mkdir -p jobout/clusterGffcompareOut

if [[ ! -f jobout/clusterGffcompareOut/gffcompare.Done ]]; then
    gffcompare -G -D \
        -r ${Annotation} \
        -o jobout/clusterGffcompareOut/gffcompare \
        ${clusterMappedGtfs[@]} \
    && touch jobout/clusterGffcompareOut/gffcompare.Done
fi

if [[ ! -f gff/gffcompare.combined.gtf \
   &&   -f jobout/clusterGffcompareOut/gffcompare.combined.gtf ]]; then
    ln -s $PWD/jobout/clusterGffcompareOut/gffcompare.combined.gtf gff/gffcompare.combined.gtf
fi

if [[ ! -f gffcompare.combined.chained.code.bb \
   &&   -f gff/gffcompare.combined.gtf ]]; then
    gtfToGenePred \
        gff/gffcompare.combined.gtf \
        /dev/stdout \
    | genePredToBed \
        /dev/stdin \
        /dev/stdout \
    | bedSort \
        /dev/stdin \
        bed/gffcompare.combined.bed.temp \
    && gawk 'BEGIN{FS=OFS="\t"}{ \
        if(ARGIND==1){ \
            h[\$1]=\$4; \
        } else { \
             if(h[\$4]=="=") \$9="0,0,0"; \
        else if(h[\$4]=="c") \$9="230,159,0"; \
        else if(h[\$4]=="j") \$9="255,0,0"; \
        else if(h[\$4]=="e") \$9="0,114,178"; \
        else if(h[\$4]=="i") \$9="0,158,115"; \
        else if(h[\$4]=="o") \$9="240,228,66"; \
        else if(h[\$4]=="x") \$9="213,94,0"; \
        else                 \$9="204,121,167"; \
        \$4=h[\$4]\$4; \
        print \$0}}' \
        jobout/clusterGffcompareOut/gffcompare.tracking \
        bed/gffcompare.combined.bed.temp \
        > bed/gffcompare.combined.bed \
    && bedToBigBed \
        bed/gffcompare.combined.bed \
        ${genomesize} \
        bed/gffcompare.combined.chained.code.bb
fi

# if [[ ! -f jobout/clusterCuffmergeOut/Done ]]; then
#     cuffmerge \
#         -p ${Threads} \
#         -g ${Annotation} \
#         -s ${genomefa} \
#         -o jobout/clusterCuffmergeOut \
#         fofn/cluster_mapped_gtf.fofn \
#     && ln -s $PWD/jobout/clusterCuffmergeOut/merged.gtf gff/Cuffmerged.gtf \
#     && touch jobout/clusterCuffmergeOut/Done
# fi
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
declare -i gffcompare_job=$(${SUBMIT_CMD} -o log -e log -N job_gffcompare -hold_jid ${cluster_jobids} < jobs/gffcompare.sh | cut -f3 -d' ')

# generate final report
echo2 "Generate final report"
# depends on 
#   ${size_dis_jobid} for length distribution
#   ${gffcompare_jobid} for quantification 
#   ${tss_tes_report} for png files
# generate ${final_html_report}
cat > jobs/html_report.sh << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bash ${MYBIN}/generate_Rmd.sh
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
declare -i final_html_report=$(${SUBMIT_CMD} -o log -e log -N job_generate_html -hold_jid ${gffcompare_jobid},${size_dis_jobid},${tss_tes_report},${coverage_report} < jobs/html_report.sh | cut -f3 -d' ')

# send notification
declare last_job=${gffcompare_job}
if [[ ! -z ${EmailAdd} ]]; then
    echo "bash ${MYBIN}/send_mail.sh \"${JobName}\" ${OutDirFull} ${EmailAdd}" \
    | declare -i notify_jobid=$(${SUBMIT_CMD} -o log -e log -N notify_done -hold_jid ${last_job} | cut -f3 -d' ')
fi