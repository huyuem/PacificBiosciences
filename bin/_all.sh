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
declare -r MODULE_VERSION=0.0.1

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
============================================
Run CCS + IsoSeq Classify + IsoAux + Reports
============================================

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
declare -a REQUIRED_PROGRAMS=('smrtshell' 'readlink' 'find' 'gmap' 'gmap_build' 'samtools' \
                            'Rscript' 'faSize' 'trim_isoseq_polyA' 'faSize')

#############################
# ARGS reading and checking #
#############################
while getopts "hvc:o:t:J:E:D" OPTION; do
    case $OPTION in
        h)  usage && exit 0 ;;
        v)  echo ${PACKAGE_NAME}::${MODULE_NAME} v${MODULE_VERSION} && exit 0;;
        c)  declare ConfigCsvFile=$(readlink -f ${OPTARG});;
        o)  declare OutputDir=${OPTARG};;
        t)  declare -i Threads=${OPTARG};; 
        J)  declare JobName="${OPTARG}";;
        E)  declare EmailAdd=${OPTARG};;
        D)  set -x;;
        *)  usage && exit 1;;
    esac
done
[[ -z ${ConfigCsvFile} ]] && echo2 "You have to provide a sample csv file with -c option" error
[[ -z $Threads ]] && declare -i Threads=${DEFAULT_NUM_THREADS}
[[ -z ${JobName} ]] && declare JobName="$(date)"

[[ -z $OutputDir ]] && OutputDir=${RANDOM}.out
declare OutDirFull=$(readlink -f ${OutputDir})
mkdir -p $OutputDir || echo2 "Don't have permission to create folder $OutputDir" error
cd $OutputDir || echo2 "Don't have permission to access folder $OutputDir" error
( touch a && rm -f a ) || echo2 "Don't have permission to write in folder $OutputDir" error

mkdir -p log jobs fofn table pdf html

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
declare -a flncfiles=()
declare -a flncsizefiles=()
declare -a genomebamfiles=()
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
    
    # build gmap index if not exist
    declare gmapindex=${INDEX_DIR}/gmap_index/${genome}
    if ! gmapIndexCheck ${INDEX_DIR}/gmap_index $genome; then
        if [[ ! -f jobs/${genome}.gmapbuild.sh ]]; then
            echo2 "cannot find gmap index ${gmapindex}, submiting a job to generate it" warning
            echo "gmap_build -d ${genome} -D ${INDEX_DIR}/gmap_index/ ${genomefa}" > jobs/${genome}.gmapbuild.sh
            declare gmap_build_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${genome}.gmapbuild < jobs/${genome}.gmapbuild.sh | cut -f3 -d' ')
        else # already jobs/${genome}.gmapbuild.sh file
            echo2 "cannot find gmap index ${gmapindex}, but looks like a job has been submit to generate it" warning
        fi
    else
        echo2 "Use existing gmap index ${gmapindex}"
        declare gmap_build_jobid=""
    fi

    # barcoded?
    if [[ ${barcode_id} == 'NA' ]]; then
        echo2 "${samplename} has no barcode"
        declare ccsname=${samplename}
        declare pbtrascript_option=""
        declare -i total_barcode_number=1
        declare split_barcode_cmd="echo"
    else # barcoded 
        echo2 "${samplename} is barcoded"
        declare ccsname=${cellid}
        declare pbtrascript_option="--primer=${barcodefile}"
        declare -i total_barcode_number=$(grep '>' ${barcodefile} | wc -l)
        let total_barcode_number/=2 # barcode files are in pair
        declare split_barcode_cmd="python ${MYBIN}/split_flnc_barcodes.py ${ccsname}.isoseq_flnc.fasta ${total_barcode_number}"
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

[[ ! -f ${ccsname}.isoseq_flnc.fasta ]] \
&& $SMRT_HOME/smrtcmds/bin/smrtshell -c "ConsensusTools.sh CircularConsensus \
    --minFullPasses 0 \
    --minPredictedAccuracy 75 \
    --numThreads ${Threads} \
    --fofn fofn/${ccsname}.bax.fofn \
    -o ${ccsname}.CCS " \
&& cat ${ccsname}.CCS/*ccs.fasta > ${ccsname}.CCS.fa \
&& ls ${ccsname}.CCS/*.ccs.h5 > ${ccsname}.reads_of_insert.fofn \
&& $SMRT_HOME/smrtcmds/bin/smrtshell -c "pbtranscript.py classify \
    --cpus ${Threads} \
    --primer_search_window 200 \
    --min_dist_from_end 200 \
    --min_seq_len 300 \
    ${pbtrascript_option} \
    --flnc ${ccsname}.isoseq_flnc.fasta \
    --nfl ${ccsname}.isoseq_nfl.fasta \
    -d ${ccsname}.classifyOut \
    ${ccsname}.CCS.fa \
    ${ccsname}.isoseq_draft.fasta" \
&& ${split_barcode_cmd}

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
        echo "mv ${ccsname}.isoseq_flnc.fasta.barcode${barcode_id} ${samplename}.isoseq_flnc.fasta" > jobs/${samplename}.rename.sh
        declare -i rename_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${samplename}.rename -hold_jid ${ccsclassify_jobid} < jobs/${samplename}.rename.sh | cut -f3 -d' ')
        declare gmap_ready=${gmap_build_jobid},${rename_jobid}
    else
        declare gmap_ready=${gmap_build_jobid},${ccsclassify_jobid}
    fi

    # 2. Run gmap alignemnt
    echo2 "Submit polyA trimming and gmap job"
    cat > jobs/${samplename}.gmap.sh << EOF
# TODO: currently trim_isoseq_polyA causes segmentfault when ran on a working node through qsub for unknown reasons
trim_isoseq_polyA -t ${Threads} ${samplename}.isoseq_flnc.fasta > ${samplename}.isoseq_flnc.trima.fa 2> log/${samplename}.isoseq_flnc.trima.log \
&& faSize -detailed -tab ${samplename}.isoseq_flnc.trima.fa > ${samplename}.isoseq_flnc.trima.sizes \
&& bash ${MYBIN}/gmap.sh ${samplename}.isoseq_flnc.trima.fa ${INDEX_DIR}/gmap_index/ ${genome} ${Threads}

EOF
    declare -i gmap_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${samplename}.gmap -hold_jid ${gmap_ready} < jobs/${samplename}.gmap.sh | cut -f3 -d' ')
    gmap_job_ids=${gmap_jobid},${gmap_job_ids}

    flncfiles+=("${samplename}.isoseq_flnc.trima.fa")
    flncsizefiles+=("${samplename}.isoseq_flnc.trima.sizes")
    genomebamfiles+=("${samplename}.isoseq_flnc.trima.${genome}.sorted.bam")
done

echo2 "Submit job to draw length distribution on the flnc files"
cat > jobs/size.sh << EOF
bash ${MYBIN}/draw_size_dis.sh ${flncsizefiles[@]}
EOF
declare -i size_dis_jobid=$(${SUBMIT_CMD} -o log -e log -N job_size -hold_jid ${gmap_job_ids} < jobs/size.sh | cut -f3 -d' ')


# send notification
declare last_job=${size_dis_jobid}
[[ ! -z ${EmailAdd} ]] \
 && echo "bash ${MYBIN}/send_mail.sh \"${JobName}\" ${OutDirFull} ${EmailAdd}" \
 | declare -i notify_jobid=$(${SUBMIT_CMD} -o log -e log -N notify_done -hold_jid ${last_job} | cut -f3 -d' ')