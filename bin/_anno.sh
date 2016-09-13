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
declare -r MODULE_NAME=Anno
declare -r MODULE_VERSION=0.0.2.160911

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
=====================================
Prepare Annotation Files for a Genome
=====================================

This module prepare annotation files for other modules to use.
It currently have two modes, 
1. obtaining everything from UCSC
2. generating a new assembly from two existing assemblies

OPTIONS:
        -h      Show usage
        -v      Print version

${REQUIRED}[ required ]
        -g      Name of the assembly, for example hg19
${OPTIONAL}[ optional ]
        -s      Sub assembly names, create new assembly from existing assemblies, like hg19 + SIRV
        -t      Number of CPU to use in each job. Default: ${DEFAULT_NUM_THREADS}
${ADVANCED}[ advanced ]
        -D      use debug mode (bash -x). Default: off
EOF
echo -e "${FONT_COLOR_RESET}"
}

##########
# Config #
##########
declare -a REQUIRED_PROGRAMS=( 'rsync' 'twoBitToFa' 'smrtshell' 'readlink' 'find' \
                            'gmap' 'gmap_build' \
                            'faSize' 'bedToGenePred' 'genePredToGtf' \
                            )
declare -a SUB_ASSEMLIES=()

#############################
# ARGS reading and checking #
#############################
while getopts "hvg:s:t:D" OPTION; do
    case $OPTION in
        h)  usage && exit 0 ;;
        v)  echo ${PACKAGE_NAME}::${MODULE_NAME} v${MODULE_VERSION} && exit 0;;
        s)  SUB_ASSEMLIES+=(${OPTARG});;
        g)  declare -x AssemblyName=$(echo "${OPTARG}" | tr '[A-Z]' '[a-z]');;
        t)  declare -i Threads=${OPTARG};; 
        D)  set -x;;
        *)  usage && exit 1;;
    esac
done

if [[ -z ${AssemblyName} ]]; then echo2 "You have to provide an assembly name with -g option" error; fi
if [[ -z $Threads ]]; then declare -i Threads=${DEFAULT_NUM_THREADS}; fi
if [[ -z "${ANNOTATION_DIR}" ]]; then echo2 "Unspecified ANNOTATION_DIR, please check your config.sh file" error; fi
if [[ ! -d "${ANNOTATION_DIR}" ]]; then echo2 "${ANNOTATION_DIR} does not exist" warning; mkdir -p ${ANNOTATION_DIR} || echo2 "cannot create folder ${ANNOTATION_DIR}" error; fi
if [[ ! -w "${ANNOTATION_DIR}" ]]; then echo2 "${ANNOTATION_DIR} is not writable" error; fi
cd "${ANNOTATION_DIR}" || echo2 "Don't have permission to access folder $ANNOTATION_DIR" error

if [[ -z ${SUB_ASSEMLIES} ]]; then REQUIRED_PROGRAMS+=('mysql'); fi

for program in "${REQUIRED_PROGRAMS[@]}"; do binCheck $program; done

#########
# Begin #
#########
declare GenomeSequence=${ANNOTATION_DIR}/${AssemblyName}.fa
declare GenomeSize=${ANNOTATION_DIR}/${AssemblyName}.sizes
declare TranscriptomeBed=${ANNOTATION_DIR}/${AssemblyName}.genes.bed12
declare TranscriptomeGp=${ANNOTATION_DIR}/${AssemblyName}.genes.gp
declare TranscriptomeGtf=${ANNOTATION_DIR}/${AssemblyName}.genes.gtf
declare RefFlat=${ANNOTATION_DIR}/${AssemblyName}.refFlat.txt
declare gmapindex=${GMAP_INDEX_DIR}/${genome}

if [[ -z ${SUB_ASSEMLIES} ]]; then 
    # no sub assemlby
    echo2 "No sub assemblies are provided, will obtain everything from UCSC" warning
    
    # genome sequence
    if [[ ! -f ${GenomeSequence} ]]; then
        echo2 "Obtain genome sequence from UCSC"
        echo2 "Try compressed gzip"
        rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/${AssemblyName}/bigZips/${AssemblyName}.fa.gz ./ \
        && gunzip ${AssemblyName}.fa.gz
        if [[ ! -f ${GenomeSequence} ]]; then
            echo2 "Try 2bits"
            rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/${AssemblyName}/bigZips/${AssemblyName}.2bit ./ \
        && twoBitToFa ${AssemblyName}.2bit ${GenomeSequence}
        fi
        
        if [[ ! -f ${GenomeSequence} ]]; then echo2 "Cannot obtain ${AssemblyName}.fa ... " error; fi
    fi

    # genome size
    if [[ ! -s ${GenomeSize} ]]; then
        echo2 "Generate size file"
        faSize -detailed -tab ${GenomeSequence} > ${GenomeSize} \
        || echo2 "Failed to generate size file" error
    fi
    
    # transcriptome
    if [[ ! -f ${TranscriptomeGtf} ]]; then
        echo2 "Trying to obtain gtf file from UCSC"
        if [[ ! -s ${HOME}/.hg.conf ]]; then
            cat > ${HOME}/.hg.conf << EOF
db.host=genome-mysql.cse.ucsc.edu
db.user=genomep
db.password=password
central.db=hgcentral
EOF
            chmod 600 ${HOME}/.hg.conf
            genePredToGtf ${AssemblyName} refFlat ${TranscriptomeGtf}
    else 
        echo2 "${TranscriptomeGtf} has already existed" warning
    fi 
    
    if [[ ! -s ${TranscriptomeGtf} ]]; then
        echo2 "Failed to obtain GTF file" error
    fi

    # refflat
    if [[ ! -s ${RefFlat} ]]; then
        mysql \
            -h genome-mysql.cse.ucsc.edu \
            -u genome \
            -D ${AssemblyName} \
            -N \
            -A \
            -e 'SELECT * FROM refFlat' \
            > ${RefFlat} \
        || echo2 "Cannot obtain refFlat file from UCSC"
    fi

else # combine assemblies
    echo2 "sub assemblies ${SUB_ASSEMLIES[@]} are provided, will try to combine them" warning

    # merge genome sequence
    declare -a subgenomefas=()
    declare -a subgenomegtfs=()
    declare -a refflatfiles=()
    for subgenome in "${SUB_ASSEMLIES[@]}"; do
        declare subgenomefa=${ANNOTATION_DIR}/${subgenome}.fa
        if [[ ! -s ${subgenomefa} ]]; then 
            echo2 "file ${subgenomefa} does not exist or it is empty" error; 
        fi

        declare subgenomegtf=${ANNOTATION_DIR}/${subgenome}.genes.gtf
        if [[ ! -s ${subgenomegtf} ]]; then 
            echo2 "file ${subgenomegtf} does not exist or it is empty" error; 
        fi
        subgenomegtfs+=(${subgenomegtf})

        declare refflatfile=${ANNOTATION_DIR}/${subgenome}.refFlat.txt
        if [[ ! -s ${refflatfile} ]]; then
            echo2 "file ${refflatfile} does not exist or it is empty" error; 
        fi
            refflatfiles+=(${refflatfile})
    done
    
    cat ${subgenomefas[@]} > ${GenomeSequence}

    # genome size
    if [[ ! -s ${GenomeSize} ]]; then
        echo2 "Generate size file"
        faSize -detailed -tab ${GenomeSequence} > ${GenomeSize} \
        || echo2 "Failed to generate size file" error
    fi

    # merge transcriptome
    cat ${subgenomegtfs[@]} > ${TranscriptomeGtf}

    # merge refflat
    cat ${refflatfiles[@]} > ${RefFlat}
fi

# build gmap index
if ! gmapIndexCheck ${GMAP_INDEX_DIR} $AssemblyName; then
    echo2 "gmap index is not found"
    declare jobfile=$(mktemp)
    cat > ${jobfile} << EOF
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mkdir -p ${GMAP_INDEX_DIR}/${AssemblyName} \
&& gmap_build -d ${AssemblyName} -D ${GMAP_INDEX_DIR} ${GenomeSequence}
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
    declare gmap_build_jobid=$(${SUBMIT_CMD} -o log -e log -N job_${AssemblyName}.gmapbuild < ${jobfile} | cut -f3 -d' ')
    echo2 "Submit job ${gmap_build_jobid} to build gmap index"
fi
