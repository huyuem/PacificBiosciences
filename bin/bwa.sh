#!/bin/bash

if [[ $# -lt 2 ]]; then
	echo "usage: $0 input.fa index_name threads=8"
	exit 1
fi

inputfile=${1}
inputname=$(basename ${inputfile})
index_name=${2}
threads=${3:-8}

function _bwaIndexCheck {
    [[ $# -ne 2 ]] && return 255;
    [[ ! -f ${1}.bwa ]] && return 1;
    [[ ! -f ${1}.pac ]] && return 1;
    [[ ! -f ${1}.ann ]] && return 1;
    [[ ! -f ${1}.amb ]] && return 1;
    [[ ! -f ${1}.sa  ]] && return 1;
    return 0;
}
	
if ! _bwaIndexCheck $index_name; then
	echo "failed to detect bwa index at $index_name"
	exit 1
fi

bwa mem \
    -t $threads \
    -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0 \
    ${index_name} \
    ${inputfile} \
| samtools view -bS - \
| samtools sort -f - bam/${inputname%f[aq]*}${index_name}.sorted.bam \
&& samtools index bam/${inputname%f[aq]*}${index_name}.sorted.bam
