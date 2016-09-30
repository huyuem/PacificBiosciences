#!/bin/bash

if [[ $# -lt 2 ]]; then
	echo "usage: $0 input.fa index_dir index_name threads=8"
	exit 1
fi

inputfile=${1}
inputname=$(basename ${inputfile})
index_dir=${2}
index_name=${3}
threads=${4:-8}

function _gmapIndexCheck {
    [[ $# -ne 2 ]] && return 255;
    [[ ! -f ${1}/${2}/${2}.chromosome ]]     && return 1;
    [[ ! -f ${1}/${2}/${2}.chromosome.iit ]] && return 1;
    [[ ! -f ${1}/${2}/${2}.chrsubset ]]      && return 1;
    [[ ! -f ${1}/${2}/${2}.contig ]]         && return 1;
    [[ ! -f ${1}/${2}/${2}.contig.iit ]]     && return 1;
    [[ ! -f ${1}/${2}/${2}.genomecomp ]]     && return 1;
    return 0;
}
	
if ! _gmapIndexCheck $index_dir $index_name; then
	echo "failed to detect gmap index at $index_dir/$index_name"
	exit 1
fi

gmap \
	${GmapOption} \
	-t ${threads} \
	-D ${index_dir} \
	-d ${index_name} \
	-f samse \
	--suboptimal-score=1 \
	-n 0 \
	-z sense_force \
	${inputfile} \
| samtools view -bS - \
| samtools sort -f - bam/${inputname%f[aq]*}${index_name}.sorted.bam \
&& samtools index bam/${inputname%f[aq]*}${index_name}.sorted.bam