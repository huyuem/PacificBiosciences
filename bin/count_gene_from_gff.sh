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
set -x
declare info="${1}"
declare refgtf="${2}"
shift 
shift
source ${info}
mkdir -p jobout/gffcompareOut
if [[ ! -f jobout/gffcompareOut/_Done ]]; then
    gffcompare -r ${refgtf} -G -o jobout/gffcompareOut/_ "$@" \
 && touch jobout/gffcompareOut/_Done
fi
declare -a geneCountToMerge=()
declare -a mrnaCountToMerge=()
echo -e "genome\ttissue\ttreatment\tsize_bin\tname\tcounts" > table/gene.counts.melted.tsv \
&& echo -e "genome\ttissue\ttreatment\tsize_bin\tname\tcounts" > table/mRNA.counts.melted.tsv \
&& echo -ne "name" > table/gene.counts.tsv \
&& echo -ne "name" > table/mRNA.counts.tsv

for i in $(seq -s " " 0 $((SampleSize-1))); do
    declare genome=${Genomes[$i]}
    declare tissue=${Tissues[$i]}
    declare treatment=${Treatments[$i]}
    declare sizebin=${SizeBins[$i]}
    declare -i j=$((i+1))
    declare gtf=${!j}
    declare gtfbase=$(basename $gtf)
    declare prefix=${gtfbase%.isoseq*}
    # fix special chars which confuse CSS selector
    for r in '+' '.'; do 
        prefix=${prefix//${r}/_};
    done

    echo -ne "\t${prefix}" >> table/gene.counts.tsv # convert . to _ because . is mistaken as class in js
    echo -ne "\t${prefix}" >> table/mRNA.counts.tsv 

    python ${MYBIN}/count_gene_from_refmap.py jobout/gffcompareOut/_.${gtfbase}.refmap gene > table/${prefix}.gene.count \
 && awk -v genome=${genome} -v tissue=${tissue} -v treatment=${treatment} -v sizebin=${sizebin} \
        'BEGIN{FS=OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%f\n", genome, tissue, treatment, sizebin, $1, $2}' table/${prefix}.gene.count >> table/gene.counts.melted.tsv \
 && python ${MYBIN}/count_gene_from_refmap.py jobout/gffcompareOut/_.${gtfbase}.refmap mRNA > table/${prefix}.mRNA.count \
 && awk -v genome=${genome} -v tissue=${tissue} -v treatment=${treatment} -v sizebin=${sizebin} \
        'BEGIN{FS=OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%f\n", genome, tissue, treatment, sizebin, $1, $2}' table/${prefix}.mRNA.count >> table/mRNA.counts.melted.tsv

    geneCountToMerge+=("table/${prefix}.gene.count")
    mrnaCountToMerge+=("table/${prefix}.mRNA.count")
done

echo >> table/gene.counts.tsv
colmerge -i ${geneCountToMerge[@]} -c 0 -t 1 -d "0.0" >> table/gene.counts.tsv

echo >> table/mRNA.counts.tsv
colmerge -i ${mrnaCountToMerge[@]} -c 0 -t 1 -d "0.0" >> table/mRNA.counts.tsv

# make HTML
cat ${PIPELINE_DIRECTORY}/html_templates/scatter_plot_gene1.html > html/gene_abundance.html
cat table/gene.counts.tsv >> html/gene_abundance.html
cat ${PIPELINE_DIRECTORY}/html_templates/scatter_plot_gene2.html >> html/gene_abundance.html

# mRNA needs gene name and length
mrna_size_from_gff $refgtf > table/mRNA_gene.sizes.tsv
cat ${PIPELINE_DIRECTORY}/html_templates/scatter_plot_mRNA1.html > html/mRNA_abundance.html
echo -ne "transcript_name\tgene_name\tlength\t" >>  html/mRNA_abundance.html
head -1 table/mRNA.counts.tsv | cut -f1 --complement >> html/mRNA_abundance.html 
gawk 'BEGIN{FS=OFS="\t"}{if(ARGIND==1) {g[$3]=$2; l[$3]=$4} else {if(g[$1]){ $1 = $1"\t"(g[$1]?g[$1]:"unknown")"\t"(l[$1]?l[$1]:"unknown"); print $0 }}}' table/mRNA_gene.sizes.tsv table/mRNA.counts.tsv >> html/mRNA_abundance.html
cat ${PIPELINE_DIRECTORY}/html_templates/scatter_plot_mRNA2.html >> html/mRNA_abundance.html

Rscript ${MYBIN}/R/abundance.R table/gene.counts.melted.tsv pdf/gene_counts.pdf
Rscript ${MYBIN}/R/abundance.R table/mRNA.counts.melted.tsv pdf/mRNA_counts.pdf