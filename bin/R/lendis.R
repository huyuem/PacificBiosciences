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

source (paste(Sys.getenv("PIPELINE_DIRECTORY"),"/bin/R/lib.R",sep=""))

pkgTest("readr")
pkgTest("ggplot2")
pkgTest("dplyr")
pkgTest("tidyr")
pkgTest("ggthemes")

argv  = commandArgs (TRUE)
file = argv[1]
outpdf = argv[2]

data = read_tsv(file, T)
num_colors = length(levels(as.factor(data$treatment)))

pdf(outpdf)
ggplot(data) + 
    geom_line(aes(size, y=..density.., colour=treatment), lwd = 1.25, stat="density") +
    scale_color_brewer(palette=(ifelse(num_colors < 9, "Set1", "Set3"))) +
    facet_grid(sizebin ~ tissue) +
    xlim(0, median(data$size) * 3) +
    xlab('length(nt) of flnc') +
    ylab('density') +
    theme_bw()

ggplot(data) + 
    geom_line(aes(size, y=..count.., colour=treatment), lwd = 1.25, stat="density") +
    scale_color_brewer(palette=(ifelse(num_colors < 9, "Set1", "Set3"))) +
    facet_grid(sizebin ~ tissue) +
    xlim(0, median(data$size) * 3) +
    xlab('length(nt) of flnc') +
    ylab('count') +
    theme_bw()

dev.off()