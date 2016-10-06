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
declare -r MODULE_NAME=Gather
declare -r MODULE_VERSION=0.0.1.161006

########################
# Variable declaration #
########################
declare    OutputDir=""
declare -a Cells=()

#########
# Usage #
#########
function usage {
cat << EOF
${PACKAGE_NAME}::${MODULE_NAME}
====================================================
Gather Multiple RS II Cells and Make a new directory
====================================================

This script make a new directory and gather all the RS II cells (biological
and technical replicates) so that a single directory can be used in the 
other pipeline.
usage:
    isoseq.sh gather \ 
        -o /path/to/new_cell \ 
        -i /path/to/cell1 \ 
        -i /path/to/cell2 \ 
        -i /path/to/cell3

OPTIONS:
        -h      Show usage
        -v      Print version

${REQUIRED}[ required ]
        -o      Output directory
        -i      Absolute path to individual cells
${ADVANCED}[ advanced ]
        -D      use debug mode (bash -x). Default: off
EOF
echo -e "${FONT_COLOR_RESET}"
}

##########
# Config #
##########
declare -a REQUIRED_PROGRAMS=('find' 'ln' 'mkdir' 'touch')

#############################
# ARGS reading and checking #
#############################
while getopts "hvi:o:D" OPTION; do
    case $OPTION in
        h)  usage && exit 0 ;;
        v)  echo ${PACKAGE_NAME}::${MODULE_NAME} v${MODULE_VERSION} && exit 0;;
        o)  OutputDir=${OPTARG};;
        i)  if [[ -d ${OPTARG} ]]; then 
                Cells+=( $(readlink -f ${OPTARG}) )
            else 
                echo2 "directory ${OPTARG} does not exist" error
            fi
            ;;
        D)  set -x;;
        *)  usage && exit 1;;
    esac
done
# make sure the user has input -o
if [[ -z $OutputDir ]]; then echo2 "Please specify the output directory to store the symble links" error; fi
# make sure the user has input -i
if [[ -z $Cells ]]; then echo2 "Please specify the input cell paths" error; fi
# create outputdir
if ! mkdir -p ${OutputDir}/Analysis_Results; then echo2 "Cannot create directory ${OutputDir}/Analysis_Results" error; fi
# access outputdir
if ! cd ${OutputDir}/Analysis_Results; then echo2 "Cannot access directory ${OutputDir}/Analysis_Results" error; fi
# test writting permission
if touch a &>/dev/null; then rm a; else echo2 "Cannot crate file in ${OutputDir}" error; fi

# gathering
for cell in "${Cells[@]}"; do
    if ! validateRSII "${cell}"; then echo2 "Cell ${cell} doesn't look like a RS II run" error; fi
    for h5 in $(find ${cell} -name "*h5"); do
        if ! ln -s ${h5}; then echo2 "cannot create symbol link to ${h5} in ${PWD}" error; fi
    done
done

echo2 "Done"