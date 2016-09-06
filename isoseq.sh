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

declare -xr PACKAGE_NAME="pacbio_isoseq_pipeline"
declare -xr CONTACT_EMAILS="bhan@pacb.com"

######
# OS #
######
if [[ `uname` == "Darwin" ]]; then
    echo "Mac is currently no supported" >&2 && exit 1;
fi

##########
# Config #
##########
declare -xr CALLED_PROGRAM_NAME=$(basename $0)
declare -xr PROGRAM_NAME=$(readlink -f ${0})
declare -xr PIPELINE_DIRECTORY=$(dirname ${PROGRAM_NAME})
declare -x  ANNOTATION_DIR=${PIPELINE_DIRECTORY}/annotation
declare -xr INDEX_DIR=${ANNOTATION_DIR}/index
declare -xr MYBIN=${PIPELINE_DIRECTORY}/bin
declare -xr PATH=${MYBIN}:$PATH

set -eu -o pipefail

################
# load modules #
################
. ${PIPELINE_DIRECTORY}/bin/color.sh
. ${PIPELINE_DIRECTORY}/bin/functions.sh
. ${PIPELINE_DIRECTORY}/config/config.sh

########
# Args #
########
if [[ $# -lt 1 ]]; then usage && exit 1; fi
declare SUBPROGRAM=$(echo ${1} | tr '[A-Z]' '[a-z]')
case $SUBPROGRAM in
  all)
    shift && bash _all.sh "$@" ;;
  *)
    echo2 "unrecognized option \"${1}\"! \nplease type \"${PACKAGE_NAME}\" without options to see all the options and usage." error;;
esac