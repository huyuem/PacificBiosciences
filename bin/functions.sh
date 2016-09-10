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

export ISO_8601='%Y-%m-%d %H:%M:%S %Z'
function echo2 {
    declare p=${2:-}
    case $p in
        error|e)
            echo -e $FONT_COLOR_RED"[`date "+$ISO_8601"`] Error: $1${FONT_COLOR_RESET}" >&2 && exit 1;;
        warning|w)
            echo -e $FONT_COLOR_MAGENTA"[`date "+$ISO_8601"`] Warning: $1${FONT_COLOR_RESET}" >&2;;
        *)
            echo -e $FONT_COLOR_GREEN"[`date "+$ISO_8601"`] $1${FONT_COLOR_RESET}" >&2;;
    esac
}
export -f echo2

function binCheck {
    local program_path=$(which $1 2>/dev/null)
    if [[ $? -ne 0 ]]; then echo2 "Required program \"$1\" is not available" error; fi
}
export -f binCheck

function fileCheck {
    if [[ ! -f "${1}" ]]; then echo2 "Required file \"${1}\" doesn't exist" error; fi
}
export -f fileCheck

function dirOrFileCheck {
    if [[ ! -d "${1}" && ! -f "${1}"  ]]; then echo2 "Required file or directory \"${1}\" doesn't exist" error; fi
}
export -f dirOrFileCheck

function isFile {
    if [[ ! -f "${1}" ]]; then return 1; fi
    return 0
}
export -f isFile

function isDir {
    if [[ ! -d "${1}" ]]; then return 1; fi
    return 0
}
export -f isDir

function indexOf {
    local element
    local i=0
    for element in "${@:2}"; do
        if [[ "$element" == "$1" ]]; then echo $i ; return 0; fi
        let i+=1
    done
    echo -1 && return 1
}
export -f indexOf

function validateRSII {
    # TODO: need more robust validation
    if [[ ! -d $1/Analysis_Results ]]; then return 1; fi
    return 0
}
export -f validateRSII

function gmapIndexCheck {
    if [[ $# -ne 2 ]]                          ;then return 255; fi
    if [[ ! -f ${1}/${2}/${2}.chromosome ]]    ;then return 1; fi
    if [[ ! -f ${1}/${2}/${2}.chromosome.iit ]];then return 1; fi
    if [[ ! -f ${1}/${2}/${2}.chrsubset ]]     ;then return 1; fi
    if [[ ! -f ${1}/${2}/${2}.contig ]]        ;then return 1; fi
    if [[ ! -f ${1}/${2}/${2}.contig.iit ]]    ;then return 1; fi
    if [[ ! -f ${1}/${2}/${2}.genomecomp ]]    ;then return 1; fi
    return 0;
}
export -f gmapIndexCheck