# author: Bo Han (bhan@pacb.com)

function usage {
    export BOLD=`echo -ne '\e[1m'`;
    export UNDERLINE=`echo -ne '\e[4m'`;
    export BLINK=`echo -ne '\e[5m'`;
    export RESET=`echo -ne '\e[0m'`;

    cat << EOF
==================
${RESET}${BOLD}${PACKAGE_NAME}${RESET}
==================
This is a pipeline to run PacBio Iso-Seq pipeline from RS II cells (primary analysis) to a final report.

[ all ]
    CCS + classify + isoaux + report

EOF
}


export ISO_8601='%Y-%m-%d %H:%M:%S %Z'
function echo2 {
case $2 in
    error|e)
        echo -e $FONT_COLOR_RED"[`date "+$ISO_8601"`] Error: $1${FONT_COLOR_RESET}" >&2 && exit 1;;
    warning|w)
        echo -e $FONT_COLOR_MAGENTA"[`date "+$ISO_8601"`] Warning: $1${FONT_COLOR_RESET}" >&2;;
    *)
        echo -e $FONT_COLOR_GREEN"[`date "+$ISO_8601"`] $1${FONT_COLOR_RESET}" >&2;;
esac
}


function binCheck {
    program_path=`which $1 2>/dev/null`
    [[ $? -ne 0 ]] && echo2 "Required program \"$1\" is not available" error
    [[ $VERBOSE -gt 0 ]] && echo2 "$1: $program_path"
}

function fileCheck {
    [[ ! -f "${1}" ]] && echo2 "Required file \"${1}\" doesn't exist" error
}

function dirOrFileCheck {
    [[ ! -d "${1}" ]] && [[ ! -f "${1}"  ]] && echo2 "Required file or directory \"${1}\" doesn't exist" error
}

function isFile {
    [[ ! -f "${1}" ]] && return 1
    return 0
}


function isDir {
    [[ ! -d "${1}" ]] && return 1
    return 0
}

function indexOf {
    local element
    local i=0
    for element in "${@:2}"; do
        [[ "$element" == "$1" ]] && echo $i && return 0
        let i+=1
    done
    echo -1 && return 1
}

function validateRSII {
    # TODO: need more robust validation
    [[ ! -d $1/Analysis_Results ]] && return 1
    return 0
}

function gmapIndexCheck {
    [[ $# -ne 2 ]] && return 255;
    [[ ! -f ${1}/${2}/${2}.chromosome ]]     && return 1;
    [[ ! -f ${1}/${2}/${2}.chromosome.iit ]] && return 1;
    [[ ! -f ${1}/${2}/${2}.chrsubset ]]      && return 1;
    [[ ! -f ${1}/${2}/${2}.contig ]]         && return 1;
    [[ ! -f ${1}/${2}/${2}.contig.iit ]]     && return 1;
    [[ ! -f ${1}/${2}/${2}.genomecomp ]]     && return 1;
    return 0;
}
	
export -f usage
export -f echo2
export -f binCheck
export -f fileCheck
export -f dirOrFileCheck
export -f isFile
export -f isDir
export -f indexOf
export -f validateRSII
export -f gmapIndexCheck