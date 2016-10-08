#!/bin/bash 

# this scripts upload file to ftp and generate UCSC genome browser custom track automatically

# const
declare -r DefaultWatsonColorRgb="0,0,255"
declare -r DefaultCrickColorRgb="255,0,0"
declare -r DefaultWidthInPixel=25
declare -r DefaultVisibility=full

# argv
if [[ $# -lt 4 ]]; then
    echo "usage:"
    echo -e "\t$0 username password ftpaddress dirname file1 file2 ... filen"
    exit 1
fi

# begin
declare -r username=${2}
declare -r password=${3}
declare -r ftpname=${1}
declare -r dirname=${4}
shift;shift;shift;shift;

# generate the public folder if not already exist, this command return 1 if the folder has already existed
lftp -u ${username},${password} ${ftpname} -e "mkdir -p ${dirname}; bye" &>/dev/null

for f in $@; do
    if ! lftp -u ${username},${password} ${ftpname} -e "put -O ${dirname} $f; bye" 1>&2; then
        echo "unable to use this ftp, please update username and password in $0"
        exit 2
    fi
    # generate track
    declare suffix=${f##*.}
    declare filename=$(basename $f)
    case $suffix in 
    bigWig|bw)
        declare Color=${DefaultCrickColorRgb} # default Crick red
        if [[ $filename == *Forward* ]]; then Color=${DefaultWatsonColorRgb}; fi # change color to blue for Watson strand
        echo "track name=${f%.*} type=bigWig visibility=${DefaultVisibility} color=$Color maxHeightPixels=${DefaultWidthInPixel}:${DefaultWidthInPixel}:8 bigDataUrl=ftp://${username}:${password}@${ftpname}/${dirname}/${filename}"
        ;;

    bam)
        declare idx=${f}.bai
        if [[ ! -f ${idx} ]]; then samtools index ${f}; fi
        lftp -u ${username},${password} ${ftpname} -e "put -O ${dirname} $idx; bye" 1>&2
        echo "track name=${f%.*} type=bam visibility=${DefaultVisibility} maxHeightPixels=${DefaultWidthInPixel}:${DefaultWidthInPixel}:8 bigDataUrl=ftp://${username}:${password}@${ftpname}/${dirname}/${filename}"
        ;;

    bigBed|bb)
        echo "track name=${f%.*} type=bigBed visibility=${DefaultVisibility} itemRgb=\"On\" maxHeightPixels=${DefaultWidthInPixel}:${DefaultWidthInPixel}:8 bigDataUrl=ftp://${username}:${password}@${ftpname}/${dirname}/${filename}"
        ;;

    *)  
        echo "Unrecognized format; cannot generate track but the file has been uploaded to:";
        echo "ftp://${username}:${password}@${ftpname}/${dirname}/${filename}"
    ;;
    esac
done
