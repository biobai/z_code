#!/bin/bash
#set -e
#set -o pipefail
version=0.1
############################################################
#  Program: extract coordinates by gene names from gtf file
#  Version: 0.1
#  Author: Youhuang Bai
###################################s#########################
LC_ALL=C
function usage() {
    echo "
usage:   q2c.sh -i query.gene -g gtf.file -o output.file
"
}
while getopts ":hi:o:g:" OPTION
do
    case "${OPTION}" in
            h)
                usage
                exit 1
                ;;
            i)
                queryGene="$OPTARG"
                ;;
            o)
                outputFile="$OPTARG"
                ;;
            g)
                gtfFile="$OPTARG"
                ;;
        esac
done
sed -i '/^[[:space:]]*$/d' $queryGene
cat $queryGene |xargs -n1 -I {} grep -w \"{}\" $gtfFile >hg19.gtf
cut -f9 hg19.gtf|cut -d ";" -f3|cut -f2 -d '"' >hg19.id
sed -i '/^[[:space:]]*$/d' hg19.id
n1=`wc -l <hg19.gtf`
n2=`wc -l <$queryGene`
if [ $n1 != $n2 ]; then
   diff -y $queryGene hg19.id >error.txt
   exit
else
   paste $queryGene <(awk 'BEGIN{FS="\t"}{print $1":"$4+1"-"$5+1}'  hg19.gtf) >$outputFile
   rm hg19.gtf hg19.id
fi
