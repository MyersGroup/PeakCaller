#!/bin/bash
set -euo pipefail

# # # # Example Usage # # # #

# sh MAPeakCaller.sh \
# --outdir 538916/ \
# --chrsizes hg38.sizes \
# --name 223180_vs_221156 \
# -a Fragment_Position_538916_223180.sorted.bed.PR1 \
# -b Fragment_Position_538916_223180.sorted.bed.PR2 \
# -i Fragment_Position_538916_221156.sorted.bed \
# --autosomes 22 \
# --pthresh 0.000001 \
# --peakminsep 250



POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -o|--outdir)
    OUTPATH="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--chrsizes)
    CHRSIZES="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--name)
    NAME="$2"
    shift # past argument
    shift # past value
    ;;
    -a)
    A="$2"
    shift # past argument
    shift # past value
    ;;
    -b)
    B="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--input)
    I="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--autosomes)
    AUTOSOMES="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--pthresh)
    PTHRESH="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--peakminsep)
    MINPEAKSEP="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters


Rscript EstimateConstants.R \
$OUTPATH \
$CHRSIZES \
$NAME \
$A \
$B \
$I \
$AUTOSOMES


Rscript DeNovoPeakCalling-SingleBase.R \
$OUTPATH \
$CHRSIZES \
$NAME \
$A \
$B \
$I \
$AUTOSOMES \
$PTHRESH \
$MINPEAKSEP \
