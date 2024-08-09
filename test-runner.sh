#!/usr/bin/env bash
set -e
#
# test-runner.sh runs a the entire ARTIC field bioinformatics pipeline using a small set
# of data (the Mayinga barcode from an Ebola amplicon library sequenced on a flongle).
#
# full data available: http://artic.s3.climb.ac.uk/run-folders/EBOV_Amplicons_flongle.tar.gz
#
# usage:
#       ./test-runner.sh [medaka|nanopolish]
#
#   specify either medaka or nanopolish to run the respective workflow of the pipeline
#
###########################################################################################
# Setup the data, commands and the testing function.

# data
inputData="./20190830_1509_MN22126_AAQ411_9efc5448_barcoded"
primerSchemes="../test-data/primer-schemes"
primerScheme="IturiEBOV/V1"
prefix="ebov-mayinga"
barcode="03"
threads=2
downloadCmd="wget https://loman-labz-public-datasets.s3.climb.ac.uk/EBOV_Amplicons_flongle_barcoded.tar.gz"
extractCmd="tar -vxzf EBOV_Amplicons_flongle_barcoded.tar.gz"


guppyplexCmd="artic guppyplex \
        --min-length 400 \
        --max-length 800 \
        --prefix ${prefix} \
        --directory ./${inputData}/fastq_pass/barcode${barcode} \
        --output ${prefix}_guppyplex_fastq_pass-NB${barcode}.fastq"

## medaka workflow specific
minionCmd_m="artic minion \
            --normalise 200 \
            --threads ${threads} \
            --scheme-directory ${primerSchemes} \
            --read-file ${prefix}_guppyplex_fastq_pass-NB${barcode}.fastq \
            --model r941_min_high_g351 \
            ${primerScheme} \
            ${prefix}"

# clair3 workflow specific
minionCmd_c="artic minion \
            --normalise 200 \
            --threads ${threads} \
            --scheme-directory ${primerSchemes} \
            --read-file ${prefix}_guppyplex_fastq_pass-NB${barcode}.fastq \
            --model r941_prom_hac_g360+g422 \
            ${primerScheme} \
            ${prefix}"

# colours
NC='\033[0m'
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'

# cmdTester is function to run a command and check for failure
function cmdTester {
    echo "###########################################################################################"
    echo -e "${BLUE}Running:${NC} $*"
    echo
    "$@"
    echo
    local status=$?
    if [ $status -ne 0 ]
    then
        echo -e "${RED}FAIL${NC}" >&2
        exit $status
    else
        echo -e "${GREEN}PASS${NC}" >&2
    fi
    echo
    return
}

###########################################################################################
# Run the tests.

# check that nanopolish or medaka is specified
if [ "$1" == "nanopolish" ] || [ "$1" == "medaka" ]; then
    echo -e "${BLUE}Starting tests...${NC}"
    echo -e "${BLUE} - using the $1 workflow${NC}"
    echo
else
    echo "please specify medaka or nanopolish"
    echo "./test-runner.sh [medaka|nanopolish]"
    exit 1
fi

# setup a tmp directory to work in
mkdir tmp && cd tmp || exit

# download the data
echo "downloading the test data..."
cmdTester $downloadCmd
cmdTester $extractCmd

# run the correct workflow
echo "running the pipeline..."
if [ "$1" == "medaka" ]
then

    # collect the reads
    cmdTester $gatherCmd

    # run the core pipeline with medaka
    cmdTester $minionCmd_m
else
    # guppyplex the reads
    cmdTester $guppyplexCmd

    # run the core pipeline with clair3
    cmdTester $minionCmd_c
fi

###########################################################################################
# Check the output and clean up.

# check output created
echo -e "${BLUE}Checking output...${NC}"
if test -f "${prefix}.consensus.fasta"
then
    echo -e "${GREEN} - consensus found${NC}"
else
    echo -e "${RED} - no consensus found${NC}"
    exit 1
fi

# cleanup
cd .. && rm -r tmp
echo -e "${BLUE}Done.${NC}"