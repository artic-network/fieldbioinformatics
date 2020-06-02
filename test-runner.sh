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
#   specify either medaka or nanopolish to run the respective workflow of the pipeline.
#
#   the medaka workflow tests guppyplex by tagging it after the original gather and demux
#   steps (in real-world cases, guppyplex would feed on guppy output instead).
#
###########################################################################################
# Setup the data, commands and the testing function.

# data
inputData="./20190830_1509_MN22126_AAQ411_9efc5448"
primerSchemes="../test-data/primer-schemes"
primerScheme="IturiEBOV/V1"
prefix="ebov-mayinga"
barcode="03"
threads=2
downloadCmd="wget http://artic.s3.climb.ac.uk/run-folders/EBOV_Amplicons_flongle.tar.gz"
extractCmd="tar -vxzf EBOV_Amplicons_flongle.tar.gz"

# pipeline commands
## nanopolish workflow specific
gatherCmd_n="artic gather \
        --min-length 400 \
        --max-length 800 \
        --prefix ${prefix} \
        --directory ${inputData} \
        --fast5-directory ${inputData}/fast5_pass"

demuxCmd_n="artic demultiplex \
            --threads ${threads} \
            ${prefix}_fastq_pass.fastq"

minionCmd_n="artic minion \
                --normalise 200 \
                --threads ${threads} \
                --scheme-directory ${primerSchemes} \
                --read-file ${prefix}_fastq_pass-NB${barcode}.fastq \
                --fast5-directory ${inputData}/fast5_pass \
                --sequencing-summary ${inputData}/lab-on-an-ssd_20190830_160932_AAQ411_minion_sequencing_run_EBOV_Amplicons_flongle_sequencing_summary.txt \
                ${primerScheme} \
                ${prefix}"

## medaka workflow specific
gatherCmd_m="artic gather \
        --min-length 400 \
        --max-length 800 \
        --prefix ${prefix} \
        --directory ${inputData} \
        --no-fast5s"

demuxCmd_m="artic demultiplex \
            --threads ${threads} \
            ${prefix}_fastq_pass.fastq"

guppyplexCmd_m="artic guppyplex \
        --min-length 400 \
        --max-length 800 \
        --prefix ${prefix} \
        --directory ./ \
        --output ${prefix}_guppyplex_fastq_pass-NB${barcode}.fastq"

minionCmd_m="artic minion \
            --normalise 200 \
            --threads ${threads} \
            --scheme-directory ${primerSchemes} \
            --read-file ${prefix}_guppyplex_fastq_pass-NB${barcode}.fastq \
            --medaka \
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
if [ "$1" == "nanopolish" ]
then

    # collect the reads
    cmdTester $gatherCmd_n

    # demultiplex
    cmdTester $demuxCmd_n

    # run the core pipeline with nanopolish
    cmdTester $minionCmd_n
else

    # collect the reads
    cmdTester $gatherCmd_m

    # demultiplex
    cmdTester $demuxCmd_m

    # guppyplex
    cmdTester $guppyplexCmd_m

    # run the core pipeline with medaka
    cmdTester $minionCmd_m
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