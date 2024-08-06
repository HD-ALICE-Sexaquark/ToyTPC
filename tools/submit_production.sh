#!/bin/bash

####### ###### #### ### ## #
#  Submit batch jobs to generate simulations
### ## #

# hard-coded options
export OUTPUT_DIR= # output directory, assigned by default within SIM_DIR, unless you modify it here
export N_EVENTS_PER_RUN=100 # number of events per run (recommended: 100)
export GEANT4_NUM_THREADS=1 # number of threads to run simultaneously during reconstruction (recommended: 4)

export SBATCH_SETUP="--partition=main"
SBATCH_SETUP+=" --cpus-per-task=${GEANT4_NUM_THREADS}"
SBATCH_SETUP+=" --time=02:00:00"
SBATCH_MEMORY=$((2000 / ${GEANT4_NUM_THREADS}))
SBATCH_SETUP+=" --mem-per-cpu=${SBATCH_MEMORY}"

function print_help() {
    echo "SCRIPT: send_production.sh"
    echo "=========================="
    echo "REQUIREMENTS:"
    echo "  * ROOT"
    echo "  * GEANT4"
    echo "USAGE:"
    echo "  ./send_production.sh --run1 <run1> --run2 <run2> --jdir <job_dir> --bz <mag_field>"
    echo "  ./send_production.sh --rn <rn> --jdir <job_dir> --bz <mag_field>"
    echo "                       where:"
    echo "                       <run1>      = run number of the first job (starting point of loop)"
    echo "                       <run2>      = run number of the last job (end of loop)"
    echo "                       <rn>        = specific run numbers, separated by comma"
    echo "                                     for example: --rn 0,1,2"
    echo "                       <job_dir>   = subdirectory (it's a string!) within the output directory, for organizational purposes"
    echo "                       <mag_field> = z-component magnitude of magnetic field (in Tesla)"
    echo "EXAMPLE:"
    echo "  ./send_production.sh --bz 0.2 --run1 0 --run2 10 --jdir 00"
    echo "  ./send_production.sh --bz 0.2 --rn 4,8,15,16 --jdir 01"
}

function process_args() {
    arr=("$@")
    ic=0
    while [[ $ic -le $((${#arr[@]}-1)) ]]; do
        if [[ "${arr[$ic]}" == "--run1" ]]; then
            RUN1=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--run2" ]]; then
            RUN2=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--rn" ]]; then
            RUN_NUMBER_LIST=${arr[$((ic+1))]}
            RUN_NUMBER_ARR=(${RUN_NUMBER_LIST//,/ }) # convert comma-separated list into array of numbers
        elif [[ "${arr[$ic]}" == "--jdir" ]]; then
            export JOB_DIR=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--bz" ]]; then
            export MAG_FIELD=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--help" ]]; then
            print_help
            exit 0
        else
            echo "ERROR: unrecognized argument ${arr[$((ic))]}, use --help to see available commands"
            exit 1
        fi
        ((ic+=2))
    done
}

argArray=("$@")
process_args "${argArray[@]}"

if [[ -z ${RUN1} ]] && [[ -z ${RUN2} ]] && [[ -z ${RUN_NUMBER_LIST} ]]; then
    echo "ERROR: run number options empty, use --help to see available commands"
    exit 1
fi

if [[ -n ${RUN1} || -n ${RUN2} ]] && [[ -n ${RUN_NUMBER_LIST} ]]; then
    echo "ERROR: you must choose between --rn OR --run1 and --run2, not combine them!"
    exit 1
fi

if [[ -z ${JOB_DIR} ]]; then
    echo "ERROR: you must choose between --rn OR --run1 and --run2, not combine them!"
    exit 1
fi

# set paths

export SIM_DIR=
if [[ "${PWD##*/}" == "tools" ]]; then
    SIM_DIR=${PWD}/..
else
    SIM_DIR=${PWD}
fi

SBATCH_SETUP+=" --output=${SIM_DIR}/slurm-logs/slurm-%J.out"
SBATCH_SETUP+=" --error=${SIM_DIR}/slurm-logs/slurm-%J.err"

if [[ -z ${OUTPUT_DIR} ]]; then
    echo "WARNING: OUTPUT_DIR empty, setting default: ${SIM_DIR}/output"
    OUTPUT_DIR=${SIM_DIR}/output
fi

# check for environment

if [[ ! $(geant4-config) ]]; then
    echo "ERROR: make sure to set GEANT4"
    exit 1
fi

if [[ -z ${ROOTSYS} ]]; then
    echo "ERROR: make sure to set ROOT"
    exit 1
fi

# check for binaries

if [[ ! -e ${SIM_DIR}/reconstruction/G4Setup_build/main ]]; then
    echo "ERROR: reconstruction/G4Setup_build/main not found. Maybe forgot to build?"
    exit 1
fi

# if output directory doesn't exist, create it
mkdir -p ${OUTPUT_DIR}
mkdir -p slurm-logs

# fill array of run numbers
if [[ ${#RUN_NUMBER_ARR[@]} -eq 0 ]] && [[ -n ${RUN1} || -n ${RUN2} ]]; then
    RUN_NUMBER_ARR=()
    for ((run=${RUN1}; run <= ${RUN2}; run++)); do
        RUN_NUMBER_ARR+=(${run})
    done
fi

# print info
echo "send_production.sh :: initiating..."
echo "send_production.sh ::"
echo "send_production.sh :: hard-coded parameters:"
echo "send_production.sh :: >> OUTPUT_DIR         = ${OUTPUT_DIR}"
echo "send_production.sh :: >> N_EVENTS_PER_RUN   = ${N_EVENTS_PER_RUN}"
echo "send_production.sh :: >> GEANT4_NUM_THREADS = ${GEANT4_NUM_THREADS}"
echo "send_production.sh ::"
echo "send_production.sh :: chosen options:"
echo "send_production.sh :: >> sim_dir   = ${SIM_DIR}"
echo "send_production.sh :: >> mag_field = ${MAG_FIELD}"
echo "send_production.sh :: >> job_dir   = ${JOB_DIR}"
echo -n "send_production.sh :: >> runs      = "
for run in ${RUN_NUMBER_ARR[@]}; do
    echo -n ''$(printf "%02d" ${run})' '
done
echo ""
echo "send_production.sh :: "

for run in ${RUN_NUMBER_ARR[@]}; do

    STR_RUN=run$(printf "%02d" ${run})

    JOB_NAME="--job-name=b${MAG_FIELD}-j${JOB_DIR}-r$(printf "%02d" ${run})"
    sbatch ${JOB_NAME} ${SBATCH_SETUP} single_run.sh ${STR_RUN}
done
