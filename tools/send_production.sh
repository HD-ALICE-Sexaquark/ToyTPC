#!/bin/bash

   ##
  ##  Master Script to send batch jobs
 ##    to generate simulations
# # # # # # # # # # # # # # # # # # # # # # #

   ##
  ## Global Variables
 ## (Hard-coded parameters)
# # # # # # # # # # # # # # #

OUTPUT_DIR= # top output directory, will be assigned later based on the value of SIM_DIR, unless you modify it here
N_EVENTS_PER_RUN=10 # number of events per run (recommended: 100)

  ##
 ## User Interface Functions
# # # # # # # # # # # # # # #

function print_help() {
    echo "SCRIPT: send_production.sh"
    echo "=========================="
    echo "REQUIREMENTS:"
    echo "  * ROOT"
    echo "  * GEANT4"
    echo "USAGE:"
    echo "  ./send_production.sh --mode <mode> --run1 <run1> --run2 <run2> --serv <serv> --job <job_dir>"
    echo "  ./send_production.sh --mode <mode> --rn <rn> --serv <serv> --job <job_dir>"
    echo "  where:"
    echo "  <mode>    = it can be:"
    echo "              * 0 : send job to the htcondor farm"
    echo "              * 1 : interactive execution in a tmux session (use with caution!)"
    echo "  <job_dir> = subdirectory (it's a string!) within the output directory, for organizational purposes"
    echo "  <rn>      = specific run numbers, separated by comma"
    echo "              for example: --rn 0,1,2"
    echo "  <run1>    = run number of the first job (starting point of loop)"
    echo "  <run2>    = run number of the last job (end of loop)"
    echo "  <serv>    = (necessary and only valid when mode is 0)"
    echo "              select which machine to use: alice-serv<serv>"
    echo "              it can be: 10, 12 or 14"
    echo "EXAMPLE:"
    echo "  ./send_production.sh --mode 0 --job 01 --run1 0 --run2 99 --serv 10"
    echo "  ./send_production.sh --mode 1 --job 00 --run1 0 --run2 10"
}

function process_args() {
    arr=("$@")
    ic=0
    while [[ $ic -le $((${#arr[@]}-1)) ]]; do
        if [[ "${arr[$ic]}" == "--mode" ]]; then
            INT_MODE=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--job" ]]; then
            JOB_DIR=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--rn" ]]; then
            RUN_NUMBER_LIST=${arr[$((ic+1))]}
            RUN_NUMBER_ARR=(${RUN_NUMBER_LIST//,/ }) # convert comma-separated list into array of numbers
        elif [[ "${arr[$ic]}" == "--run1" ]]; then
            RUN1=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--run2" ]]; then
            RUN2=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--serv" ]]; then
            SERV=${arr[$((ic+1))]}
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

  ##
 ## TMUX management
# # # # # # # # # #

function prepare_and_exec_tmux() {

    # create execution file
    run_file="${STR_RUN}.sh"

    # (note to self: duplicated code, it's the same as in prepare_execution_script...)
    echo "#!/bin/bash"                                                         > ${run_file}
    echo ""                                                                   >> ${run_file}
    ### run injector
    echo "for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do"                >> ${run_file}
    echo '    event_id=$(printf "%03d" ${event})'                                   >> ${run_file}
    echo '    MC_CSV=event${event_id}_mc.csv'                                       >> ${run_file}
    echo '    MC_LOG=event${event_id}_mc.log'                                       >> ${run_file}
    echo '    root -l -b -q '\''GenBox.C("'\''${MC_CSV}'\''")'\'' &> ${MC_LOG}'     >> ${run_file}
    echo '    echo "'${run_file}' :: injecting event ${event_id}"'                  >> ${run_file}
    echo "done"                                                                     >> ${run_file}
    echo "cat event*_mc.log > ${STR_RUN}_mc.log"                                    >> ${run_file}
    echo "rm event*_mc.log"                                                         >> ${run_file}
    echo ""                                                                         >> ${run_file}
    ### do reconstruction
    echo "for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do"                >> ${run_file}
    echo '    event_id=$(printf "%03d" ${event})'                                   >> ${run_file}
    echo '    MC_CSV=event${event_id}_mc.csv'                                       >> ${run_file}
    echo '    TRAJ_CSV=event${event_id}_traj.csv'                                   >> ${run_file}
    echo '    ITS_CSV=event${event_id}_its.csv'                                     >> ${run_file}
    echo '    TPC_CSV=event${event_id}_tpc.csv'                                     >> ${run_file}
    echo '    RECO_LOG=event${event_id}_reco.log'                                   >> ${run_file}
    echo '    echo "'${run_file}' :: reconstructing event ${event_id}"'             >> ${run_file}
    echo '    ./main ${MC_CSV} ${TRAJ_CSV} ${ITS_CSV} ${TPC_CSV} &> ${RECO_LOG}'    >> ${run_file}
    echo 'done'                                                                     >> ${run_file}
    echo "cat event*_reco.log > ${STR_RUN}_reco.log"                                >> ${run_file}
    echo "rm event*_reco.log"                                                       >> ${run_file}
    echo ""                                                                         >> ${run_file}
    ### merge csv files
    echo './merge_files.sh $(readlink -f .)'                                        >> ${run_file}

    chmod a+x ${run_file}

    # create tmux session
    THIS_SESSION="${JOB_DIR}_${STR_RUN}"
    tmux new -d -s ${THIS_SESSION}

    # prepare environment
    tmux send -t ${THIS_SESSION} "source set_env.sh" ENTER
    tmux send -t ${THIS_SESSION} "./${run_file} 2> stderr.log 1> stdout.log" ENTER
    tmux send -t ${THIS_SESSION} "exit" ENTER
}

  ##
 ## Jobs management
# # # # # # # # # #

function prepare_execution_script() {

    # create execution file
    run_file="${STR_RUN}.sh"

    echo "#!/bin/bash"                                                               > ${run_file}
    echo ""                                                                         >> ${run_file}
    ### set environment
    echo "source set_env.sh"                                                        >> ${run_file}
    echo ""                                                                         >> ${run_file}
    ### run injector
    echo "for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do"                >> ${run_file}
    echo '    event_id=$(printf "%03d" ${event})'                                   >> ${run_file}
    echo '    MC_CSV=event${event_id}_mc.csv'                                       >> ${run_file}
    echo '    MC_LOG=event${event_id}_mc.log'                                       >> ${run_file}
    echo '    root -l -b -q '\''GenBox.C("'\''${MC_CSV}'\''")'\'' &> ${MC_LOG}'     >> ${run_file}
    echo '    echo "'${run_file}' :: injecting event ${event_id}"'                  >> ${run_file}
    echo "done"                                                                     >> ${run_file}
    echo "cat event*_mc.log > ${STR_RUN}_mc.log"                                    >> ${run_file}
    echo "rm event*_mc.log"                                                         >> ${run_file}
    echo ""                                                                         >> ${run_file}
    ### do reconstruction
    echo "for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do"                >> ${run_file}
    echo '    event_id=$(printf "%03d" ${event})'                                   >> ${run_file}
    echo '    MC_CSV=event${event_id}_mc.csv'                                       >> ${run_file}
    echo '    TRAJ_CSV=event${event_id}_traj.csv'                                   >> ${run_file}
    echo '    ITS_CSV=event${event_id}_its.csv'                                     >> ${run_file}
    echo '    TPC_CSV=event${event_id}_tpc.csv'                                     >> ${run_file}
    echo '    RECO_LOG=event${event_id}_reco.log'                                   >> ${run_file}
    echo '    echo "'${run_file}' :: reconstructing event ${event_id}"'             >> ${run_file}
    echo '    ./main ${MC_CSV} ${TRAJ_CSV} ${ITS_CSV} ${TPC_CSV} &> ${RECO_LOG}'    >> ${run_file}
    echo 'done'                                                                     >> ${run_file}
    echo "cat event*_reco.log > ${STR_RUN}_reco.log"                                >> ${run_file}
    echo "rm event*_reco.log"                                                       >> ${run_file}
    echo ""                                                                         >> ${run_file}
    ### merge csv files (note to self: can be improved...)
    echo './merge_files.sh $(readlink -f .)'                                        >> ${run_file}
    echo 'tmp_dir=$(readlink -f .)'                                                 >> ${run_file}
    echo 'tmp_dir=${tmp_dir##*/}'                                                   >> ${run_file}
    echo 'mv -v ${tmp_dir}_mc.csv '${STR_RUN}'_mc.csv'                              >> ${run_file}
    echo 'mv -v ${tmp_dir}_traj.csv '${STR_RUN}'_traj.csv'                          >> ${run_file}
    echo 'mv -v ${tmp_dir}_its.csv '${STR_RUN}'_its.csv'                            >> ${run_file}
    echo 'mv -v ${tmp_dir}_tpc.csv '${STR_RUN}'_tpc.csv'                            >> ${run_file}

    chmod a+x ${run_file}
}

function create_and_send_job() {

    # create job file
    job_file="${STR_RUN}.job"

    echo "executable            = ${STR_RUN}.sh"                                     > ${job_file}
    echo "output                = stdout.log"                                       >> ${job_file}
    echo "error                 = stderr.log"                                       >> ${job_file}
    echo "log                   = log.condor"                                       >> ${job_file}
    echo "should_transfer_files = YES"                                              >> ${job_file}
    echo "transfer_input_files  = set_env.sh,GenBox.C,main,merge_files.sh"          >> ${job_file}
    echo 'requirements          = (Machine == "alice-serv'${SERV}'")'               >> ${job_file}
    echo "queue"                                                                    >> ${job_file}

    # print info about output directory and copied files (debug purposes)
    echo "send_production.sh :: displaying content of ${PWD}:"
    for file in *.*; do echo "send_production.sh :: -- ${file}"; done

    # send job
    echo "send_production.sh :: sending ${job_file}..."
    condor_submit ${job_file}

    echo "send_production.sh ::"
}

  ##
 ## Process Input
# # # # # # # # # #

argArray=("$@")
process_args "${argArray[@]}"

  ##
 ## Check for command-line errors
# # # # # # # # # # # # # # # # # #

if [[ -z ${INT_MODE} ]]; then
    echo "ERROR: --mode option empty, use --help to see available commands"
    exit 1
fi

if [[ ${INT_MODE} -eq 0 ]] && [[ -z ${SERV} ]]; then
    echo "ERROR: --serv option empty, use --help to see available commands"
    exit 1
fi

if [[ -z ${RUN1} ]] && [[ -z ${RUN2} ]] && [[ -z ${RUN_NUMBER_LIST} ]]; then
    echo "ERROR: run number options empty, use --help to see available commands"
    exit 1
fi

if [[ -n ${RUN1} || -n ${RUN2} ]] && [[ -n ${RUN_NUMBER_LIST} ]]; then
    echo "ERROR: you must choose between --rn OR --run1 and --run2, not combine them!"
    exit 1
fi

  ##
 ## Set SIM_DIR and OUTPUT_DIR
# # # # # # # # # # # # # # # # # #

if [[ "${PWD##*/}" == "tools" ]]; then
    SIM_DIR=${PWD}/..
else
    SIM_DIR=${PWD}
fi

if [[ -z ${OUTPUT_DIR} ]]; then
    echo 'WARNING: OUTPUT_DIR empty, setting default: ${SIM_DIR}/output'
    OUTPUT_DIR=${SIM_DIR}/output
fi

  ##
 ## Check for environment
# # # # # # # # # # # # # #

if [[ ! $(geant4-config) ]]; then
    echo "ERROR: make sure to set GEANT4"
    exit 1
fi

if [[ -z ${ROOTSYS} ]]; then
    echo "ERROR: make sure to set ROOT"
    exit 1
fi

  ##
 ## Check for binaries
# # # # # # # # # # # #

if [[ ! -e ${SIM_DIR}/reconstruction/G4Setup_build/main ]]; then
    echo "ERROR: reconstruction/G4Setup_build/main not found. Maybe forgot to compile?"
    exit 1
fi

  ##
 ## Main
# # # # #

# if output directory doesn't exist, create it
mkdir -p ${OUTPUT_DIR}

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
echo "send_production.sh :: >> OUTPUT_DIR       = ${OUTPUT_DIR}"
echo "send_production.sh :: >> N_EVENTS_PER_RUN = ${N_EVENTS_PER_RUN}"
echo "send_production.sh ::"
echo "send_production.sh :: chosen options:"
echo "send_production.sh :: >> sim_dir  = ${SIM_DIR}"
echo "send_production.sh :: >> mode     = ${INT_MODE}"
echo "send_production.sh :: >> job_dir  = ${JOB_DIR}"
echo -n "send_production.sh :: >> runs     = "
for run in ${RUN_NUMBER_ARR[@]}; do
    echo -n ''$(printf "%02d" ${run})' '
done
echo ""
echo "send_production.sh :: "

# create set_env.sh in case we're sending jobs to the HTCondor farm
if [[ ! -e ${SIM_DIR}/set_env.sh ]]; then

    ROOT_ENV_SCRIPT=${ROOTSYS}/bin/thisroot.sh
    G4_ENV_SCRIPT=$(geant4-config --prefix)/bin/geant4.sh

    echo "#!/bin/bash"                > ${SIM_DIR}/set_env.sh
    echo ""                          >> ${SIM_DIR}/set_env.sh
    echo "source ${ROOT_ENV_SCRIPT}" >> ${SIM_DIR}/set_env.sh
    echo "source ${G4_ENV_SCRIPT}"   >> ${SIM_DIR}/set_env.sh
    chmod a+x ${SIM_DIR}/set_env.sh
fi

# start loop over runs
for run in ${RUN_NUMBER_ARR[@]}; do

    # define run number
    STR_RUN=run$(printf "%02d" ${run})

    # define the specific output dir for this run
    RUN_DIR=${OUTPUT_DIR}/${JOB_DIR}/${STR_RUN}
    mkdir -p ${RUN_DIR}

    # enter output dir
    cd ${RUN_DIR}

    # bring necessary files
    cp ${SIM_DIR}/set_env.sh .
    cp ${SIM_DIR}/injector/GenBox.C .
    cp ${SIM_DIR}/reconstruction/G4Setup_build/main .
    cp ${SIM_DIR}/tools/merge_files.sh .

    if [[ ${INT_MODE} -eq 1 ]]; then
        prepare_and_exec_tmux
    else
        prepare_execution_script
        create_and_send_job
    fi

    # go back to sim dir
    cd ${SIM_DIR}
done

echo "send_production.sh :: all sent."
