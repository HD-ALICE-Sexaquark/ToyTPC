#!/bin/bash

   ##
  ##  Master Script to send batch jobs
 ##    to generate simulations
# # # # # # # # # # # # # # # # # # # # # # #

   ##
  ## Global Variables
 ## (Hard-coded parameters)
# # # # # # # # # # # # # # #

N_EVENTS_PER_RUN=10 # number of events per run (recommended: 1000)

  ##
 ## User Interface Functions
# # # # # # # # # # # # # # #

function print_help() {
    echo "SCRIPT: send_production.sh"
    echo "=========================="
    echo "REQUIREMENTS:"
    echo "  * ROOT"
    echo "  * GEANT4"
    echo "  * PYTHIA8"
    echo "USAGE:"
    echo "  ./send_production.sh --mode <mode> --sim_dir <sim_dir> --out_dir <out_dir> --run1 <run1> --run2 <run2> --serv <serv>"
    echo "  ./send_production.sh --mode <mode> --sim_dir <sim_dir> --out_dir <out_dir> --rn <rn> --serv <serv>"
    echo "  where:"
    echo "  <mode>        = it can be:"
    echo "                  * 0 : send job to the HTCondor farm"
    echo "                  * 1 : run processes in current terminal session"
    echo "  <sim_dir>     = main directory"
    echo '                  (default: if ${PWD} == tools/, then sim_dir = "tools/.."'
    echo '                            if not, then sim_dir = ${PWD})'
    echo "  <out_dir>     = output directory"
    echo '                  (default: <sim_dir>/output")'
    echo "  <rn>          = specific run numbers, separated by comma"
    echo "                  for example: --rn 0,1,2"
    echo "  <run1>        = run number of the first job (starting point of loop)"
    echo "  <run2>        = run number of the last job (end of loop)"
    echo "  <serv>        = (only valid and mandatory when mode == 0) select which machine to use: alice-serv<serv>"
    echo "                  it can be: 10, 12, 13 or 14"
    echo "                  IMPORTANT: make sure to send a max of 8 jobs per serv"
    echo "EXAMPLES:"
    echo "  ./send_production.sh --mode 1 --rn 14"
    echo "  ./send_production.sh --mode 0 --run1 1 --run2 8 --serv 14"
}

function process_args() {
    arr=("$@")
    ic=0
    while [[ $ic -le $((${#arr[@]}-1)) ]]; do
        if [[ "${arr[$ic]}" == "--mode" ]]; then
            INT_MODE=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--rn" ]]; then
            RUN_NUMBER_LIST=${arr[$((ic+1))]}
            RUN_NUMBER_ARR=(${RUN_NUMBER_LIST//,/ }) # convert comma-separated list into array of numbers
        elif [[ "${arr[$ic]}" == "--run1" ]]; then
            RUN1=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--run2" ]]; then
            RUN2=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--serv" ]]; then
            SERV=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--sim_dir" ]]; then
            SIM_DIR=${arr[$((ic+1))]}
        elif [[ "${arr[$ic]}" == "--out_dir" ]]; then
            OUTPUT_DIR=${arr[$((ic+1))]}
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
 ## Event Generation
# # # # # # # # # # #

function generate_events() {
    # generate pp collisions
    PYTHIA_CFG=config_pp.cmnd
    PYTHIA_LOG=${str_run}_mc.log
    # the loop is done internally
    echo "send_production.sh :: >> generating events with PYTHIA8"
    ./main_GenCollision --n ${N_EVENTS_PER_RUN} --config ${PYTHIA_CFG} --output . &> ${PYTHIA_LOG} &
    echo "send_production.sh :: " # empty line
}

  ##
 ## Reconstruction
# # # # # # # # # #

function do_reconstruction() {

    N_NEEDED_MC_FILES=${N_EVENTS_PER_RUN}

    # wait until all files are ready
    while true; do
        # count mc files
        N_MC_FILES=$(ls -1 event*_mc.csv 2> /dev/null | wc -l)
        if [[ ${N_MC_FILES} -eq ${N_NEEDED_MC_FILES} ]]; then
            break
        fi
        # wait 1 second until next iteration
        sleep 1
    done

    # start reconstruction process
    for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do
        # set filenames
        str_event="$(printf "%03d" ${event})"
        # - input
        MC_CSV=event${str_event}_mc.csv
        # - output
        TRAJ_CSV=event${str_event}_traj.csv
        ITS_CSV=event${str_event}_its.csv
        TPC_CSV=event${str_event}_tpc.csv
        # - log
        RECO_LOG=event${str_event}_reco.log

        # run sequentially, geant4 takes care of the parallelization
        echo "send_production.sh :: >> reconstructing event ${str_event}"
        ./main ${MC_CSV} ${TRAJ_CSV} ${ITS_CSV} ${TPC_CSV} &> ${RECO_LOG}
    done

    # merge all log files into one
    RUN_RECO_LOG=${str_run}_reco.log
    cat event*_reco.log > ${RUN_RECO_LOG}
    rm event*_reco.log

    echo "send_production.sh :: " # empty line
}

  ##
 ## Jobs management
# # # # # # # # # #

function prepare_execution_script() {

    # copy env scripts
    ## alienv printenv AliDPG/latest > set_env.sh
    cp ${SIM_DIR}/set_env.sh .

    # create run/command file (within OUTDIR)
    run_file="${str_run}.sh"
    echo "#!/bin/bash"                                                                 > ${run_file}
    echo ""                                                                           >> ${run_file}
    ### environment
    echo "source set_env.sh"                                                          >> ${run_file}
    echo ""                                                                           >> ${run_file}
    ### inject bkg
    echo "for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do"                  >> ${run_file}
    echo '    str_event=$(printf "%03d" ${event})'                                    >> ${run_file}
    echo '    BKG_CSV=event${str_event}_bkg.csv'                                      >> ${run_file}
    echo '    BKG_LOG=event${str_event}_bkg.log'                                      >> ${run_file}
    echo '    root -l -b -q '\''GenBox.C('${BKG_OPT}', "'\''${BKG_CSV}'\''")'\'' &> ${BKG_LOG}' >> ${run_file}
    echo '    echo "'${run_file}' :: sending bkg ${str_event}"'                       >> ${run_file}
    echo "done"                                                                       >> ${run_file}
    echo "cat event*_bkg.log > ${str_run}_bkg.log"                                    >> ${run_file}
    echo "rm event*_bkg.log"                                                          >> ${run_file}
    ### inject signal
    if [[ ${ONLY_BKG} -eq 0 ]]; then
        echo "for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do"                          >> ${run_file}
        echo '    str_event=$(printf "%03d" ${event})'                                            >> ${run_file}
        echo '    SIGNAL_CSV=event${str_event}_sig.csv'                                           >> ${run_file}
        echo '    SIGNAL_LOG=event${str_event}_sig.log'                                           >> ${run_file}
        echo '    root -l -b -q '\''GenSexaquarkReaction.C("'\''${SIGNAL_CSV}'\''")'\'' &> ${SIGNAL_LOG}' >> ${run_file}
        echo '    echo "'${run_file}' :: sending signal ${str_event}"'                            >> ${run_file}
        echo 'done'                                                                               >> ${run_file}
        echo "cat event*_sig.log > ${str_run}_sig.log"                                            >> ${run_file}
        echo "rm event*_sig.log"                                                                  >> ${run_file}
    fi
    ### do reconstruction
    echo "for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do"                              >> ${run_file}
    echo '    str_event=$(printf "%03d" ${event})'                                                >> ${run_file}
    if [[ ${ONLY_BKG} -eq 0 ]]; then
        echo '    SIGNAL_CSV=event${str_event}_sig.csv'                                           >> ${run_file}
    else
        echo '    SIGNAL_CSV=0'                                                                   >> ${run_file}
    fi
    echo '    BKG_CSV=event${str_event}_bkg.csv'                                                  >> ${run_file}
    echo '    RECO_CSV=event${str_event}_reco.csv'                                                >> ${run_file}
    echo '    RECO_LOG=event${str_event}_reco.log'                                                >> ${run_file}
    echo '    echo "'${run_file}' :: running reco ${str_event}"'                                  >> ${run_file}
    echo '    ./exampleB2a '${BKG_OPT}' ${SIGNAL_CSV} ${BKG_CSV} ${RECO_CSV} 1 &> ${RECO_LOG}'    >> ${run_file}
    echo 'done'                                                                                   >> ${run_file}

    chmod a+x ${run_file}
}

function create_and_send_job() {

    # create job file (within OUTDIR)
    job_file="${str_run}.job"
    echo "executable            = ${str_run}.sh"                                              > ${job_file}
    echo "output                = stdout.log"                                                >> ${job_file}
    echo "error                 = stderr.log"                                                >> ${job_file}
    echo "log                   = log.condor"                                                >> ${job_file}
    echo "should_transfer_files = YES"                                                       >> ${job_file}
    if [[ ${ONLY_BKG} -eq 0 ]]; then
        echo "transfer_input_files  = set_env.sh,GenBox.C,GenSexaquarkReaction.C,exampleB2a" >> ${job_file}
    else
        echo "transfer_input_files  = set_env.sh,GenBox.C,exampleB2a"                        >> ${job_file}
    fi
    echo 'requirements          = (Machine == "alice-serv'${SERV}'")'                        >> ${job_file}
    echo "queue"                                                                             >> ${job_file}

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

if [[ -z ${SIM_DIR} ]]; then
    echo "WARNING: --sim_dir option empty, setting default"
    if [[ "${PWD##*/}" == "tools" ]]; then SIM_DIR="${PWD}/..";
    else SIM_DIR=${PWD};
    fi
fi

if [[ -z ${OUTPUT_DIR} ]]; then
    echo "WARNING: --out_dir option empty, setting default"
    OUTPUT_DIR=${SIM_DIR}/output  # output directory (recommended: a location with storage)
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

# if [[ -z ${PYTHIA8} ]]; then
    # echo "ERROR: make sure to set PYTHIA8"
    # exit 1
# fi

  ##
 ## Check for binaries
# # # # # # # # # # # #

if [[ ! -e ${SIM_DIR}/generator/main_GenCollision ]]; then
    echo "ERROR: generator/main_GenCollision not found. Maybe forgot to compile?"
    exit 1
fi

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
echo "send_production.sh :: >> N_EVENTS_PER_RUN = ${N_EVENTS_PER_RUN}"
echo "send_production.sh ::"
echo "send_production.sh :: chosen options:"
echo "send_production.sh :: >> sim_dir  = ${SIM_DIR}"
echo "send_production.sh :: >> out_dir  = ${OUTPUT_DIR}"
echo "send_production.sh :: >> mode     = ${INT_MODE}"
if [[ ${INT_MODE} -eq 0 ]]; then
    echo "send_production.sh :: >> serv     = ${SERV}"
fi
echo -n "send_production.sh :: >> runs     = "
for run in ${RUN_NUMBER_ARR[@]}; do
    echo -n ''$(printf "%03d" ${run})' '
done
echo ""
echo "send_production.sh :: "

# create set_env.sh in case we're sending jobs to the HTCondor farm
if [[ ${INT_MODE} -eq 0 ]] && [[ ! -e ${SIM_DIR}/set_env.sh ]]; then

    ROOT_ENV_SCRIPT=${ROOTSYS}/bin/thisroot.sh
    G4_ENV_SCRIPT=$(geant4-config)/bin/geant4.sh

    echo "#!/bin/bash"                > ${SIM_DIR}/set_env.sh
    echo ""                          >> ${SIM_DIR}/set_env.sh
    echo "source ${ROOT_ENV_SCRIPT}" >> ${SIM_DIR}/set_env.sh
    echo "source ${G4_ENV_SCRIPT}"   >> ${SIM_DIR}/set_env.sh
    chmod a+x ${SIM_DIR}/set_env.sh
fi

# start loop
for run in ${RUN_NUMBER_ARR[@]}; do
    # define run number
    str_run=run$(printf "%03d" ${run})

    # define an create output dir
    # and enter it (important to do so, because condor requires input files to be in the same dir)
    RUN_DIR=${OUTPUT_DIR}/${str_run}
    mkdir -p ${RUN_DIR}

    # enter run dir
    cd ${RUN_DIR}

    # bring necessary files
    cp ${SIM_DIR}/generator/config_pp.cmnd .
    cp ${SIM_DIR}/generator/main_GenCollision .
    cp ${SIM_DIR}/reconstruction/G4Setup_build/main .
    cp ${SIM_DIR}/tools/merge_files.sh .

    if [[ ${INT_MODE} -eq 1 ]]; then
        echo "send_production.sh :: starting run #$(printf "%03d" ${run})"
        echo "send_production.sh :: " # empty line
        generate_events
        do_reconstruction
        ./merge_files.sh ${RUN_DIR}
    else
        prepare_execution_script
        create_and_send_job
    fi

    # remove copied files
    rm config_pp.cmnd
    rm main_GenCollision
    rm main
    rm merge_files.sh

    # go back to sim dir
    cd ${SIM_DIR}
done

echo "send_production.sh :: " # empty line
echo "send_production.sh :: all done."
