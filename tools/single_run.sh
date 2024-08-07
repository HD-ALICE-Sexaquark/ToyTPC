#!/bin/bash

if [[ -z ${SIM_DIR} || -z ${OUTPUT_DIR} || -z ${MAG_FIELD} || -z ${JOB_DIR} || -z ${N_EVENTS_PER_RUN} || -z ${GEANT4_NUM_THREADS} ]]; then
    echo "single_run.sh :: ERROR :: options must be defined!"
    exit 1
fi

if [[ $# -ne 1 ]]; then
    echo "single_run.sh :: ERROR :: provide a run name!"
    exit 1
fi

STR_RUN=$1

# define paths
RUN_DIR=${OUTPUT_DIR}/${JOB_DIR}/${STR_RUN}
mkdir -p ${RUN_DIR}

# bring necessary files
cp ${SIM_DIR}/injector/GenBox.C ${RUN_DIR}
cp ${SIM_DIR}/reconstruction/G4Setup_build/main ${RUN_DIR}
cp ${SIM_DIR}/tools/merge_files.sh ${RUN_DIR}

cd ${RUN_DIR}

# print run info
RUN_LOG=${STR_RUN}_cfg.log
echo "${STR_RUN}" > ${RUN_LOG}
echo ">> SIM_DIR               = ${SIM_DIR}" >> ${RUN_LOG}
echo ">> OUTPUT_DIR            = ${OUTPUT_DIR}" >> ${RUN_LOG}
echo ">> N_EVENTS_PER_RUN      = ${N_EVENTS_PER_RUN}" >> ${RUN_LOG}
echo ">> GEANT4_NUM_THREADS    = ${GEANT4_NUM_THREADS}" >> ${RUN_LOG}
echo ">> MAG_FIELD             = ${MAG_FIELD}" >> ${RUN_LOG}
echo ">> TRIGGER_CONDITION_OFF = ${TRIGGER_CONDITION_OFF}" >> ${RUN_LOG}

# run injector
for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do
    event_id=$(printf "%03d" ${event})
    MC_CSV=event${event_id}_mc.csv
    MC_LOG=event${event_id}_mc.log
    root -l -b -q 'GenBox.C("'${MC_CSV}'")' &> ${MC_LOG}
    echo "single_run.sh :: injecting event ${event_id}"
done
cat event*_mc.log > ${STR_RUN}_mc.log
rm event*_mc.log

# do reconstruction
for ((event=0; event < ${N_EVENTS_PER_RUN}; event++)); do
    event_id=$(printf "%03d" ${event})
    MC_CSV=event${event_id}_mc.csv
    TRAJ_CSV=event${event_id}_traj.csv
    ITS_CSV=event${event_id}_its.csv
    TPC_CSV=event${event_id}_tpc.csv
    RECO_LOG=event${event_id}_reco.log
    TRIGGER_CONDITION=1
    if [[ ${TRIGGER_CONDITION_OFF} -eq 1 ]]; then
        TRIGGER_CONDITION=0
    fi
    echo "single_run.sh :: reconstructing event ${event_id}"
    ./main ${MC_CSV} ${TRAJ_CSV} ${ITS_CSV} ${TPC_CSV} ${TRIGGER_CONDITION} ${MAG_FIELD} &> ${RECO_LOG}
done
cat event*_reco.log > ${STR_RUN}_reco.log
rm event*_reco.log

# merge csv files
./merge_files.sh
echo "single_run.sh :: finished merging files"
