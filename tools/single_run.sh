#!/bin/bash

if [[ -z ${OUTPUT_DIR} || -z ${JOB_DIR} || -z ${N_EVENTS_PER_RUN} ]]; then
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
    echo "single_run.sh :: reconstructing event ${event_id}"
    ./main ${MC_CSV} ${TRAJ_CSV} ${ITS_CSV} ${TPC_CSV} &> ${RECO_LOG}
done
cat event*_reco.log > ${STR_RUN}_reco.log
rm event*_reco.log

# merge csv files (note to self: can be improved...)
./merge_files.sh $(readlink -f .)
tmp_dir=$(readlink -f .)
tmp_dir=${tmp_dir##*/}
mv -v ${tmp_dir}_mc.csv ${STR_RUN}_mc.csv
mv -v ${tmp_dir}_traj.csv ${STR_RUN}_traj.csv
mv -v ${tmp_dir}_its.csv ${STR_RUN}_its.csv
mv -v ${tmp_dir}_tpc.csv ${STR_RUN}_tpc.csv
