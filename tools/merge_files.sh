#!/bin/bash

   ###
  ###  Merge output CSV files
 ##  by adding the event number as an additional column
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

RUN_DIR=$1
RUN_STR=${RUN_DIR##*/}

FILE_TYPES_ARR=("mc" "traj" "its" "tpc")

for file_type in ${FILE_TYPES_ARR[@]}; do

    n_files=$(ls -1 ${RUN_DIR}/event*_${file_type}.csv 2> /dev/null | wc -l)

    for ((event_number = 0; event_number < n_files; event_number++)); do

        event_str=$(printf "%03d" ${event_number})
        single_file=${RUN_DIR}/event${event_str}_${file_type}.csv
        merged_file=${RUN_DIR}/${RUN_STR}_${file_type}.csv

        if [[ ${event_number} -eq 0 ]]; then
            awk -F',' 'BEGIN {OFS=FS} {NF++; for (i=NF; i>1; i--) $i = $(i-1); $1 = '${event_number}'} 1' ${single_file} > ${merged_file}
        else
            awk -F',' 'BEGIN {OFS=FS} {NF++; for (i=NF; i>1; i--) $i = $(i-1); $1 = '${event_number}'} 1' ${single_file} >> ${merged_file}
        fi

        rm event${event_str}_${file_type}.csv

    done

done
