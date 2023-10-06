#!/bin/bash

# default value
COLLISION_TYPE="pp"
# if first arg non-empty, then choose: hi (heavy-ion) or pp (proton-proton)
if [[ -n $1 ]]; then
    COLLISION_TYPE=${1}
fi

make main_fct
./main_fct config_${COLLISION_TYPE}.cmnd &> output_fct.log
