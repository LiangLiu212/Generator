#!/bin/bash

current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export INCLXX_DIR=$(dirname "$current_dir")
export PATH=${current_dir}:${PATH}
export LD_LIBRARY_PATH=${INCLXX_DIR}/lib:${LD_LIBRARY_PATH}
export INCLXX_DATA_DIR=${INCLXX_DIR}/share
