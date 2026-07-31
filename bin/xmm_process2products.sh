#!/bin/bash

filename=${1}
name=${2}
source=$(echo ${name// /_})

if [ -z "$1" ]
  then
    echo "No filename provided, exiting..."
    return 1
  fi

if [ -z "$2" ]
  then
    echo "No target name provided, exiting..."
    return 1
  fi

if [ -z "$3" ]
  then
    echo "No directory provided, using '.'"
    directory=.
else
  directory=${3}
fi


echo "filename='${filename}' source='${source}' directory='${directory}'"

if [ ! -d "$directory" ]; then
  echo "$directory does not exist, aborting."
  return 1
fi

bname=$(basename "$filename" | cut -d. -f1)


# 'heainit'/'sasinit' are aliases defined in ~/.bashrc; aliases are not
# expanded in non-interactive shells, so calling them directly here is a
# silent no-op ("command not found") and SAS/HEASoft never end up on PATH
# for this script or any of the stages it sources below. Source the init
# scripts directly instead.
: "${HEADAS:=/home/majestix/hdd/tools/heasoft-6.34/x86_64-pc-linux-gnu-libc2.39}"
: "${SAS_DIR:=/home/majestix/hdd/tools/sas/xmmsas_20230412_1735}"
: "${SAS_CCFPATH:=/home/majestix/hdd/tools/XMM_CCF}"
export HEADAS SAS_DIR SAS_CCFPATH

echo "Using HEADAS=${HEADAS}"
echo "Using SAS_DIR=${SAS_DIR}"
echo "Using SAS_CCFPATH=${SAS_CCFPATH}"

source "${HEADAS}/headas-init.sh"
source "${SAS_DIR}/setsas.sh"

# Must come after sourcing setsas.sh: setsas.sh caches the caller's PATH in
# SAS_PREV_PATH on its first run and restores PATH from that cache on every
# later invocation (e.g. if 'sasinit' already ran earlier in this shell
# session). Exporting these before setsas.sh would get silently overwritten,
# leaving xmm_process.py and friends off PATH ("command not found").
export PATH=$PATH:/home/majestix/hdd/tools/xmmpy/scripttools
export PATH=$PATH:/home/majestix/hdd/tools/xmmpy/bin
export PYTHONPATH=$PYTHONPATH:/home/majestix/hdd/tools

echo "Recognized parameters; filename='${filename}', target='${name}', directory='${directory}'"
ofn=${directory}/xmm_process_${bname}.sh
echo "Writing XMM-processing script to ${ofn}"
lfn=${ofn%.sh}.log
xmm_process.py $filename -target="${name}" --script=$ofn
echo "Running script written to '${ofn}' and log to '${lfn}'"
source ${ofn}  2>&1 | tee ${lfn}

ofn=${directory}/${bname}/xmm_source_regions_${source}.sh
lfn=${ofn%.sh}.log
echo "Writing source region script to ${ofn}"
xmm_source_regions.py ${directory}/${bname} "${name}" --script=${ofn}
source ${ofn}  2>&1 | tee ${lfn}

ofn=${directory}/${bname}/xmm_source_products_${source}.sh
lfn=${ofn%.sh}.log
echo "Writing source product script to ${ofn}"
xmm_source_products.py ${directory}/${bname} "${name}" --script=${ofn}
source ${ofn}  2>&1 | tee ${lfn}
