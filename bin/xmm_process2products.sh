#!/bin/bash

obsID=${1}
name=${2}
source=$(echo ${name// /_})

if [ -z $1 ] 
  then
    echo "No filename provided, exiting..."
    return 1
  fi

if [ -z $2 ] 
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


echo "obsID='${obsID}' source='${source}' directory='${directory}'"

if [ ! -d "$directory" ]; then
  echo "$directory does exist, aborting."
  exit 1
fi

filename=${1}
target=${2}
bname=$(basename "$filename" | cut -d. -f1)


export PATH=$PATH:/home/majestix/hdd/tools/xmmpy/scripttools
export PATH=$PATH:/home/majestix/hdd/tools/xmmpy/bin
export PYTHONPATH=$PYTHONPATH:/home/majestix/hdd/tools

heainit
sasinit

echo "Recognized parameters; filename='${filename}', target='${target}', 'directory='${directory}'"
ofn=${directory}/xmm_process_${bname}.sh
echo "Writing XMM-processing script to ${ofn}"

lfn=${ofn%.sh}.log
xmm_process.py $filename -target=${target} --script=$ofn
echo "Running script written to '${ofn}' and log to '${lfn}'"
source ${ofn}  2>&1 | tee ${lfn}

ofn=${directory}/${bname}/xmm_source_regions_${source}.sh
lfn=${ofn%.sh}.log
echo "Writing source region script to ${ofn}"
xmm_source_regions.py ${directory}/${bname} "${name}" --script=${ofn}
source ${ofn}  2>&1 | tee ${lfn}

ofn=${directory}/${bname}/xmm_source_products_${source}.sh
echo "Writing source product script to ${ofn}"
xmm_source_products.py ${directory}/${bname} "${name}" --script=${ofn}
source ${ofn}  2>&1 | tee ${lfn}
