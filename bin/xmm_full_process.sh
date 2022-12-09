#!/bin/bash

obsID=${1}
name=${2}
source=$(echo ${name// /_})

if [ -z $1 ] 
  then
    echo "No obsID provided, exiting..."
    exit 1
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


export PATH=$PATH:/home/majestix/hdd/tools/xmmpy/scripttools
export PYTHONPATH=$PYTHONPATH:/home/majestix/hdd/tools


ofn=${directory}/xmm_retrieve_${obsID}.sh
echo "Writing retrieve script to ${ofn}"
lfn=${ofn%.sh}.log
xmm_retrieve.py $obsID ${directory} --script=$ofn
echo "Running script written to '${ofn}' and log to '${lfn}'"
source ${ofn}  2>&1 | tee ${lfn}

ofn=${directory}/${obsID}/xmm_source_regions_${source}.sh
lfn=${ofn%.sh}.log
echo "Writing source region script to ${ofn}"
xmm_source_regions.py ${directory}/${obsID} "${name}" --script=${ofn}
source ${ofn}  2>&1 | tee ${lfn}

ofn=${directory}/${obsID}/xmm_source_products_${source}.sh
echo "Writing source product script to ${ofn}"
xmm_source_products.py ${directory}/${obsID} "${name}" --script=${ofn}
source ${ofn}  2>&1 | tee ${lfn}
