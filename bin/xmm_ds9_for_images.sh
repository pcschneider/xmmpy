#!/bin/bash

arg1=$1
arg2=$2
pythoncall=/home/majestix/hdd/python/bin/python3.7
xmmpy=/home/majestix/hdd/tools/xmmpy

$pythoncall ${xmmpy}/bin/xmm_ds9_for_images.py "$arg1" "$arg2"
