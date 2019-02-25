#!/bin/bash
#
if [ $# -ne 2 ]; then
  echo "Error usage: $0 in.mrcs keyword.mrcs"
fi
#
source activate eman-env
imod $1 $2 --unstack

