#!/bin/bash
#usage $0 includedir VARNAME

if [ "$#" != 2 ]
then
	echo "Error: invalid usage $0 $@" >&2
	exit 1
fi
CHARMINC="$1"
VAR="$2"

. $CHARMINC/conv-config.sh $CHARMINC

eval "printf \"%s\\n\" \"\${${VAR}}\""
