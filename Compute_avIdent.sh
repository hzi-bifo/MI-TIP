#!/bin/bash
for f in $(cat $1)
do
	av_ident=$(trimal -sident -in $f |grep '\#Mean Percentage of identity\:'| awk -F ':\\s+' '{print $2}')
	echo -e $f"\t"$av_ident
done;
