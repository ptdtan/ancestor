#!/bin/bash
#./findmax_scaffold <scaffold_links>
LINKS=$1
awk -v max=0 '{if($3>max) {max=$2; want=$1}}END{print max"\t" want}' $LINKS
