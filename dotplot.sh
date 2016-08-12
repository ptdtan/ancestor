s#!/bin/bash
DELTA=$1
grep ">" $DELTA > chroms
sed -i 's/>//g' chroms
echo "#!/bin/bash" >> plot.sh
awk -F $" " '{print "mummerplot -r "$1" -q "$2" -p "$1"vs"$2" --png "${DELTA}"}' chroms > plot.sh
