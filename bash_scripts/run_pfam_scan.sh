#!/bin/bash

# Usage: bash run_pfam_scan.sh dir

# Dependencies
# 1. HMMER software (http://hmmer.janelia.org/)
# including pfam_scan.pl (part of HMMER) Move in same directory as this script or set path at command string.
# 2. Pfam database (http://pfam.xfam.org/)
# 3. File names should be consistent with Phytozome and include Species_*_protein.fa

# specify directory in which you want to search for all *protein.fa files and annotate them with pfamscan
IN_DIR="$1" 

# set location of Pfam db
Pfam=

# set number of jobs to run in parallel
NUMCORES=1

# record current date and time and create a commands file.
TSTAMP=$(date +%Y%m%d-%H%M%S)
CMDFILE=commands.${TSTAMP}.txt

protfiles=$(find $IN_DIR -name '*protein.fa')
for FILE in $protfiles
do
    echo "FILE: $FILE"
    basename=$(basename $FILE _protein.fa)
    echo "BASENAME: $basename"
    dirname=$(dirname $FILE)
    echo "DIR: $dirname"
    species=$(echo $basename | cut -f1 -d"_")
    echo "SPECIES: $species"
    if [ ! -e "$dirname/pfam/" ]; then
	mkdir $dirname/pfam
        echo "mkdir $dirname/pfam"
    fi
    CMDSTR+="'time perl pfam_scan.pl -e_seq 1 -e_dom 1 -pfamB -as -outfile $dirname/pfam/${basename}_pfamscan-$(date +%m-%d-%Y).out -cpu 8 -fasta $FILE -dir $Pfam"
    CMDSTR+="'"$'\n'

done
echo "$CMDSTR" > $CMDFILE
xargs="cat $CMDFILE | xargs -P ${NUMCORES} -n 1 sh -c"
echo $xargs
eval $xargs
