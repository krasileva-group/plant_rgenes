#!/bin/bash
IN_DIR="$1"
TSTAMP=$(date +%Y%m%d-%H%M%S)
CMDFILE=addspeciesname_commands.${TSTAMP}.txt
SCRIPTS=/tgac/workarea/group-tg/projects/Disease_resistance/scripts
NUMCORES=3

fastafiles=$(find $IN_DIR -name '*fasta')

for FILE in $fastafiles

do

    echo "FILE: $FILE"
    basename=$(basename $FILE .fasta)
    echo "BASENAME: $basename"
    dirname=$(dirname $FILE)
    echo "DIR: $dirname"
    species=$(echo $basename | cut -f1 -d"_")
    echo "SPECIES: $species"

     CMDSTR="sed 's/>/>$species_/' $FILE > $dirname/$basename.addspeciesname.fasta"  

        echo "$CMDSTR"
	echo "$CMDSTR" >> $CMDFILE
	eval $CMDSTR

done
