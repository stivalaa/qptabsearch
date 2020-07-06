#!/bin/sh
#
# File:    build_skolnick_db.sh
# Author:  Alex Stivala
# Created: September 2008
#
# build_skolnick_db.sh - build tableaux database for Skolnick data set
#
# Usage: build_skolnick_db.sh outdir 
#
#   outdir is name of diretory which is created, and each tableau 
#   in ASCII format for use with tsrchd_sparse etc. is
#   created as a separate file in that directory, in format for input
#   for use with qptabmatch_allpairs.py for example
#
# builds database of tableaux, using pytableaycreate.py,
# for the Skolnick data set (see Lancia et al.
# RECOMB 2001) using listing of the skolnick.tar.gz file supplied as
# supplementary material for Pelta et al 2008 BMC Bioinformatics 9:161
# http://modo.ugr.es/files/jrgonzalez/MSVNS4MaxCMO/datasets/skolnick.tar.gz,
# embedded in here (note edit: 1AW2_A rather than 1AW2)
#
#

# root of divided PDB hierarchy
PDBROOT=/local/charikar/pdb/pdb


# list of PDB ids in Skolnick data set
SKOLNICK_LIST="1BTM 1BYO_B 1TMH 1AW2_A 1IER 3CHY 1NIN 2B3I 1RCD 1TRI 3YPI 1AMK 8TIM 1DPS 1NTR 1QMP_C 1FHA 4TMY_B 1RN1_C 1PLA 1HTI 1YDV 1BCF 1RN1_A 2PLT 1TRE 1QMP_A 1QMP_D 2PCY 1QMP_B 1BYO_A 1DBW 1B9B 4TMY_A 1NAT 1B71 1RN1_B 1KDI 1BAW 1B00"


if [ $# -ne 1 ]; then
    echo "Usage: $0  outdir" 2>&1
    exit 1
fi
outdir=$1

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

# pytableaucreate.py options
tabopts="-35 -f -t dssp -p none"

for i in $SKOLNICK_LIST
do
    pdb=`echo $i | tr A-Z a-z`
    if [ `expr index $pdb _` -ne 0 ]; then
        # get chainid from e.g. 1BYO_B
        chainid=`expr substr $pdb 6 1`
        chainopt="-c $chainid"
        pdbid=`expr substr $pdb 1 4`_${chainid}
    else
        chainopt=""
        pdbid=`expr substr $pdb 1 4`
    fi
    pdb=`expr substr $pdb 1 4`
    div=`expr substr $pdb 2 2`
    pdbfile=/local/charikar/pdb/pdb/${div}/pdb${pdb}.ent.gz
    cat /dev/null > ${outdir}/${pdbid}.tableaudistmatrix
    pytableaucreate.py ${tabopts} ${chainopt} ${pdbfile} >>${outdir}/${pdbid}.tableaudistmatrix
    # append distance matrix, removing identifier on first line
    pytableaucreate.py -d ${tabopts} ${chainopt} ${pdbfile} | awk 'NR > 1' >>${outdir}/${pdbid}.tableaudistmatrix 
done

