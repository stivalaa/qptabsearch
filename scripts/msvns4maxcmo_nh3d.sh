#!/bin/sh
###############################################################################
#
# msvns4maxcmo_nh3d.sh - run MSVNS4MaxCMO on Nh3D data set
#
# File:    msvns4maxcmo_nh3d.sh
# Author:  Alex Stivala
# Created: September 2008
#
# Run MSVNS4MaxCMO 
# the Nh3D data set (Thiruv et al 2005 BMC Struct. Biol. 5:12)
# with the 73 queries defined in Pelta et al 2008 BMC Bioinformatics 9:161
#
# Usage:
#     msvns4maxcmo_nh3d.sh indir outdir
#
#     indir is directory containing the contact matrices built with
#     build_nh3d_cm.sh
#
#     outdir is diretory to place corresponding output from MSVNS4MaxCMO
#     created if it does not exist
#     WARNNG: .out files in outdir overwritten if they exist
#
# Environment variables:
#
#   PATH must contain the location of MSVNS4MaxCMO
#
# $Id: msvns4maxcmo_nh3d.sh 1955 2008-10-06 23:48:52Z astivala $
# 
###############################################################################


# List of query CATH identifiers, from the Additional File 1 spreadsheet
# for Pelta et al 2008
QUERY_LIST="1.10.1040 1.10.1320 1.10.533 1.10.645 1.20.1280 1.20.210 1.20.5 1.20.840 2.10.25 2.10.260 2.10.270 2.10.90 2.170.16 2.170.230 2.170.290 2.170.40 2.30.110 2.30.18 2.30.230 2.30.29 2.30.40 2.40.155 2.40.160 2.40.180 2.40.340 2.40.50 2.60.130 2.60.260 2.60.420 2.60.90 2.70.100 2.70.180 2.70.220 2.70.98 3.10.105 3.10.170 3.10.270 3.10.330 3.10.400 3.20.120 3.20.140 3.20.19 3.20.70 3.20.90 3.30.1530 3.30.1690 3.30.240 3.30.559 3.30.560 3.30.60 3.30.990 3.40.1210 3.40.1380 3.40.225 3.40.720 3.60.100 3.60.120 3.60.20 3.60.40 3.60.90 3.90.1280 3.90.1300 3.90.1350 3.90.1580 3.90.510 3.90.850 4.10.1080 4.10.1090 4.10.220 4.10.260 4.10.480 4.10.540 4.10.790"

if [ $# -ne 2 ]; then
    echo "Usage: $0 indir outdir" 2>&1
    exit 1
fi

indir=$1
outdir=$2

if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

for aid in ${QUERY_LIST}
do
    afile=${indir}/${aid}.cm_a7.0
    outfile=${outdir}/${aid}.out
    cat /dev/null > ${outfile}
    for bfile in ${indir}/*.cm_a7.0
    do
        bid=`basename ${bfile} .cm_a7.0`
        zbid=`echo ${bid} | tr -d .`
        score=`MSVNS4MaxCMO 3 10 7467 ${afile} ${bfile} ${aid} ${bid} | awk 'NR == 2'`
        printf '%8s %10d\n' ${zbid} ${score} >> ${outfile}
    done
done

