#!/bin/bash
###############################################################################
#
# qptabmatchstructs.sh - run the QP tableau matching on two structures
#
# File:    qptabmatchstructs.sh
# Author:  Alex Stivala
# Created: August 2008
#
# Run the QP tableau matching program on two structures in PDB format.
# Outputs to stdout a PyMOL script (.pml file) containing a visualzation
# of the matches (maximally similar substructures coloured) with the match
# socre as a comment, and, if the -s option is given, writes a PDB
# file with structures superimposed in cwd.
#
# Usage:
#     qptabmatchstructs.sh [-skq] [-e sse_num_list]
#                                        struct1.pdb struct2.pdb > output.pml
#
#     -s: write PDB file with superposition to cwd, filenme is 
#         queryid_dbid.pdb. WARNING: overwrites if file exists.
#
#     -k: use the tableau HH and KK codes for antiprallel/parallel
#         strands in same sheet.
#
#     -q: do not use the ordering constraint (allow nonsequential  matchings)
#
#     -u: randomly permute the rows+columns of the struct1 tableau+distmatrix
#         (This is for testing non-sequential tableau matching; so isn't
#         useful except with -q)
#
#     -e sse_num_list: list of SSE sequential numbers to select from
#                      struct1 reather than whole structure
#
# Then in PyMOL use @output.pml to run the PyMOL script.
#
# Uses the Python scripts pytableaucreate.p to create tableaux for input
# to the FORTRAN tsrchd_sparse program, and Python scripts 
# soln2ssemap.py and ssemap2pml.py to process the output of tsrchd_sparse
# into PyMOL script.
#
# Environment variables:
#
#   PATH must contain the location of the Python scripts, ie where this
#   script itself is and the ptgraph/ directory with pytableaucreate.py etc.,
#   and the location of tsrchd_sparse.
#   The dssp program must also be in the PATH.
#
#   PYTHONPATH must contain the directory containing the ptsecstruct.py 
#   and other Python modules used by the Python scripts.
#
# Note you would not use this script to do a large number of comparisons
# such as a database search - for that you would run tsrchd_sparse (or other)
# for a tableau against database of tableaux previously created with
# pytableaucrate.py and buildtableauxdb.py, and soln2ssemap.py | ssemap2pml.py
# pipline on output for desired matches (e.g. top n scores obtained by sorting
# output of tsrchd). This script is a convenient way of running a single
# matching between two structures and examining the result correspondence
# between SSEs.
# 
# $Id: qptabmatchstructs.sh 2042 2008-11-26 23:55:39Z astivala $
# 
###############################################################################


# write tableau and distance matrix for tsrchd to stdout
# Parameters:
#     pdbfile - filename of PDB file
#     use_hk - if 1, use the HH and KK codes
#     extra_opts - extra options to add to pytableaucreate.py
writetableau() {
    pdbfile=$1
    use_hk=$2
    extra_opts="$3"
    tabopts="-b -35 -f -t dssp -p none ${extra_opts}"
    if [ $use_hk -eq 1 ]; then
        tabopts="${tabopts} -k"
    fi
    pytableaucreate.py ${tabopts} ${pdbfile}
}



writepdb=0
use_hk=0
use_ordering=1
sse_num_list=''
sse_num_list_opt=''
randomly_permute=0

while getopts 'skqe:u' opt
do
    case $opt in
    s) writepdb=1
    ;;
    k) use_hk=1
    ;;
    q) use_ordering=0
    ;;
    e) sse_num_list="$OPTARG"
       sse_num_list_opt="-s ${sse_num_list}"
    ;;
    u) randomly_permute=1
    ;;
    ?)
    echo "Usage: $0 [-squk] [-e sse_num_list] struct1.pdb struct2.pdb" >&2
    exit 1
    ;;
    esac
done
shift $(($OPTIND - 1))

if [ $# -ne 2 ]; then
    echo "Usage: $0 [-squk] [-e sse_num_list] struct1.pdb struct2.pdb" >&2
    exit 1
fi

struct1=$1
struct2=$2


tmpfile1=/tmp/qptabsh1$$
tmpfile2=/tmp/qptabsh2$$
tmpfile3=/tmp/qptabsh3$$
tmpfile4=/tmp/qptabsh4$$

writetableau ${struct2} ${use_hk} "" > ${tmpfile2}
echo ${tmpfile2} > ${tmpfile1}  # filename of tableaux database
if [ $use_ordering -eq 0 ]; then
    echo "T F T" >> ${tmpfile1}     # options: type,order,output
else
    echo "T T T" >> ${tmpfile1}     # options: type,order,output
fi
extra_opts="${sse_num_list_opt}"
if [ $randomly_permute -ne 0 ]; then
    extra_opts="${extra_opts} -u"
    writetableau ${struct1} ${use_hk} "${extra_opts}" > ${tmpfile4}
    # get permutation from first line of tableaucreate output
    permutation=`awk 'NR == 1 {print $3}' ${tmpfile4}`
    awk 'NR > 1'  ${tmpfile4} >> ${tmpfile1}
else
    writetableau ${struct1} ${use_hk} "${extra_opts}" >> ${tmpfile1}
fi

qsize=`awk 'NR == 3 {print $2}' < ${tmpfile1}`

if [ -z ${sse_num_list} ]; then
    remap_stage="cat"  # no need to remap SSE nums: stdin to stdout unchanged
else
    remap_stage="ssesubsetremap.py ${sse_num_list}"
fi
if [ $randomly_permute -ne 0 ]; then
    perm_remap_stage="ssepermutationremap.py ${permutation}"
else
    perm_remap_stage="cat"
fi

if [ $writepdb -eq 1 ]; then
    tsrchd_sparse < ${tmpfile1} | soln2ssemap.py -q ${qsize} | ${perm_remap_stage}  | ${remap_stage} | tee ${tmpfile3} |  ssemap2pml.py -s -u ${struct1} -b ${struct2}
    superimposessemap.py -u ${struct1} -b ${struct2} -o . < ${tmpfile3} >/dev/null
else
    tsrchd_sparse < ${tmpfile1} | soln2ssemap.py -q ${qsize} | ${perm_remap_stage} | ${remap_stage} | ssemap2pml.py -s -u ${struct1} -b ${struct2}
fi

rm -f ${tmpfile1} ${tmpfile2} ${tmpfile3} ${tmpfile4}

