#!/usr/bin/env python
###############################################################################
#
# pairwise2matrix.py - convert list of pairwise scores to distance matrix
#
# File:    pairwise2matrix.py
# Author:  Alex Stivala
# Created: September 2008
#
#
# Given list of pairwise scores from qptabmatch_allpairs.py output,
# convert to distance matrix by reformatting as matrix and normalizing.
#
# Usage:
#     pairwise2matrix.py db_directory < scoresfile > matrixfile
#
# Input is from stdin, tab-delimited in the format
# 
#   pdbid1 pdbid2 score
#
# Also reads .tableuxdistmatrix files from db_directory, as used by
# qptabmatch_allpairs.py, to get orders of tableaux (number of SSEs)
# require for normalization.
#
# Output is to stdout, a (square symmetric) matrix of normalized scores
# between each structure. It is in format suitable for reading in R with
# read.table(filename, header=TRUE)
# i.e. the first line is space-delimited identifiers for each column,
# then the matrix (space delimited)
#
# Requires Numeric library (builds Numeric.array)
#
# $Id: pairwise2matrix.py 2703 2009-07-27 06:01:05Z astivala $
# 
###############################################################################

import sys,os,glob
import numpy.oldnumeric as Numeric

from norms import norm1,norm2,norm3

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " <db_directory>\n")
    sys.exit(1)


def main():
    """
    main for pairwise2matrix.py
    """
    if len(sys.argv) != 2:
        usage(os.path.basename(sys.argv[0]))

    db_directory = sys.argv[1]

    # get numbers of SSEs from .tableauxdistmatrix files in db_diurectory
    # build dict of {name : size}
    size_dict = {}
    input_list = glob.glob(os.path.join(db_directory, '*.tableaudistmatrix'))
    for dbfile in input_list:
        idline = open(dbfile).readline()
        qid = idline[:8].lstrip().rstrip()
        qsize = int(idline[8:])
        size_dict[qid] = qsize
 
    scoredict = {}
    for line in sys.stdin:
        if line[0] == '#':
            continue
        s = line.split('\t')
        scoredict[(s[0],s[1])] = float(s[2])


    # get list of unique names (identifiers)
    allnames = [k[0] for k in scoredict.keys() ] + \
               [k[1] for k in scoredict.keys() ]
    namedict = {}
    for name in allnames:
        namedict[name] = name
    names = namedict.values()
    n = len(names)
    assert(len(scoredict) == n*(n-1)/2) # n(n-1)/2 comparisons for n things
    
    # build dictionary mapping identfieres (names) to indicies for matrix
    idict = dict([(b,a) for (a,b) in enumerate(names)])

    # build the distance matrix
    distmatrix = Numeric.zeros((n,n),'d')
    for ((name1,name2),score) in scoredict.iteritems():
        i = idict[name1]
        j = idict[name2]
        normscore = norm3(score, size_dict[name1], size_dict[name2])
        distmatrix[i,j] = normscore
        distmatrix[j,i] = normscore

    namestr = reduce(lambda a,b : a+' '+b,  names)
    print namestr
    for i in xrange(n):
        for j in xrange(n):
            sys.stdout.write('%.4f ' % distmatrix[i,j])
        sys.stdout.write('\n')

if __name__ == "__main__":
    main()

