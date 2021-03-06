#!/usr/bin/env python
###############################################################################
#
# normalize_msvns4maxcmo.py - normalize MSVNS4MACMO output scores by size
#
# File:    normalize_msvns4maxcmo.py
# Author:  Alex Stivala
# Created: September 2008
#
#
# Given list of scores from MSVNS4MAXCMO (Pelta et al 2008)
# run with msvns4maxcmo_allall.py etc.
# normalize by protein size.
#
# Usage:
#     normalize_msvns4maxcmo.py normtype queryid db_directory < scoresfile > matrixfile
#
# normtype is the normalization function: 1,2, or 3
# queryid is the identifier for the query structure, all the scores are
# assumed to be for matching this to the identifier on each line of the input.
# Input is from stdin, tab-delimited in the format
# 
#   pdbid score
#
# Also reads .cm_a7.0 files from db_directory, as used by MSVNS4MAXCMO
# to get number of contacts in the contact map.
#
# Output is to stdout, in the same format as input i.e.
#
#  pdbid normscore
#
# where normscore is the normalized score.
#
# $Id: normalize_msvns4maxcmo.py 1964 2008-10-08 06:20:53Z astivala $
# 
###############################################################################

import sys,os,glob

from norms import norm1,norm2,norm3

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " <normtype> <queryid> <db_directory>\n")
    sys.exit(1)


def main():
    """
    main for normalize_tabmatch.py
    """
    if len(sys.argv) != 4:
        usage(os.path.basename(sys.argv[0]))

    normtype = int(sys.argv[1])
    queryid = sys.argv[2]
    db_directory = sys.argv[3]
    if normtype < 1 or normtype > 3:
        sys.stderr.write('normtype must be 1,2, or 3\n')
        usage(sys.argv[0])

    # get numbers of SSEs from .tableauxdistmatrix files in db_diurectory
    # build dict of {name : size}
    size_dict = {}
    input_list = glob.glob(os.path.join(db_directory, '*.cm_a7.0'))
    for dbfile in input_list:
        splitfile = os.path.basename(dbfile).split('.')
        qid = reduce(lambda a,b:a+b, splitfile[:len(splitfile)-2]).upper()
        if not qid.isdigit(): # if not 'compressed' CATH id (periods removed)
            qid = reduce(lambda a,b:a+'.'+b, splitfile[:len(splitfile)-2]).upper()
        firstline = open(dbfile).readline()
        num_residues = int(firstline)
        num_contacts = 0
        lineno = 0
        fh = open(dbfile)
        for line in fh:
            if lineno == 0:
                lineno += 1
                continue # first line is number of residues
            if line.lstrip().rstrip() == '':
                lineno += 1
                continue
            sline = line.split()
            if sline[0].isdigit() and sline[1].isdigit():
                num_contacts += 1
            lineno += 1
        fh.close()
        # size is number of contacts in contact map, not number of residues
        size_dict[qid] = num_residues
        # FIXME NB  this is SUPPOSED to be num_contacts, but actually
        # num_residues works better and is closer to Pelta et al 2008 
        # published results when measured for accuracy with AUC
#        size_dict[qid] = num_contacts
 
    # build list of (name , score)
    scorelist = []
    commentlines = []
    for line in sys.stdin:
        if line[0] == '#':
            commentlines.append(line)
            continue
        s = line.split()
        scorelist.append((s[0], float(s[1])))

    querysize = size_dict[queryid.upper()]

    sys.stdout.write('# normalization type ' + str(normtype) + '\n')
    sys.stdout.write('# generated by ' + os.path.basename(sys.argv[0]) +'\n')
    sys.stdout.write('# from:\n')
    for line in commentlines:
        sys.stdout.write(line)
    sys.stdout.write('#\n')

    for (pdbid, score) in scorelist:
        dsize = size_dict[pdbid.upper()]
        if normtype == 1:
            normscore = norm1(score, querysize, dsize)
        elif normtype == 2:
            normscore = norm2(score, querysize, dsize)
        elif normtype == 3:
            normscore = norm3(score, querysize, dsize)
        else:
            raise ValueError('unknown norm type ' + str(normtype) + '\n')
        sys.stdout.write('%s\t%20.8f\n' % (pdbid, normscore))
    

if __name__ == "__main__":
    main()

