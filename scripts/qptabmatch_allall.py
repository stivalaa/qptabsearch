#!/usr/bin/env python
###############################################################################
#
# qptabmatch_allall.py - run the QP tableau matching with all tableauxdistmatrix
#                        files in a directory against a db (file of multiple
#                        tableaux+distmatrices)
#
# File:    qptabmatch_allall.py
# Author:  Alex Stivala
# Created: September 2008
#
#
# Run QP tableau match on all tableaux in a directory against database
# of tableaux such 
# as that created by build_fischer_db.sh. Each file in the directory
# has a .tableaudistmatrix suffix and contains the header line with
# identifeir and dimension, then tabelau and SSE distance matrix
# (both lower triangle fortran format
# for tsrchrd_sparse etc.)
# The db is is many such files concantentaed together (with blank line
# between each)
#
# This is _allall since for Fischer data set the db is actually just
# all the tableaux in directory, so it is doing
# duplicate comparisons so for n (=68) tableux it does
# n*n (=4624) comparisons.
# 
#
# Usage:
#     qptabmatch_allall.py query_directory dbfile results_directory
#
# query_directory is the directory containing .tableaudistmatrix files,
# as built with build_fischer_db.sh for example.
#
# db_file is the 'database' file of tableaux+distmatrices.
#
# results_dirctory is a directory to write the output to.
# Each query (.tableauxdistmatrix file) results in one file created
# with .out suffix
# in the output directory, containing the results from that query
# against the db.
# Each file is created by tsrchd_sparse (see tsrchd.f for output format):
# each line is query identifier and score (whitespace delimited).
# WARNING: these files overwritten if they exist.
# results_directory is created if it does not exist.
#
# Environment variables:
#
#   PATH must contain the location of tsrchd_sparse.
#
#
# $Id: qptabmatch_allall.py 1894 2008-09-15 08:37:50Z astivala $
# 
###############################################################################

import sys,os,glob

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " <query_directory> <db_file> <results_directory>\n")
    sys.exit(1)

    
def main():
    """
    main for qptabmatch_allall.py
    """
    if len(sys.argv) != 4:
        usage(os.path.basename(sys.argv[0]))

    query_directory = sys.argv[1]
    db_file  = sys.argv[2]
    results_directory = sys.argv[3]

    if not os.path.exists(results_directory):
        os.mkdir(results_directory)
    elif not os.path.isdir(results_directory):
        sys.stderr.write('%s is not a directory\n' % results_directory)
        sys.exit(1)

    input_list = glob.glob(os.path.join(query_directory, '*.tableaudistmatrix'))
    i = 0
    while i < len(input_list):
        qfile = input_list[i]
        qid = open(qfile).readline()[:8].lstrip().rstrip()
        outfile = os.path.join(results_directory, 
                        os.path.splitext(os.path.basename(qfile))[0] + '.out' )
        tsrchd_in = os.popen('tsrchd_sparse >'+ outfile, 'w')
        tsrchd_in.write(db_file + '\n')   # name of db file
        tsrchd_in.write('T T F\n') # LTYPE LORDER LSOLN
        tsrchd_in.write(open(qfile).read()) # tableau+distmatrix of qfile
        tsrchd_in.close()
        i += 1

if __name__ == "__main__":
    main()

