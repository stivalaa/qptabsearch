#!/usr/bin/env python
###############################################################################
#
# msvns4maxcmo_allall.py - run MSVNS4MaxCMO with all contact maps in dir 
#
# File:    msvns4maxcmo_allall.py
# Author:  Alex Stivala
# Created: September 2008
#
#
# Run MSVNS4MaxCMO (Pelta et al 2008 BMC Bioinformatics 9:161)
# http://modo.ugr.es/jrgonzalez/msvns4maxcmo
# aginst all contact maps in directory such
# as that created by build_fischer_cm.sh.
# between each)
#
# This is _allall since for Fischer data set the db is actually just
# all the contact maps in directory, so it is doing
# duplicate comparisons so for n (=68) tableux it does
# n*n (=4624) comparisons.
# 
#
# Usage:
#     msvns4maxcmo_allall.py query_directory results_directory
#
# query_directory is the directory containing contact map files,
# as built with build_fischer_cm.sh for example.
#
#
# results_dirctory is a directory to write the output to.
# Each query (contact map file) results in one file created
# with .out suffix
# in the output directory, containing the results from that query
# against all the others.
# The output  is  in the format with identifier and MSVNS4MaxCMO
# score on each line separated by whitespace.
# WARNING: these files overwritten if they exist.
# results_directory is created if it does not exist.
#
# Environment variables:
#
#   PATH must contain the location of MSVNS4MaxCMO
#
# $Id: msvns4maxcmo_allall.py 1955 2008-10-06 23:48:52Z astivala $
# 
###############################################################################

import sys,os,glob

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " <query_directory>  <results_directory>\n")
    sys.exit(1)

    
def main():
    """
    main for msvns4maxcmo.py
    """
    if len(sys.argv) != 3:
        usage(os.path.basename(sys.argv[0]))

    query_directory = sys.argv[1]
    results_directory = sys.argv[2]

    if not os.path.exists(results_directory):
        os.mkdir(results_directory)
    elif not os.path.isdir(results_directory):
        sys.stderr.write('%s is not a directory\n' % results_directory)
        sys.exit(1)

    input_list = glob.glob(os.path.join(query_directory, '*.cm_a7.0'))
    i = 0
    while i < len(input_list):
        afile = input_list[i]
        aid = os.path.basename(afile).partition('.')[0]
        outfile = os.path.join(results_directory, aid+'.out')
        outfh = open(outfile, 'w')
        j = 0
        while j < len(input_list):
            bfile = input_list[j]
            bid = os.path.basename(bfile).partition('.')[0]
            (msvns_in, msvns_out) = \
              os.popen2(['MSVNS4MaxCMO', '3', '10', '7467', afile, bfile,
                         aid, bid])
            msvns_in.close()
            lineno = 1
            for line in msvns_out:
              if lineno == 2:
                  score = int(line)
                  break
              lineno += 1
            msvns_out.close()
            outfh.write('%8s %10d\n' % (bid, score))
            j += 1
        outfh.close()
        i += 1

if __name__ == "__main__":
    main()

