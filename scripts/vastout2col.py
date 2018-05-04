#!/usr/bin/env python
#
# File:    vastout2col.sh
# Author:  Alex Stivala
# Created: November 2008
#
# vastout2col.sh - Convert VAST .gibbs output format to 2 column
#                   format as output by tsrchd_sparse etc. which can
#                   be processed with tsevalfn.py etc.
#
# Usage: vastout2col.sh < domain.gibbs
# 
# Output has two columns, database id and VAST Pcli score
#
# Output is to stdout.
#
# Uses the output format from  VAST (Gibrat et al 1996; Madej et al 1995),
# available from
# http://migale.jouy.inra.fr/outils/mig/vast/
#
# $Id: vastout2col.py 2032 2008-11-23 05:00:03Z astivala $
#

import os,sys

value_header = False
dbid = None

for line in sys.stdin:
    splitline = line.split()
    if len(splitline) > 1 and splitline[1] == 'Nclique=':
        dbid = splitline[0]
        value_header = False
    elif splitline[0] == 'Nres' and splitline[6] == 'Pcli':
        value_header = True
    elif value_header:
        Pcli = splitline[6]
        sys.stdout.write('%s    %s\n' % (dbid, Pcli))
        value_header = False

