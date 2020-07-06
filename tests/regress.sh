#!/bin/sh
#
# simple regression test for DLPHS1 and DLPHS2 FORTRAN modules
#
# Alex Stivala 6July2008
#
# Usage: regress <Executable> <Valid_output_file>
#
# e.g. regress testp1 testp1.correct.out
#
# Exit status 0 if OK, nonzero on failure
#
# $Id: regress.sh 1690 2008-07-16 04:35:47Z astivala $
#

if [ $# -ne 2 ]; then
    echo "Usage: $0 <Executable> <Valid_output_file>" >&2
    exit 1
fi

program=$1
validout=$2

outfile=/var/tmp/regress$$.out

./${program}  > ${outfile}

if [ $? -ne 0 ]; then
    echo "$0: $program FAILED"
    exit 2
fi

diff ${validout} ${outfile} >/dev/null

if [ $? -ne 0 ]; then
    echo "$0: regression test $program FAILED"
    rc=3
else
    echo "$0: regression test $program OK" 
    rc=0
fi
rm ${outfile}
exit $rc

