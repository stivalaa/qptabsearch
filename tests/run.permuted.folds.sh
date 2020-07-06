#!/bin/sh
#
# run.permuted.folds.sh
#
# run tsrchd_sparse
# on the permuted queries
# saving results and stderr (for $TIME)
# in cwd (so run from results/permuted/ dir)
# and running one at a time.
#
# $Id: run.permuted.folds.sh 2148 2009-03-28 00:55:29Z astivala $
#

TIME=/usr/bin/time
TSRCHD=${HOME}/phd/qptabsearch/src/tsrchd_sparse
INPUT_DIR=${HOME}/phd/qptabsearch/data/permuted

for queryfile in ${INPUT_DIR}/*.input
do
   query=`basename ${queryfile} .input`
   $TIME ${TSRCHD} < ${queryfile} > ${query}.out 2> ${query}.err
done


