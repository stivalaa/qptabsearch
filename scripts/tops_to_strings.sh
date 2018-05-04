#!/bin/sh
#
# File:    tops_to_strings.sh
# Author:  Alex Stivala
# Created: March 2009
#
# tops_to_strings.sh - build TOPS strings from database of TOPS files
#
# Usage: tops_to_strings.sh tops_db_dir 
#
#   top_db_dir is the name of the directory containing TOPS file
#   (as built with build_tops_files)
#
# The TOPS strings, one per line, are written to stdout
#
# $Id: tops_to_strings.sh 2122 2009-03-23 22:39:40Z astivala $

# location of tops_comparison directory, contains jars/translation.jar etc.
TOPS_COMPARISON_ROOT=/local/charikar/astivala/biosoftware/tops_comparison

if [ $# -ne 1 ]; then
    echo "Usage: $0 tops_db_dir" 2>&1
    exit 1
fi

tops_db_dir=$1


for topsfile in `find $tops_db_dir -maxdepth 1 -name \*.tops`
do
    sid=`basename $topsfile .tops`
    java -cp ${TOPS_COMPARISON_ROOT}/jars/translation.jar tops.translation.Tops2String $topsfile $sid
done
