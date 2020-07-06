#!/bin/sh
#
# run different options on 2JMK (serially, using only one cpu at a time)
# saving results and stderr (for time)
#

TIME=/usr/bin/time

$TIME src/tsrchd_sparse < data/2jmk.tf.input > results/2jmk.tsrchd.tf.out 2> results/2jmk.tsrchd.tf.err

$TIME src/tsrchd_sparse < data/2jmk.tt.input > results/2jmk.tsrchd.tt.out 2> results/2jmk.tsrchd.tt.err

