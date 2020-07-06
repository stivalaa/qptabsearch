#!/bin/sh
#
# run different options on d1tttb1 (serially, using only one cpu at a time)
# saving results and stderr (for time)
#

TIME=/usr/bin/time

$TIME src/tsrchn_sparse < data/d1tttb1.omega.input > results/d1tttb1.tsrchn.tt.out 2> results/d1tttb1.tsrchn.tt.err

$TIME src/tsrchn_sparse < data/d1tttb1.omega.ff.input > results/d1tttb1.tsrchn.ff.out 2> results/d1tttb1.tsrchn.ff.err

$TIME src/tsrchn_sparse < data/d1tttb1.omega.tf.input > results/d1tttb1.tsrchn.tf.out 2> results/d1tttb1.tsrchn.tf.err

$TIME src/tsrchd_sparse < data/d1tttb1.input > results/d1tttb1.tsrchd.tt.out 2> results/d1tttb1.tsrchd.tt.err

$TIME src/tsrchd_sparse < data/d1tttb1.ff.input > results/d1tttb1.tsrchd.ff.out 2> results/d1tttb1.tsrchd.ff.err

$TIME src/tsrchd_sparse < data/d1tttb1.tf.input > results/d1tttb1.tsrchd.tf.out 2> results/d1tttb1.tsrchd.tf.err

