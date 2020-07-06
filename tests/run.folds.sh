#!/bin/sh
#
# run discrete (tableaux) version with sparse matrix operations on the
# 8 example folds,  saving results and stderr (for $TIME)
# and running only one at a time.
#
# WARNING: runs for a long time, overwrites results/full/<name>.out files
#
TIME=/usr/bin/time

$TIME src/tsrchd_sparse < data/d1ubia_.input > results/full/d1ubia_.tsrchd.tt.out 2> results/full/d1ubia_.tsrchd.tt.err

$TIME src/tsrchd_sparse < data/d1tttb1.input > results/full/d1tttb1.tsrchd.tt.out 2> results/full/d1tttb1.tsrchd.tt.err

$TIME src/tsrchd_sparse < data/d1ae6h1.input > results/full/d1ae6h1.tsrchd.tt.out 2> results/full/d1ae6h1.tsrchd.tt.err

$TIME src/tsrchd_sparse < data/d1bhne_.input > results/full/d1bhne_.tsrchd.tt.out 2> results/full/d1bhne_.tsrchd.tt.err

$TIME src/tsrchd_sparse < data/d1h6rb_.input > results/full/d1h6rb_.tsrchd.tt.out 2> results/full/d1h6rb_.tsrchd.tt.err

$TIME src/tsrchd_sparse < data/d2phlb1.input > results/full/d2phlb1.tsrchd.tt.out 2> results/full/d2phlb1.tsrchd.tt.err

$TIME src/tsrchd_sparse < data/d1tima_.input > results/full/d1tima_.tsrchd.tt.out 2> results/full/d1tima_.tsrchd.tt.err

$TIME src/tsrchd_sparse < data/d1f6dc_.input > results/full/d1f6dc_.tsrchd.tt.out 2> results/full/d1f6dc_.tsrchd.tt.err

