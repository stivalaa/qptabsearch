#!/bin/sh
#
# run discrete (tableaux) version with sparse matrix operations on the
# 8 example folds,  saving results and stderr (for $TIME)
# and running only one at a time.
#
# WARNING: runs for a long time, overwrites results/<name>.soln.out files
#
TIME=/usr/bin/time

$TIME src/tsrchd_sparse < data/d1ubia_.soln.input > results/d1ubia_.tsrchd.tt.soln.out 2> results/d1ubia_.tsrchd.tt.soln.err

$TIME src/tsrchd_sparse < data/d1tt.solntb1.soln.input > results/d1tt.solntb1.tsrchd.tt.soln.out 2> results/d1tt.solntb1.tsrchd.tt.soln.err

$TIME src/tsrchd_sparse < data/d1ae6h1.soln.input > results/d1ae6h1.tsrchd.tt.soln.out 2> results/d1ae6h1.tsrchd.tt.soln.err

$TIME src/tsrchd_sparse < data/d1bhne_.soln.input > results/d1bhne_.tsrchd.tt.soln.out 2> results/d1bhne_.tsrchd.tt.soln.err

$TIME src/tsrchd_sparse < data/d1h6rb_.soln.input > results/d1h6rb_.tsrchd.tt.soln.out 2> results/d1h6rb_.tsrchd.tt.soln.err

$TIME src/tsrchd_sparse < data/d2phlb1.soln.input > results/d2phlb1.tsrchd.tt.soln.out 2> results/d2phlb1.tsrchd.tt.soln.err

$TIME src/tsrchd_sparse < data/d1tima_.soln.input > results/d1tima_.tsrchd.tt.soln.out 2> results/d1tima_.tsrchd.tt.soln.err

$TIME src/tsrchd_sparse < data/d1f6dc_.soln.input > results/d1f6dc_.tsrchd.tt.soln.out 2> results/d1f6dc_.tsrchd.tt.soln.err

