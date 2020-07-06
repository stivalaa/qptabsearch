#!/bin/sh
#
# run discrete (tableaux) version with sparse matrix operations on the
# 8 example folds,  saving results and stderr (for $TIME)
# and running only one at a time.
#
# WARNING: runs for a long time, overwrites results/<name>.out files
#
TIME=/usr/bin/time

$TIME src/tsrchd_sparse < data/d1ubia_.tf.input > results/d1ubia_.tsrchd.tf.out 2> results/d1ubia_.tsrchd.tf.err

$TIME src/tsrchd_sparse < data/d1tttb1.tf.input > results/d1tttb1.tsrchd.tf.out 2> results/d1tttb1.tsrchd.tf.err

$TIME src/tsrchd_sparse < data/d1ae6h1.tf.input > results/d1ae6h1.tsrchd.tf.out 2> results/d1ae6h1.tsrchd.tf.err

$TIME src/tsrchd_sparse < data/d1bhne_.tf.input > results/d1bhne_.tsrchd.tf.out 2> results/d1bhne_.tsrchd.tf.err

$TIME src/tsrchd_sparse < data/d1h6rb_.tf.input > results/d1h6rb_.tsrchd.tf.out 2> results/d1h6rb_.tsrchd.tf.err

$TIME src/tsrchd_sparse < data/d2phlb1.tf.input > results/d2phlb1.tsrchd.tf.out 2> results/d2phlb1.tsrchd.tf.err

$TIME src/tsrchd_sparse < data/d1tima_.tf.input > results/d1tima_.tsrchd.tf.out 2> results/d1tima_.tsrchd.tf.err

$TIME src/tsrchd_sparse < data/d1f6dc_.tf.input > results/d1f6dc_.tsrchd.tf.out 2> results/d1f6dc_.tsrchd.tf.err

