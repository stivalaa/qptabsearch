#!/bin/sh
#
# run different options on d1ubia_ (serially, using only one cpu at a $TIME)
# saving results and stderr (for $TIME)
#
TIME=/usr/bin/time

$TIME ../src/tsrchn_sparse < ../data/d1ubia_.omega.input > ../results/d1ubia_.tsrchn.tt.out 2> ../results/d1ubia_.tsrchn.tt.err

$TIME ../src/tsrchn_sparse < ../data/d1ubia_.omega.ff.input > ../results/d1ubia_.tsrchn.ff.out 2> ../results/d1ubia_.tsrchn.ff.err

$TIME ../src/tsrchn_sparse < ../data/d1ubia_.omega.tf.input > ../results/d1ubia_.tsrchn.tf.out 2> ../results/d1ubia_.tsrchn.tf.err

$TIME ../src/tsrchd_sparse < ../data/d1ubia_.input > ../results/d1ubia_.tsrchd.tt.out 2> ../results/d1ubia_.tsrchd.tt.err

$TIME ../src/tsrchd_sparse < ../data/d1ubia_.ff.input > ../results/d1ubia_.tsrchd.ff.out 2> ../results/d1ubia_.tsrchd.ff.err

$TIME ../src/tsrchd_sparse < ../data/d1ubia_.tf.input > ../results/d1ubia_.tsrchd.tf.out 2> ../results/d1ubia_.tsrchd.tf.err

