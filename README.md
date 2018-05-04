Tableau-based protein substructure search using quadratic programming
# Tableau-based protein substructure search using quadratic programming

Imported from https://sites.google.com/site/alexdstivala/home/qpprotein

This software can be used freely for any purpose, modified, redistributed, etc.
with no restrictions. However we would appreciate it if you acknowledge
your use of it, and in particular if you would cite our paper
in any publication that makes use of it.

QP Tableau Search is also available for online use via the [Pro-origami web server](http://munk.csse.unimelb.edu.au/pro-origami).

## Tableau databases

The tableau databases were built with 
buildtableauxdb.py -t dssp -3 -5 -p none and converted to the ASCII
format available here with convdbnumeric2ascii.py and
convdbpacked2ascii.py for the numeric and tableaux formats
respectively.  The distance matrices were built with 
buildtableauxdb.py -t dssp -3 -5 -p none -d and the combined
tableaux and distance matrix database in ASCII format created with
the convdb2.py script.
These scripts are included with the source code.


## Randomly permuted tableaux

In the "Non-linear matchings" section of our paper, we describe the use
of random permutations of tableaux as an artificial test to verify the
capability of our method to find non-linear matchings, i.e. sets of
correspondences between SSEs in which the sequential order of corresponding
SSEs is not preserved. Here we provide some more details on this method
for those who might like to use to evaluate other non-linear structural
alignment techniques.

As described in the paper (and citations therein) a tableau is a square
symmetric matrix in which each row represents the orientation of one SSE
relative to every other SSE. Normally, the rows (and columns) are in the
same order as the SSEs in the protein sequence, from N to C terminus.
In order to generate randomly permuted tableaux, we instead generate
the tableaux by randomly permuting the sequence of SSEs, so that rather
than the rows being ordered according to the sequence, they are ordered
according to the random permutation that was generated.

The -u option on the pytableaucreate.py script included 
with the source code is used to perform this task. It also outputs the
permutation used, which is necessary so that any alignment based on the
permuted tableau can be mapped back to the actual SSEs in the structure
which the row and column in the tableau represents. The 
ssepermutationremap.py script accomplishes this task, and is
used by the qptabmatchstructs.sh script when the -u option
is used to randomly permute the query structure. This script illustrates
the pipeline used to perform an alignment (possibly with a permuted tableau)
and visualize it using [PyMOL](http://www.pymol.org).
The genpermutedqueries.sh script was used to generate the queries
described in the paper.

All of the scripts mentioned here are provided with the source code download,
and all contain internal documentation as to their purpose and usage. 

## Faster implementation

As described in our paper, all the results published there were with
tsrchd_sparse
the sparse matrix implementation of the QP solver using
the UMFPACK solver in
[SuiteSparse](http://www.cise.ufl.edu/research/sparse/SuiteSparse/).
However, if you have the
[Intel Math Kernel Library](http://software.intel.com/en-us/intel-mkl/)
(MKL) version 10.1, which contains the 
[PARDISO](http://www.pardiso-project.org/) sparse linear solver,
you can use the tsrchd_pardiso implementation, which is approximately
70% faster than tsrchd_sparse.

## Reference

If you use our software, data, or results in your research, please cite:

Stivala, A., Wirth, A. and Stuckey, P.,
[Tableau-based protein substructure search using quadratic programming](http://www.biomedcentral.com/1471-2105/10/153/abstract)BMC Bioinformatics 2009, 10:153

