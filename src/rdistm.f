*=======================================================================
* File:    rdistm.f
* Author:  Alex Stivala
* Created: August 2008
*
* Read numeric SSE distance matrix.
*
* Each row is a row of the distance matrix, lower triangle
* only (since it is symmetric).
* The diagonal entries are meaningless (self-distance)
* and are included instead to specify the SSE type, with
* the following codes:
*
* 0.000 beta strand
* 1.000 alpha helix
* 2.000 pi helix
* 3.000 3_10 helix
*
* Each entry in matrix is in Angstroms format
* F6.3 with a space between each on a line, and one line
* per row of matrix.
*
* E.g.:
*
*   0.000 
*   4.501  0.000 
*  11.662 10.386  1.000 
*  16.932 17.644  9.779  3.000 
*  10.588 13.738 11.815 10.527  0.000 
*  15.025 18.692 17.143 15.341  6.466  0.000 
*  15.298 17.276 16.276 20.075 13.264 11.610  3.000 
*   7.549 11.072 12.248 12.446  4.583  9.903 15.689  0.000 
*
* 
* This subroutine is used to read the distance matrix from the 
* database, where the identifier and order and tableau have already
* been read by subroutine RDTABD.
*
* $Id: rdistm.f 1879 2008-09-09 06:11:24Z astivala $
*=======================================================================

      subroutine rdistm(nunit, dmat, ldo, n)

      implicit none
*
* rdistm - read a numeric SSE distance matrix
*
*
*     .. Scalar Arguments ..
      integer nunit,ldo,n
*     ..
*     .. Array Arguments
      double precision dmat(ldo, *)
*     ..
*
*
* Arguments
* =========
*
* nunit  (input) INTEGER
*        Unit number to read numeric distance matrix from
*
* dmat  (output) DOUBLE PRECISION array, dimension (ldo, ldo)
*        Distance matrix read from the file. Must be dimensioned large
*        enough to hold dimension (n,n) data.
*        The full symmetric matrix is stored.
*
* ldo   (input) INTEGER
*       Leading dimension of distance matrix array.
*
* n     (input) INTEGER
*       Order of distance matrix (dimension of matrix).
*       n <= ldo
*
*=======================================================================
*

*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      integer i,j
      real dist
      character line*2048
*     ..
*     .. Local Arrays ..
*     ..
*     .. Intrinsic Functions ..
      intrinsic dble
*     ..
*     .. External Subroutines and Functions ..
*     ..
*     .. Executable Statements ..
*     
      do 50 i = 1, n
         read(nunit,'(a2048)') line
         do 40 j = 1, i
            read(line((j-1)*7+1:(j-1)*7+7), '(f6.3)') dist
            dmat(i,j) = dble(dist)
            dmat(j,i) = dmat(i,j)
 40      continue
 50   continue

 99   continue
      return
      end

