*=======================================================================
* File:    rdtabn.f
* Author:  Alex Stivala
* Created: July 2008
*
* Read tableau database entry in numeric (Omega matrix) format.
*
* The format of the 'database' is a text file with an entry for each
* structure. The first line of an entry is the identifier and
* order of tableau (i.e. dimension of square Omega matrix), then
* each subsequent row is a row of the Omega matrix, lower triangle
* only (since it is symmetric).
* The diagonal entries are meaningless (self-angle) in tableaux,
* and are included instead to specify the SSE type, with
* the following codes:
*
* 0.000 beta strand
* 1.000 alpha helix
* 2.000 pi helix
* 3.000 3_10 helix
*
* Width of identifier is 8 chars, blank padded on right,
* width of order is 4 digits, blank padded on left.
* There is a single space between identifier and order.
* Each entry in Omega matrix is in radians in [-pi, pi] format 
* F6.3 with a space between each on a line, and one line
* per row of matrix.
*
* E.g.:
*
* D1UBIA_    6
*  0.000 
*  2.650  0.000
* -1.170  2.150  1.000
*  2.040 -1.140  2.080  0.000
* -1.260  1.560 -1.110  2.990  0.000
* -0.590  2.100 -1.230  2.570 -0.720  0.000
* 
* $Id: rdtabn.f 1879 2008-09-09 06:11:24Z astivala $
*=======================================================================

      subroutine rdtabn(nunit, omega, ldo, n, name)

      implicit none
*
* rdtdbn - read a numeric tableau (Omega matrix)
*
*
*     .. Scalar Arguments ..
      integer nunit,ldo,n
      character name*8
*     ..
*     .. Array Arguments
      double precision omega(ldo, *)
*     ..
*
*
* Arguments
* =========
*
* nunit  (input) INTEGER
*        Unit number to read tableau numeric matrix from
*
* omega (output) DOUBLE PRECISION array, dimension (ldo, ldo)
*        Omega matrix read from the file. Must be dimensioned large
*        enough to hold dimension (n,n) data.
*        The full symmetric matrix is stored.
*
* ldo   (input) INTEGER
*       Leading dimension of omega array.
*
* n     (output) INTEGER
*       Order of omega (dimension of matrix).
*       if n > ldo then the read failed as the tableau was too large
*       for the supplied array.
*
* name  (output) CHARACTER*8
*       Identifier (PDB or SCOP id) of the tableau
*    
*=======================================================================
*

*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      integer i,j
      real angle
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
      read(nunit, 10, err=99, end=99) name,n
 10   format (a8,1x,i4)
      if (n .gt. ldo) then
         write(*,20) name,n
 20      format ('Tableau id ',a8,' order ',i4,' is too large.')
         return
      endif
      do 50 i = 1, n
         read(nunit,'(a2048)') line
         do 40 j = 1, i
            read(line((j-1)*7+1:(j-1)*7+7), '(f6.3)') angle
            omega(i,j) = dble(angle)
            omega(j,i) = omega(i,j)
 40      continue
 50   continue

 99   continue
      return
      end

