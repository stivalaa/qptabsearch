*=======================================================================
* File:    rdtabd.f
* Author:  Alex Stivala
* Created: July 2008
*
* Read tableau database entry in discrete tableaux format.
*
* The format of the 'database' is a text file with an entry for each
* structure.
* The first line of an entry is the identifier and
* order of tableau (i.e. dimension of square array), then
* each subsequent row is a row of the tableau, lower triangle
* only (since it is symmetric).
* The diagonal entries are meaningless (self-angle) in tableaux,
* and are included instead to specify the SSE type, with
* the following codes:
*
* e     beta strand
* xa    alpha helix
* xi    pi helix
* xg    3_10 helix
*
* Width of identifier is 8 chars, blank padded on right,
* width of order is 4 digits, blank padded on left.
* There is a single space between identifier and order.
* Each entry in tableau is two characters, with a space betwen
* each on a line, and one line
* per row of matrix.
*
* E.g.:
*
* /local/charikar/astivala/tableauxdb/astral/tableauxdb.ascii
*  T F
* D1UBIA_    6
* e  
* OT e  
* LE RT xa 
* RT LE RT e  
* LE RD LE OT e  
* PE RT LE OT PE e  
*
* The tableau is followed by a SSE distance matrix, not read
* by this subroutine, rather handled by a call to subroutine RDISTM.
*
* 
* $Id: rdtabd.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================

      subroutine rdtabd(nunit, tab, ldo, n, name)

      implicit none
*
* rdtdbn - read a (discrete) tableau
*
*
*     .. Scalar Arguments ..
      integer nunit,ldo,n
      character name*8
*     ..
*     .. Array Arguments
      character*2 tab(ldo, *)
*     ..
*
*
* Arguments
* =========
*
* nunit  (input) INTEGER
*        Unit number to read tableau numeric matrix from
*
* tab (output) CHARACTER*2 array, dimension (ldo, ldo)
*        Tableau read from the file. Must be dimensioned large
*        enough to hold dimension (n,n) data.
*        The full symmetric tableau is stored.
*
* ldo   (input) INTEGER
*       Leading dimension of tab array.
*
* n     (output) INTEGER
*       Order of tab (dimension of array).
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
      character line*2048
      character tcode*2
*     ..
*     .. Local Arrays ..
*     ..
*     .. Intrinsic Functions ..
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
            read(line((j-1)*3+1:(j-1)*3+2), '(a2)') tcode
            tab(i,j) = tcode
            tab(j,i) = tcode
 40      continue
 50   continue

 99   continue
      return
      end

