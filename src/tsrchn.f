*=======================================================================
* File:    tsrchn.f
* Author:  Alex Stivala
* Created: July 2008
*
* FORTRAN-77 implementation of tableau searching numeric (ie Omega
* matrix rather than actual tableau) using the DLSOQP implementation
* of Prof. Yinyu Ye's interior point QP solver.
*
* The 'database' to search is an ASCII file of numeric tableaux
* (Omega matrices) in format described in rdtabn.f.
*
* The results are printed to stdout as 
*
* name score
*
* just as in tabesarchqpml.py, the original Python implementation.
*
* This file is a program not a subroutine.
* Both the name of the database file to read, and the actual
* query tableau are read from stdin. 
* The first line is the name
* of the database file.
* The second line is for options. There are currently 3 logical
* options, for SSE type constraint (only allow SSEs of same type ot
* match) and ordering constraint (disallow out of sequence order 
* matches). The 3rd is the option to output solution x vector values,
* not just score.
* They are single character logical values (T or F).
* First is type, second is order, separated by one space.
*
* The subsequent lines are a single tableau in the same format as
* each entry in the database i.e.:
*
* The first line of an entry is the identifier and
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
* /local/charikar/astivala/tableauxdb/astral/tableauxdb.numeric.ascii
*  T F F
* D1UBIA_    6
*  0.000 
*  2.650  0.000
* -1.170  2.150  1.000
*  2.040 -1.140  2.080  0.000
* -1.260  1.560 -1.110  2.990  0.000
* -0.590  2.100 -1.230  2.570 -0.720  0.000
*
*
* 
* $Id: tsrchn.f 2064 2009-02-24 06:11:28Z astivala $
*=======================================================================

      program tsrchn

      implicit none

*     .. Parameters ..
*     largest tableau (omega matrix) dimension allowed
      integer maxdim
      parameter(maxdim = 110)
*     unit number of standard input
      integer inunit
      parameter(inunit = 5)
*     unit number for database file
      integer dbunit
      parameter(dbunit=10)
*     ..
*     .. Local Scalars ..
      integer qn,dbn,info,j
      character dbfile*1024,qid*8,dbid*8
      logical ltype,lorder,lsoln
      double precision score
*     ..
*     .. Local Arrays ..
      double precision qtab(maxdim,maxdim),dbtab(maxdim,maxdim)
      double precision soln(maxdim*maxdim+maxdim+maxdim)
      
*     ..
*     .. Intrinsic Functions ..
!      intrinsic Time,CTime
*     ..
*     .. External Subroutines and Functions ..
      external tmatn,rdtabn

*     ..
*     .. Executable Statements ..
*
      read (*,'(a1024)') dbfile
      read (*,'(l1,1x,l1,1x,l1)') ltype,lorder,lsoln
      call rdtabn(inunit, qtab, maxdim, qn, qid)

      write(*, 5) ltype,lorder
 5    format('# TSRCHN LTYPE = ',l1,' LORDER = ',l1, ' LSOLN = ',l1)
      write(*, 6) qid
 6    format('# QUERY ID = ',a8)
      write(*, 7) dbfile
 7    format('# DBFILE = ',a80)
!      write(*,8) CTime(Time())
! 8    format('# ',a24)

      open(dbunit, file=dbfile, status='OLD')
 10   if (.true.) then
         call rdtabn(dbunit, dbtab, maxdim, dbn, dbid)
         call tmatn(qtab, maxdim, qn, dbtab, maxdim, dbn, 
     $        ltype, lorder,
     $        score, soln, info)
         if (info .ne. 0) then
            write(*,20) dbid,info
 20         format(a8,1x,'ERROR, INFO = ',i4)
         else
            write(*,30) dbid,score
 30         format(a8,f12.4)
            if (lsoln) then
               write(*,40), (soln(j), j = 1, qn*dbn+qn+dbn)
 40            format (f8.4)
            endif
         endif
         read(dbunit, '(a80)', err=999, end=999)
         goto 10
      endif
 999  continue
      close(dbunit)
      end
      
