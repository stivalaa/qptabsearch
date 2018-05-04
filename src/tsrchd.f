*=======================================================================
* File:    tsrchd.f
* Author:  Alex Stivala
* Created: July 2008
*
* FORTRAN-77 implementation of tableau searching discrete (tableau)
* the SOLVQP implementation of Prof. Yinyu Ye's interior point QP
* solver.
*
* The 'database' to search is an ASCII file of  tableaux
* (Omega matrices) in format described in rdtabd.f.
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
* matches). The third is to output not just the scores but also solution
* vector values.
* They are single character logical values (T or F).
* First is type, second is order, third is solution output,
* separated by one space.
*
* The subsequent lines are a single tableau in the same format as
* each tableau entry in the database i.e.:
*
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
* Following the tableau is the distance matrix.
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
* 
* E.g.:
* 
* /local/charikar/astivala/tableauxdb/astral/tableauxdistmatrixdb.ascii
*  T T F
* D1UBIA_    8
* e  
* OT e  
* LE RT xa 
* PD OS RD xg 
* RT LE RT LS e  
* LE RD LE LS OT e  
* RT LS LS RD PE OS xg 
* PE RT LE RD OT PE RT e  
*  0.000 
*  4.501  0.000 
*  1.662 10.386  1.000 
* 16.932 17.644  9.779  3.000 
* 10.588 13.738 11.815 10.527  0.000 
* 15.025 18.692 17.143 15.341  6.466  0.000 
* 15.298 17.276 16.276 20.075 13.264 11.610  3.000 
*  7.549 11.072 12.248 12.446  4.583  9.903 15.689  0.000 
*
* 
* $Id: tsrchd.f 2064 2009-02-24 06:11:28Z astivala $
*=======================================================================

      program tsrchd

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
      integer qn,dbn,info
      character dbfile*1024,qid*8,dbid*8
      logical ltype,lorder,lsoln
      double precision score
      integer i,j
*     ..
*     .. Local Arrays ..
      character*2 qtab(maxdim,maxdim),dbtab(maxdim,maxdim)
      double precision qdmat(maxdim,maxdim),dbdmat(maxdim,maxdim)
      double precision soln(maxdim*maxdim+maxdim+maxdim)
*     ..
*     .. Intrinsic Functions ..
!      intrinsic Time,CTime
*     ..
*     .. External Subroutines and Functions ..
      external tmatd,rdtabd,rdistm

*     ..
*     .. Executable Statements ..
*
      read (*,'(a1024)') dbfile
      read (*,'(l1,1x,l1,1x,l1)') ltype,lorder,lsoln
      call rdtabd(inunit, qtab, maxdim, qn, qid)
      call rdistm(inunit, qdmat, maxdim, qn)
      
      write(*, 5) ltype,lorder,lsoln
 5    format('# TSRCHD LTYPE = ',l1,' LORDER = ',l1, ' LSOLN = ',l1)
      write(*, 6) qid
 6    format('# QUERY ID = ',a8)
      write(*, 7) dbfile
 7    format('# DBFILE = ',a80)
!      write(*,8) CTime(Time())
! 8    format('# ',a24)

      open(dbunit, file=dbfile, status='OLD')
 10   if (.true.) then
         call rdtabd(dbunit, dbtab, maxdim, dbn, dbid)
         call rdistm(dbunit, dbdmat, maxdim, dbn)
         call tmatd(qtab, maxdim, qn, dbtab, maxdim, dbn, 
     $        ltype, lorder, qdmat, dbdmat,
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
      
