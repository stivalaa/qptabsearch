*=======================================================================
* File:    irtsrch.f
* Author:  Alex Stivala
* Created: March 2010
*
* FORTRAN-77 implementation of IR Tableau, tableau searching by
* vector cosine similarity as described by
*
* Zhang, Bailey, Konagurthu, Ramamohanarao 2010 "A fast indexing approach
* to protein structure comparison" BMC Bioinformatics 11(Suppl 1):S46
* 8th Asia-Pacific Bioinformatics Conference (APBC 2010)
*
*
* Usage: irtsrch < inputfile
*
* The 'database' to search is an ASCII file where each line contains
* the structure identifier and the 32-element integer vector representing
* the tableau for that structure reduced to the "bag of words" vector
* as described in the paper cited above.
* This is built by the convdbpacked2vector.py script.
*
*
* The results are printed to stdout as 
*
* name score
*
* This file is a program not a subroutine.
*
* Both the name of the database file to read, and the actual
* query tableau vector are read from stdin. 
* The first line is the name
* of the database file.
* The next line is the query vector in the same format as an entry
* in the database, i.e. the identifier then the 32-element vector.
*
* E.g.:
*
* /home/alexs/tableauxdb/ASTRAL-sel-gs-bib-95-1.75/tableaux.vectordb.ascii
* d1rroa_ 1 6 5 9 2 2 1 10 0 1 4 2 1 0 3 2 0 0 1 0 1 0 2 1 0 0 0 0 1 0 0 0
*
*
* Note the vectors are all integer but read as double precision just
* so we can directly call BLAS functions DDOT etc. with them.
* 
* $Id: irtsrch.F 3240 2010-01-18 03:28:54Z alexs $
*=======================================================================

      program irtsrch

      implicit none

*     .. Parameters ..
*     the vectors are all 32-element
      integer vecdim
      parameter(vecdim = 32)
*     unit number of standard input
      integer inunit
      parameter(inunit = 5)
*     unit number for database file
      integer dbunit
      parameter(dbunit=10)

*     ..
*     .. Local Scalars ..
      integer info
      character dbfile*1024,qid*8,dbid*8
      double precision score,qnorm
      integer i
*     ..
*     .. Local Arrays ..
      double precision qvec(vecdim),dbvec(vecdim)
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. External Subroutines and Functions ..
*     BLAS level 1 subroutines/functions
      double precision ddot,dnrm2
      external ddot,dnrm2
*     ..
*     .. Executable Statements ..
*
      read (*,'(a1024)') dbfile
      read (*, *) qid, (qvec(i), i = 1, vecdim)
      write(*, 5)
 5    format('# IRTSRCH')
      write(*, 6) qid
 6    format('# QUERY ID = ',a8)
      write(*, 7) dbfile
 7    format('# DBFILE = ',a80)
      qnorm = dnrm2(vecdim, qvec, 1)
      open(dbunit, file=dbfile, status='OLD')
 10   if (.true.) then
         read(dbunit,*,err=999,end=999) dbid, (dbvec(i), i = 1, vecdim)
*        compute score as cosine similarity between the two vectors
         score = ddot(vecdim, qvec, 1, dbvec, 1) / 
     $        (qnorm * dnrm2(vecdim, dbvec,1))
         
         write(*,30) dbid,score
 30      format(a8,f12.10)
         goto 10
      endif
 999  continue
      close(dbunit)
      end
      
