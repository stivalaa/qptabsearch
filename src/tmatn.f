*=======================================================================
* File:    tmatn.f
* Author:  Alex Stivala
* Created: July 2008
*
* FORTRAN-77 implementation of tableau matching numeric (ie Omega
* matrix rather than actual tableau) using the SPSOQP implementation
* of Prof. Yinyu Ye's interior point QP solver.
*
* See detailed documentation in tmatd.f header comments.
*
* $Id: tmatn.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================

      subroutine tmatn(omega1,ldo1,n1,omega2,ldo2,n2,
     $     ltype,lorder,dmat1,dmat2,score,soln,info)

      implicit none
*
* tmatn - find matching score of two Omega matrices using relaxed QP
*
*   The tableaux are specified here as double precision Omega matrices.
*

*
*     .. Scalar Arguments ..
      integer ldo1,n1,ldo2,n2,info
      logical ltype,lorder
      double precision score
*     ..
*     .. Array Arguments ..
      double precision   omega1(ldo1,*),omega2(ldo2,*)
      double precision   dmat1(ldo1,*),dmat2(ldo2,*)
      double precision   soln(*)
*     ..
*

*
* Arguments
* =========
*
* omega1 (input) DOUBLE PRECISION array, dimension (n1,n1)
*        Omega matrix for one structure. Symmetric.
*
* ldo1   (input) INTEGER
*        Leading dimension of omega1 array. ldo1 >= n1
*
* n1     (input) INTEGER
*        Dimension of omega1 matrix.
*
* omega2 (input) DOUBLE PRECISION array, dimension (n2,n2)
*        Omega matrix for second structure. Symmetric.
*
* ldo2   (input) INTEGER
*        Leading dimension of omega2 array. ldo2 >= n2
*
* n2     (input) INTEGER
*        Dimension of omega2 matrix.
*
* ltype  (input) LOGICAL
*        if true, disalow matches between SSEs that are not the same 
*        type. This uses the entry on the main diagonal as the indicator
*        of SSE type, and penalizes matches between entries with
*        different types.
*
* lorder (input) LOGICAL
*        if true, penalize matches between SSEs not maintaining sequence
*        order between the tableaux i.e. if i < k and j >= l for i,k
*        indices in omega1 and j,l indices in omega2.
*
* dmat1  (input) DOUBLE PRECISION array, dimension (n1,n1)
*        SSE distance matrix for one structure. symmetric
*
* dmat2  (input) DOUBLE PRECISION array, dimension (n2,n2)
*        SSE distance matrix for second structure. symmetric
*
* score  (output) DOUBLE PRECISION
*        score of matching omega1 with omega2 using relaxed QP.
*
* soln   (output) DOUBLE PRECISION vector, dimension (n1*n2+n1+n2)
*        solution x vector of matching omega1 with omega2 using relaxed QP
*
* info   (output) INTEGER
*        on exit, status of the computation
*        =  0 : successful exit
*        = -1 : maximum order of Q matrix exceeded
*        = -2 : maximum number of constraints exceeded
*        = -3 : SOLVQP failed
*
*=======================================================================
*

*     .. Parameters ..

*     largest tableau (omega matrix) dimension allowed
      integer maxdim1,maxdim2
      parameter(maxdim1 = 110, maxdim2 = 110)
      integer maxbdim, maxcdim
      parameter (maxcdim = maxdim1*maxdim2+maxdim1+maxdim2)
      parameter (maxbdim = maxdim1+maxdim2)
*     Threshold for penalizing difference in distances between SSEs
*     (Angstroms) - if difference in distances between SSEs in the two
*     structures exceeds this then penalize that match
      double precision mxssed
      parameter (mxssed = 4.0d0)
*     ..
*     .. Local Scalars ..
      integer i,j,k,l,qdim,arows
*     ..
*     .. Local Arrays ..
*     Q: symmetric objective matrix
*     A: constraint left-hand matrix
*     b: constraint rhs column vector
*     c: objective column vector
*     x: solution vector, with objective value as last element
*     y: dual solution vector (Langrangian multiplier)
      double precision Q(maxcdim, maxcdim),
     $     A(maxbdim, maxcdim),
     $     b(maxbdim),
     $     c(maxcdim),
     $     x(maxcdim + 1),
     $     y(maxbdim)
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. External Subroutines and Functions ..
      external dcopy
      external solvqp
      double precision tscorn
      external tscorn

*     ..
*     .. Executable Statements ..
*     
      qdim = n1*n2 + n1 + n2
      arows = n1 + n2
      if (qdim .gt. maxcdim) then
         write(*,*) 'Max order of Q matrix exceeded'
         info = -1
         return
      endif
      if (arows .gt. maxbdim) then
         write(*,*) 'Max number of constraints exceeded'
         info = -2
         return
      endif

*     zero Q and A and c. c will stay zero.
      do 20 j = 1, qdim
         c(j) = 0.0d0
         do 10 i = 1, qdim
            Q(i, j) = 0.0d0
            if (i .le. arows) then
               A(i, j) = 0.0d0
            endif
 10      continue
 20   continue

*     set b to ones
      do 30 j = 1, arows
         b(j) = 1.0d0
 30   continue

*     setup the constraint matrix A
      do 50 j = 1, n2
         do 40 i = 1, n1
            A(i, (i-1)*n2+j) = 1.0d0
 40      continue
 50   continue
      do 70 j = 1, n1
         do 60 i = 1, n2
            A(n1+i, (j-1)*n2+i) = 1.0d0
 60      continue
 70   continue
      do 80 i = 1, n1+n2
         A(i, n1*n2+i) = 1.0d0
 80   continue

C       write(*,99),( (A(i,j), j=1,maxcdim), i=1,maxbdim)
C  99   format (111(f8.4))
         
*     setup the objective matrix Q
*     NB values in Q negated since we are 
*     minimizing not maximizing (c is zero so irrelevant there)

      do 200 i = 1, n1
         do 190 j = 1, n2
            do 180 k = 1, n1
               do 170 l = 1, n2
*                 diagonal enries are SSE type not angle
                  if (i .ne. k .and. j .ne. l) then
                     if (ltype .and. 
     $                   (omega1(i,i) .ne. omega2(j,j) .or.
     $                    omega1(k,k) .ne. omega2(l,l))) then
*                       penalize matches between SSEs of different type
                        Q((i-1)*n2 + j, (k-1)*n2 + l) = 0.0d0
                     elseif (abs(dmat1(i,k) - dmat2(j,l)) .gt. mxssed)
*                       penalize matches where difference between SSE 
*                       distances exceeds the threshold
     $               then
                        Q((i-1)*n2 + j, (k-1)*n2 + l) = 0.0d0
                     else
                        if (lorder .and.
     $                       ( (i .lt. k .and. j .gt. l) .or.
     $                         (i .gt. k .and. j .lt. l) )) then
*                       penalize SSEs that are out of order
                           Q((i-1)*n2 + j, (k-1)*n2 + l) = 1.0d0
                        else
                           Q((i-1)*n2 + j, (k-1)*n2 + l) = -1.0d0 *
     $                          tscorn(omega1(i,k), omega2(j,l))
                        endif
                     endif
                  endif
 170           continue
 180        continue
 190     continue
 200  continue


C       write(*,205),( (Q(i,j), j=1,maxcdim), i=1,maxcdim)
C  205  format (111(f8.4))

      call solvqp(Q,maxcdim,A,maxbdim,arows,qdim,b,c,x,y,info)
      if (info .ne. 0) then
         write(*,210) info
 210     format('SOLVQP failed, info = ',i4)
         info = -3
         return
      endif

      score = x(qdim+1)
      call dcopy(qdim, x, 1, soln, 1)
      info = 0
      return
      end

