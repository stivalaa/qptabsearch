*=======================================================================
* File:    d1phs1.f
* Author:  Alex Stivala, based on spphase1.m by Prof. Yinyu Ye.
* Created: July 2008
*
* This is a reimplementation using FORTRAN-77 and LAPACK
* of spphase1.m from Prof. Yinyu Ye
* http://www.stanford.edu/~yyye/matlab/spphase1.m
*
* Note that the original MATLAB code uses MATLAB sparse matrix
* facilities, implemented in MATALB with the CHOLMOD routines
* from Tim Davis et al. LAPACK is for dense systems, this implementation
* uses just LAPACK so does not use sparseness (hence no 'sp' prefix
* on name, 'dl' is for 'dense LAPACK').
*
* This appears to be faster than the sparse version (spphs1.f) which
* uses UMFPACK though.
*
* $Id: dlphs1.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================

      subroutine dlphs1(A, lda, m, n, avec, b, x, y, z, ob, alpha,
     $     dx, E, lde, F, ldf, ipiv, u, v, umzv, obz, 
     $     work, lwork, info1)

      implicit none
*
* dlphs1 - phase 1 of dlsoqp: find a feasible point
*
*     .. Scalar Arguments ..
      integer            lda,m,n,lde,ldf,lwork,info1
      double precision   z,ob,alpha
*     ..
*     .. Array Arguments ..
      double precision   A(lda,*), x(*), y(*), avec(*), b(*)
      double precision   dx(*),E(lde, *),F(ldf, *),
     $     u(*), v(*), umzv(*), obz(*), work(*)
      integer            ipiv(*)

*     ..
*
* The phase 1 procedure called by dlsoqp.
* For more details see dlsoqp.f.
*
* Arguments
* =========
*
* A      (input) DOUBLE PRECISION array, dimension (m, n)
*        the constraint left-hand matrix
*
* lda    (input) INTEGER
*        leading dimenion of A. lda >= m
*
* m      (input) INTEGER
*        number of rows in matrix A
*
* n      (input) INTEGER
*        number of columns in matrix A
*
* avec   (input) DOUBLE PRECISION vector, dimension (m)
*        b - A1_{n\times1}
*
* b      (input) DOUBLE PRECISION vector, dimension(m)
*        constraint right-hand column vector
*
* x      (input/output) DOUBLE PRECISION vector, dimension (n+1)
*        solution vector with ob value appended as last element
*
* y      (output) DOUBLE PRECISION vector, dimension (m)
*        optimal dual solution (Lagrangian multiplier)
*
* z      (input/output) DOUBLE PRECISION
*        objective lower bound
*
* ob     (input/output) DOUBLE PRECISION
*        objective value
*
* alpha  (input/output) DOUBLE PRECISION
*        solution descent parameter
*
* dx     (output) DOUBLE PRECISION, dimension (n)
*        workspace
*
* E      (output) DOUBLE PRECISION, dimension (m+n, m+n)
*        workspace
*
* lde    (input) INTEGER
*        leading dimension of E. lde >= m+n
*
* F      (output) DOUBLE PRECISION, dimension (m+n, 2)
*        workspace
*
* ldf    (input) INTEGER
*        leading dimension of F. ldf >= m+n
*
* ipiv   (output) INTEGER, dimension (m+n)
*        workspace
*
* u      (output) DOUBLE PRCISION vector, dimension (n+2)
*        workspace
*
* v      (output) DOUBLE PRCISION vector, dimension (n+2)
*        workspace
*
* umzv   (output) DOUBLE PRECISION vector, dimension (n+2)
*        workspace
*
* obz    (output) DOUBLE PRECISION vector, dimension (n+2)
*        workspace
*
* work   (output) DOUBLE PRECISION vector, dimension (lwork)
*        workspace for DSYSV
*
* lwork  (input) INTEGER
*        length of work vector. lwork >= (m+n)*NB where NB is optimal
*        blocksize for DSYTRF.
*
* info1  (output) INTEGER
*        = 0  : successful
*       <> 0  : info code from LAPACK
*
*=======================================================================
*

*     .. Parameters ..

*     ..
*     .. Local Scalars ..
      integer i,j,order,info,nrhs,incr,uvlen,xlen
      double precision w1,w2,obrec,obrec2,adoty2,blasa,blasb,nora,vn2rec
      character trans
      character uplo
*     ..
*     .. Local Arrays ..

*     ..
*     .. Intrinsic Functions ..
      intrinsic min
*     ..
*     .. External Subroutines and Functions ..
      external drecip,demvv
      double precision dvmin,dvmax
      external dvmin,dvmax
*     LAPACK subroutines
      external dlange,dpotrf,dpotrs,dpocon,dsysv
*     BLAS level 1 subroutines/functions
      external daxpy,dscal,dcopy
      double precision ddot
      external ddot
*     BLAS level 2 subroutines
      external dgemv
*     ..
*     .. Executable Statements ..
*
*     ==================================================================
*
*     solve the scaled least squares against two vectors.
*
*     E*yy = F where
*
*     E (m+n, m+n)  =    ( DD (n, n)     A' (n, m) )
*                        ( A  (m, n)     0  (m, m) )
*
*     F (m+n, 2)    =    ( dx (n, 1)     0 (n, 1) )
*                        ( 0  (m, 1)     a (m, 1) )
*
*     dx (n, 1) is elementwise 1/x
*     DD (n, n) has elementwise dx**2 on diagonal
*
*     E is square, symmetric, and sparse. 
*     It may be positive definite, which can efficiently be tested
*     for with Cholesky factorization, and solved by that means
*     more efficiently than the general case using LU decomposition.
*
*     ==================================================================

      incr = 1
      call dcopy(n, x, incr, dx, incr)
*     seems to be no BLAS subroutine for elementwise reciprocal, use own
      call drecip(n, dx)

*    symmetric matrix so store only upper triangle and use DSYSV
*    (tried using packed format and DPPSV but that is slower)

      do 30 j = 1, n
         do 20 i = 1, j
            if (i .eq. j) then
               E(i,j) = dx(i)**2
            else
               E(i,j) = 0.0d0
            endif
 20      continue
 30   continue
      do 60 j = 1, m
*        note using INCX=lda for A to get rows
         call dcopy(n, A(j, 1), lda, E(1,n+j), incr)
 60   continue
      do 80 j = n+1, n+m
         do 70 i = n+1, j
            E(i,j) = 0.0d0
 70      continue
 80   continue


      do 100 i = 1, n
         F(i, 1) = dx(i)
 100  continue
      do 110 i = n+1, n+m
         F(i, 1) = 0.0d0
 110  continue
      do 120 i = 1, n
         F(i, 2) = 0.0d0
 120  continue
      do 130 i = 1, m
         F(n+i, 2) = avec(i)
 130  continue

C       write(*,999),( (E(i,j), j=1,m+n), i=1,m+n)
C  999  format (21(f4.1))

C       uplo = 'U'
C       order = m + n
C       call dpotrf(uplo, order, E, order, info)
C       if (info .ne. 0) then
C          write(*,*),'dlphs1 E matrix is not positive definite'
C       endif
*     TODO: use DPOTRS to solve if positive definite

      order = m + n
      nrhs = 2
      uplo = 'U'
      call dsysv(uplo, order, nrhs, E, lde, ipiv, F, ldf, 
     $           work, lwork, info)
      if (info .ne. 0) then
         write(*,*) 'dlphs1 DSYSV failed'
         info1 = info
         return
      endif
*     The solution is in F. We will use column vectors of F in place,
*     using just the m-vectors starting at row n+1 (i.e. the vectors
*     of dimension m at F(n+1, j) for j=1,2) 
*     (referred to as y1, y2, respectively).
      obrec = 1.0d0 / ob
      obrec2 = 1.0d0 / ob**2
      adoty2 = ddot(m, avec, incr, F(n+1, 2), incr)
      w1 = (obrec - ddot(m, avec, incr, F(n+1, 1), incr))  /
     $     (obrec2 - adoty2)
      w2 = 1.0d0 / (obrec2 - adoty2)
*     y1 = y1 - w1*y2
      blasa = -1.0d0 * w1
      call daxpy(m, blasa, F(n+1, 2), incr, F(n+1, 1), incr)
*     y2 = -w2*y2
      blasa = -1.0d0 * w2
      call dscal(m, blasa, F(n+1, 2), incr)
      

      w1 = ddot(m, b, incr, F(n+1, 1), incr)
      w2 = ddot(m, b, incr, F(n+1, 2), incr)
      blasa = 1.0d0 / (1.0d0 + w1)
      call dscal(m, blasa, F(n+1, 1), incr)
*     y2 = y2 - w2*y1
      blasa = -1.0d0 * w2
      call daxpy(m, blasa, F(n+1, 1), incr, F(n+1, 2), incr)
*
*     build vector u dimension (n+2) where
*
*     u(1:n) = x(1:n).*(-y2'*A)' 
*                 note (-y2'*A)' = A'*(-y2) 
*                 computed as the latter with DGEMV
*     u(n+1) = x(n+1)*(1-y2'*avec)
*     u(n+2) = w2/(1+w1)
*
*     and vector v dimesion (m+2) where
*
*     v(1:n) = x(1:n).*(y1'*A)'
*                 note (y1'*A)' = A'*y1
*                 computed as the latter with DGEMV
*     v(n+1) = x(n+1)*(y1'*avec)
*     v(n+2) = 1/(1+w1)
*
      trans = 'T'
      blasa = -1
      blasb = 0
      call dgemv(trans,m,n,blasa,A,lda,F(n+1, 2),incr,blasb,u,incr)
      call demvv(n, x, u)
      u(n+1) = x(n+1) * (1.0d0 - ddot(m, F(n+1, 2), incr, avec, incr))
      u(n+2) = w2/(1+w1)
      
      trans = 'T'
      blasa = 1
      blasb = 0
      call dgemv(trans,m,n,blasa,A,lda,F(n+1, 1),incr,blasb,v,incr)
      call demvv(n, x, v)
      v(n+1) = x(n+1) * ddot(m, F(n+1, 1), incr, avec, incr)
      v(n+2) = 1/(1+w1)
      
*
*     update the dual and the objective lower bound
*
      
*     umzv  = u - z*v
      uvlen = n + 2
      call dcopy(uvlen, u, incr, umzv, incr)
      blasa = -1.0d0 * z
      call daxpy(uvlen, blasa, v, incr, umzv, incr)
      if (dvmin(uvlen, umzv) .ge. 0.0d0) then
*        y = y2 + z*y1 
         call dcopy(m, F(n+1, 2), incr, y, incr)
         call daxpy(m, z, F(n+1, 1), incr, y, incr)
         z = ddot(m, b, incr, y, incr)
      endif

*     
*     find the descent direction
*

      do 200, i = 1, n+2
         obz(i) = (ob - z) / (n+2)
 200  continue
*     u = u - z*v - obz
      blasa = -1.0d0 * z
      call daxpy(uvlen, blasa, v, incr, u, incr)
      blasa = -1.0d0
      call daxpy(uvlen, blasa, obz, incr, u, incr)
      nora = dvmax(uvlen, u)

*     
*     update the solution along the descent direction
*
      
      if (nora .eq. u(n+1)) then
         alpha = 1.0d0
      end if
*     v = ones(n+2,1) - (alpha/nora)*u
      do 300 i = 1, n+2
         v(i) = 1.0d0
 300  continue
      blasa = -1.0d0 * (alpha/nora)
      call daxpy(uvlen, blasa, u, incr, v, incr)
*     x = ( x .* v(1:n+1) ) / v(n+2)
      xlen = n + 1
      call demvv(xlen, v, x)
      vn2rec = 1.0d0 / v(n+2)
      call dscal(xlen, vn2rec, x, incr)

      info1 = 0
      return

      end

