*=======================================================================
* File:    spphs1.f
* Author:  Alex Stivala, based on spphase1.m by Prof. Yinyu Ye.
* Created: July 2008
*
* This is a reimplementation using FORTRAN-77 and UMFPACK
* of spphase1.m from Prof. Yinyu Ye
* http://www.stanford.edu/~yyye/matlab/spphase1.m
*
* This uses sparse matrix routines in UMFPACK (v5.2.0)
* http://www.cise.ufl.edu/research/sparse/umfpack/
*
* $Id: spphs1.f 2073 2009-02-27 00:03:39Z astivala $
*=======================================================================

      subroutine spphs1(A, lda, m, n, avec, b, x, y, z, ob, alpha,
     $     dx, Ep, Ei, Ex, nzmax, F, ldf, u, v, umzv, obz,
     $     y1, y2,
     $     symbolic,
     $     info)

      implicit none
*
* spphs1 - phase 1 of spsoqp: find a feasible point
*
*     .. Scalar Arguments ..
      integer            lda,m,n,nzmax,ldf
      double precision   z,ob,alpha
      integer*8          symbolic

*     ..
*     .. Array Arguments ..
      double precision   A(lda,*), x(*), y(*), avec(*), b(*)
      double precision   dx(*),Ex(*),F(ldf, *),
     $     u(*), v(*), umzv(*), obz(*), y1(*), y2(*)
      integer            Ep(*),Ei(*)

*     ..
*
* The phase 1 procedure called by spsoqp.
* For more details see spsoqp.f.
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
* Ep     (output) INTEGER vector, dimension >= m+n+1
*        workspace
*        sparse storage for matrix E for UMFPACK
*
* Ei     (output) INTEGER vector, dimensino (nzmax)
*        workspace
*        sparse storage for matrix E for UMFPACK
*
* Ex     (output) DOUBLE PRECICISION vector, dimension (nzmax)
*        workspace
*        sparse storage for matrix E for UMFPACK
*
* nzmax  (input) INTEGER
*        maximum allowed number of nonzero values in E. 
*        Dimensions Ei and Ex.
*
* F      (output) DOUBLE PRECISION, dimension (m+n, 2)
*        workspace
*
* ldf    (input) INTEGER
*        leading dimension of F. ldf >= m+n
*
* y1     (output) DOUBLE PRECISION vector, dimension(m+n)
*        workspace
*        solution from UMFPACK with rhs column 1 of F
*
* y2     (output) DOUBLE PRECISION vector, dimension(m+n)
*        workspace
*        solution from UMFPACK with rhs column 2 of F
*
* u      (output) DOUBLE PRECISION vector, dimension (n+2)
*        workspace
*
* v      (output) DOUBLE PRECISION vector, dimension (n+2)
*        workspace
*
* umzv   (output) DOUBLE PRECISION vector, dimension (n+2)
*        workspace
*
* obz    (output) DOUBLE PRECISION vector, dimension (n+2)
*        workspace
*
* symbolic (output) INTEGER*8
*        handle for UMFPACK symbolic pre-order and analsysis
*
* info   (output) INTEGER
*        = 0  : successful
*        = 1  : max number of nonzero values exceeded, increase nzmax
*       other  : info code from UMFPACK
*
*=======================================================================
*

*     .. Parameters ..

*     ..
*     .. Local Scalars ..
      integer i,j,order,info,incr,uvlen,xlen,ixE
      integer*8 numeric,sys
      double precision w1,w2,obrec,obrec2,adoty2,blasa,blasb,nora,vn2rec
      character trans
C      character uplo
*     ..
*     .. Local Arrays ..
      double precision control(20),infou(90)
*     ..
*     .. Intrinsic Functions ..
      intrinsic min,int
*     ..
*     .. External Subroutines and Functions ..
      external drecip,demvv
      double precision dvmin,dvmax
      external dvmin,dvmax
*     UMFPACK subroutines
      external umf4def,umf4sym,umf4num,umf4solr,umf4fsym,umf4fnum
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
*     But it is not necessarily
*     positive definite so we cannot always use CHOLMOD to efficiently
*     solve the system. We therefore use UMFPACK with sparse matrix
*     E stored in packed column (Harwell-Boeing) format.
*     In this format, the nonzero entries are stored (column major)
*     in Ex. For each column the row indices of the nonzero entries
*     are stored in Ei. For each column k, Ep(k) and Ep(k+1)-1 are
*     respecitively the first and last indices of row numbers in Ei
*     and values in Ex for that column.
*     NB UMFPACK uses 0-based indexing in Ep and Ei.
*     No duplicate row indices are allowed, and row indicies for each
*     column must be in ascending order.
*
*     If E were positive definite, which can efficiently be tested
*     for with Cholesky factorization, it could be solved by that means
*     more efficiently than the general case using LU decomposition.
*
*     TODO: use CHOLMOD to test for postive definiteness and use if
*           it is.
*
*     ==================================================================


*     set default parameters for UMFPACK
      call umf4def(control)

      incr = 1
      call dcopy(n, x, incr, dx, incr)
*     seems to be no BLAS subroutine for elementwise reciprocal, use own
      call drecip(n, dx)

*     ixE is the current index (FORTRAN 1-based) in Ei and Ex.
*     We will subtract one from the values stored in the arrays
*     for UMFAPCK zero-based convention.
      ixE = 1
      Ep(1) = 0
      do 30 j = 1, n
         Ei(ixE) = j - 1
         Ex(ixE) = dx(j)**2
         ixE = ixE + 1
         do 20 i = 1, m
            if (A(i,j) .ne. 0.0d0) then
               Ei(ixE) = n+i - 1
               Ex(ixE) = A(i, j)
               ixE = ixE + 1
               if (ixE .gt. nzmax) then
                  write(*,*) 'Max numbero of nonzero values exceeded.'
                  info = 1
                  return
               endif
            endif
 20      continue
         Ep(j+1) = ixE - 1
 30   continue
      do 50 j = 1, m
         do 40 i = 1, n
            if (A(j,i) .ne. 0.0d0) then
               Ei(ixE) = i - 1
               Ex(ixE) = A(j,i)
               ixE = ixE + 1
               if (ixE .gt. nzmax) then
                  write(*,*) 'Max number of nonzero values exeeded.'
                  info = 1
                  return
               endif
            endif
 40      continue
         Ep(n+j+1) = ixE - 1
 50   continue

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


      order = m + n

*     pre-order and symbolic analysis
      call umf4sym(order, order, Ep, Ei, Ex, symbolic, control, infou)
      if (infou(1) .lt. 0) then
         write(*,*) 'SPPHS1 UMF4SYM failed, info = ',infou(1)
         info = int(infou(1))
         return
      endif
*     numeric factorization
      call umf4num(Ep, Ei, Ex, symbolic, numeric, control, infou)
*     NB we do not free the symbolic handle as we reuse the symbolic
*     factorization done here back in the caller spsoqp, where
*     the matrix has the same structure.
      if (infou(1) .lt. 0) then
         write(*,*) 'SPPHS1 UMF4NUM failed, info = ', infou(1)
         info = int(infou(1))
         return
      endif
*     solve E*y1=F(:,1) with iterative refinement
      sys = 0
      call umf4solr(sys, Ep, Ei, Ex, y1, F, numeric, control, infou)
      if (infou(1) .lt. 0) then
         write(*,*) 'SPPHS1 UMF4SOLR [1] failed, info = ', infou(1)
         info  = int(infou(1))
         return
      endif
*     solve E*y2=F(:,2) with iterative refinement
      sys = 0
      call umf4solr(sys, Ep, Ei, Ex, y2, F(1,2), numeric, control,
     $     infou)
      call umf4fnum(numeric)
      if (infou(1) .lt. 0) then
         write(*,*) 'SPPHS1 UMF4SOLR [2] failed, info = ', infou(1)
         info  = int(infou(1))
         return
      endif


*     We will use the solution vectors y1, y2 in place,
*     using just the m-vectors starting at row n+1 (i.e. the vectors
*     of dimension m at y1(n+1) and y2(n+1))
      obrec = 1.0d0 / ob
      obrec2 = 1.0d0 / ob**2
      adoty2 = ddot(m, avec, incr, y2(n+1), incr)
      w1 = (obrec - ddot(m, avec, incr, y1(n+1), incr))  /
     $     (obrec2 - adoty2)
      w2 = 1.0d0 / (obrec2 - adoty2)
*     y1 = y1 - w1*y2
      blasa = -1.0d0 * w1
      call daxpy(m, blasa, y2(n+1), incr, y1(n+1), incr)
*     y2 = -w2*y2
      blasa = -1.0d0 * w2
      call dscal(m, blasa, y2(n+1), incr)
      

      w1 = ddot(m, b, incr, y1(n+1), incr)
      w2 = ddot(m, b, incr, y2(n+1), incr)
      blasa = 1.0d0 / (1.0d0 + w1)
      call dscal(m, blasa, y1(n+1), incr)
*     y2 = y2 - w2*y1
      blasa = -1.0d0 * w2
      call daxpy(m, blasa, y1(n+1), incr, y2(n+1), incr)
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
      call dgemv(trans,m,n,blasa,A,lda,y2(n+1),incr,blasb,u,incr)
      call demvv(n, x, u)
      u(n+1) = x(n+1) * (1.0d0 - ddot(m, y2(n+1), incr, avec, incr))
      u(n+2) = w2/(1+w1)
      
      trans = 'T'
      blasa = 1
      blasb = 0
      call dgemv(trans,m,n,blasa,A,lda,y1(n+1),incr,blasb,v,incr)
      call demvv(n, x, v)
      v(n+1) = x(n+1) * ddot(m, y1(n+1), incr, avec, incr)
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
         call dcopy(m, y2(n+1), incr, y, incr)
         call daxpy(m, z, y1(n+1), incr, y, incr)
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

      info = 0
      return

      end

