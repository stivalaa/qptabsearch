*=======================================================================
* File:    spphs2.f
* Author:  Alex Stivala, based on spphase1.m by Prof. Yinyu Ye.
* Created: July 2008
*
* This is a reimplementation using FORTRAN-77 and UMFPACK
* of spphase1.m from Prof. Yinyu Ye
* http://www.stanford.edu/~yyye/matlab/spphase2.m
*
* This uses sparse matrix routines in UMFPACK (v5.2.0)
* http://www.cise.ufl.edu/research/sparse/umfpack/
*
* $Id: spphs2.f 2072 2009-02-26 23:55:28Z astivala $
*=======================================================================

      subroutine spphs2(Q, ldq, A, lda, m, n, c, alpha, beta, lambda,
     $     toler, obvalue, ob,
     $     x, y, Ep, Ei, Ex, nzmax, u, h, xxvec, gg, XX, ldxx, 
     $     AA, ldaa,  wvec, symbolic, info2)

      implicit none
*
* spphs2 - phase 2 of dsoqp:
*          Repeatedly solve an ellipsoid constrained QP problem
*          by solving a linear system until a positive solution is
*          found.
*
*     .. Scalar Arguments ..
      integer            ldq,lda,m,n,nzmax,ldxx,ldaa,info2
      double precision   alpha,beta,lambda,toler,obvalue,ob
      integer*8          symbolic
*     ..
*     .. Array Arguments ..
      integer            Ep(*),Ei(*)
      double precision   Q(ldq,*),A(lda,*),c(*),x(*),y(*),
     $     Ex(*),u(*),h(*),xxvec(*),gg(*),XX(ldxx,*),
     $     AA(ldaa, *), wvec(*)

*     ..
*
* The phase 1 procedure called by spsoqp.
* For more details see spsoqp.f.
*
* Arguments
* =========
*
* Q      (input) DOUBLE PRECISION array, dimension (n, n)
*        The symmetric objective matrix
*
* 
* ldq    (input) INTGER
*        leading dimension of Q. ldq >= n
*
* A      (input) DOUBLE PRECISION array, dimension (m, n)
*        the constraint left-hand matrix
* 
* lda    (input) INTGER
*        leading dimension of A. lda >= m
*
* m      (input) INTEGER
*        number of rows in matrix A
*
* n      (input) INTEGER
*        number of columns in matrix A
*
* c      (input) DOUBLE PRECISION vector, dimension(n)
*        objective column vector
*
* alpha  (input) DOUBLE PRECISION
*        solution descent parameter
*
* beta   (input) DOUBLE PRECISION
*        step size 0 < beta < 1
*
* lambda (input/output) DOUBLE PRECISION
*        solution descent parameter
*
* toler  (input) DOUBLE PRECISION
*        relative stopping tolerance
*
* obvalue (input) DOUBLE PRECISION
*        previous objective value
*
* ob     (input/output) DOUBLE PRECISION
*        objective value
*
* x      (input/output) DOUBLE PRECISION vector, dimension (n)
*        solution vector
*
* y      (output) DOUBLE PRECISION vector, dimension (m)
*        optimal dual solution (Lagrangian multiplier)
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
* u      (output) DOUBLE PRECISION vector, dimension (m+n)
*        workspace
*        solution from UMFPACK
*
* h      (output) DOUBLE PRECISION vector, dimension (m+n)
*        workspace
*        rhs vector UMFPACK
*
* xxvec  (output) DOUBLE PRECISION vector, dimension (n)
*        workspace
*
* gg     (output) DOUBLE PRECISION vector, dimension (n)
*        workspace
*
* XX     (output) DOUBLE PRECISION, dimension (n, n)
*        workspace
* 
* ldxx   (input) INTEGER
*        leading dimension of XX. ldxx >= n
*
* AA     (output) DOUBLE PRECISION, dimension (m, n)
*        workspace
*
* ldaa   (input) INTEGER
*        leading dimension of AA. ldaa >= m.
*
* wvec    (output) DOUBLE PRECISION, dimension (n)
*        workspace     
*
* symbolic (input/output) INTEGER*8
*        handle for UMFPACK symbolic pre-order and analsysis
*
* info2  (output) INTEGER
*        = 0  : successful
*        = 1  : too many nonzero values, increase nzmax
*        = 999: the problem seems unbounded
*       other : info code from LAPACK
*
*=======================================================================
*

*     .. Parameters ..
*     eps is the machine epsilon for double precision
*     this is for GNU FORTRAN (g77) on Linux x86/64
*     MAY NEED CHANGING ON DIFFERENT SYSTEMS/COMPILERS
      double precision eps
      parameter (eps = 2.2204d-16)
      double precision infty
      parameter (infty=1.79769d+308)
      double precision neginf
      parameter (neginf=-1.79769d+308)
*     ..
*     .. Local Scalars ..
      integer i,j,incr,order,ixE
      double precision w1,w2,blasa,blasb,nora,valgo
      character uplo
      integer*8 numeric,sys
*     ..
*     .. Local Arrays ..
      double precision control(20),infou(90)
*     ..
*     .. Intrinsic Functions ..
      intrinsic min,abs,int
*     ..
*     .. External Subroutines and Functions ..
      external demvv
      double precision dvmin
      external dvmin
*     UMFPACK subroutines
      external umf4def,umf4sym,umf4num,umf4solr,umf4fsym,umf4fnum
*     BLAS level 1 subroutines/functions
      external daxpy,dscal,dcopy
      double precision ddot
      external ddot
*     BLAS level 2 subroutines
      external dsymv
*     BLAS level 3 subroutines

*     ..
*     .. Executable Statements ..
*

*     ==================================================================
*
*     Repeatedly solve the ellipsoid constrained QP problem by solving
*     the linear system
*
*     E*u = h where
*
*     E (m+n, m+n)  =    ( XXlam (n, n)   AA' (n, m) )
*                        ( AA    (m, n)    0  (m, m) )
*
*     h (m+n)       =    ( dx (n, 1) )
*                        ( 0  (m, 1) )
*
*     dx (n, 1) is elementwise -x.*gg where gg = Q*x+c
*     XXlam = XX + lambda*eye(n,n) where XX = dx * Q * dx
*                                  where dx(n) is x along diagonal
*
*     until a positive solution is found.
*
*     XX (and therefore XXlam) is square, symmetric. But not really
*     sparse; it has a kind of 'grid' structure but may be e.g. 58%
*     nonzero.  
*
*     E is square, symmetric, and sparse (but not very sparse
*     as in the other cases in DPSOQP,SPPHS1, since the top left
*     submatrix is XX which is not sparse). E.g. may be 43% nonzero.
*
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
*     Inside the loop, only the top left submatrix of E is updated
*     (lambda is modified), other 3 submatrices are constant.
*     So we can build the matrix in sparse (Harwell-Boeing) format
*     outside the loop, and inside the loop need only update
*     the n diagonal elements in top left n by n submatrix by
*     adding new value of lambda to diagonal elements of XX.
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

      lambda = (1.0d0 - beta) * lambda

      valgo = 0.0d0

*     always using increment of 1 in BLAS calls
      incr = 1

*     gg = Q*x + c
*     Q is symmetric  (using upper tri arbitrarily here, see TODO below)
*     stored in gg not directly in F(1:n) as used again later 
      call dcopy(n, c, incr, gg, incr)
      uplo = 'U'
      blasa = 1.0d0
      blasb = 1.0d0
      call dsymv(uplo, n, blasa, Q, ldq, x, incr, blasb, gg, incr)


*     dx = diagonal matrix with x along main diagonal
*     AA = A*dx
*     XX = dx*Q*dx
*     since XX is diagonal we can compute 
*     XX <- A*dx and 
*     XX <- Q*dx just by
*     scaling each column by x(j) for j=1,n
*     then finally we can compute
*     XX <- dx*XX
*     by scaling each ROW by x(i) for i = 1,n
*     with only BLAS level 1 subroutines DCOPY and DSCAL,
*     although the last one is row-wise so inefficient memory access
*     in FORTRAN.
      do 10 j = 1,n
         call dcopy(m, A(1, j), incr, AA(1, j), incr)
         call dscal(m, x(j), AA(1, j), incr)
         call dcopy(n, Q(1, j), incr, XX(1, j), incr)
         call dscal(n, x(j), XX(1, j), incr)
 10   continue
      do 20 i = 1,n
*        NB using INCX=ldxx to scale rows not columns of XX
         call dscal(n, x(i), XX(i, 1), ldxx)
 20   continue
             

*     build constant part of rhs vector h
*     h(n+1:n+m) = 0
      do 190 i = n+1, n+m
         h(i) = 0.0d0
 190  continue


 200  if (valgo .le. 0.0d0) then

*     (For now we are building the whole sparse E every iteration in
*     the loop)
*        TODO:
*        update E by seting diaongal elements of top left n by n submatrix
*        to diagonal of XX + lambda. (We previously made sure these
*        elements are present in the sparse storage format, even if zero).
*     This doesn't matter much as this loop typically only has one
*     iteration anyway - the expense is multiple calls of SPPHS2
*     from SPSOQP.

*     ixE is the current index (FORTRAN 1-based) in Ei and Ex.
*     We will subtract one from the values stored in the arrays
*     for UMFAPCK zero-based convention.
*     We make sure we always store the diagonal elements of XX in 
*     top left submatrix of E, even if zero, since we update them in loop.
         ixE = 1
         Ep(1) = 0
         do 250 j = 1, n
            do 240 i = 1, n
               if (i .eq. j) then 
                  Ei(ixE) = i - 1
                  Ex(ixE) = XX(i, j) + lambda
                  ixE = ixE + 1
                  if (ixE .gt. nzmax) then
                     write(*,*) 
     $                    'Max number of nonzero values exceeded.'
                     info2 = 1
                     return
                  endif
               elseif  (XX(i, j) .ne. 0.0d0) then
                  Ei(ixE) = i - 1
                  Ex(ixE) = XX(i, j)
                  ixE = ixE + 1
                  if (ixE .gt. nzmax) then
                     write(*,*) 
     $                    'Max number of nonzero values exceeded.'
                     info2 = 1
                     return
                  endif
               endif
 240        continue
            do 245, i = 1, m
               if (AA(i,j) .ne. 0.0d0) then
                  Ei(ixE) = n+i - 1
                  Ex(ixE) = AA(i,j)
                  ixE = ixE + 1
                  if (ixE .gt. nzmax) then
                     write(*,*) 
     $                    'Max number of nonzero values exceeded.'
                     info2 = 1
                     return
                  endif
               endif
 245        continue
            Ep(j+1) = ixE - 1
 250     continue
         do 270 j = 1, m
            do 260 i = 1, n
               if (AA(j,i) .ne. 0.0d0) then
                  Ei(ixE) = i - 1
                  Ex(ixE) = AA(j,i)
                  ixE = ixE + 1
                  if (ixE .gt. nzmax) then
                     write(*,*)'Max number of nonzero values exceeded.'
                     info2 = 1
                     return
                  endif
               endif
 260        continue
            Ep(n+j+1) = ixE - 1
 270     continue
         
         
*        update nonconstant part of rhs vector h
*        h(1:n) = -x .* gg
         call dcopy(n, x, incr, h, incr)
         call demvv(n, gg, h)
         blasa = -1.0d0
         call dscal(n, blasa, h, incr)

         order = m + n

*        pre-order and symbolic analysis
*        We do this once on the first iteration of loop in first call
*        of this function, it is then reused on subsequent 
*        iterations + calls as the structure of the matrix does not change
*        and analysis is an expensive operation
         if (symbolic .eq. 0) then
            call umf4sym(order,order,Ep,Ei,Ex,symbolic,control,infou)
            if (infou(1) .lt. 0) then
               write(*,*) 'SPPHS2 UMF4SYM failed, info = ',infou(1)
               info2 = int(infou(1))
               return
            endif
         endif
*        numeric factorization
         call umf4num(Ep, Ei, Ex, symbolic, numeric, control, infou)
         if (infou(1) .lt. 0) then
            write(*,*) 'SPPHS2 UMF4NUM failed, info = ', infou(1)
            info2 = int(infou(1))
            return
         endif
*        solve E*u=h with iterative refinement
         sys = 0
         call umf4solr(sys, Ep, Ei, Ex, u, h, numeric, control, infou)
         call umf4fnum(numeric)
         if (infou(1) .lt. 0) then
            write(*,*) 'SPPHS2 UMF4SOLR failed, info = ', infou(1)
            info2  = int(infou(1))
            return
         endif

*        The solution is in u
         
*        xxvec = x + x .* u(1:n)
         call dcopy(n, x, incr, xxvec, incr)
         call demvv(n, u, xxvec)
         blasa = 1.0d0
         call daxpy(n, blasa, x, incr, xxvec, incr)

         valgo = dvmin(n, xxvec)
         if (valgo .gt. 0.0d0) then
*           ob = (xxvec' * Q * xxvec)/2 + c'*xxvec
            uplo = 'U'
            blasa = 1.0d0
            blasb = 0.0d0
*           blasb is zero so wvec can be uninitialized on entry to DSYMV
            call dsymv(uplo, n, blasa, Q, ldq, xxvec, incr, blasb, 
     $                 wvec, incr)
            ob = ddot(n, xxvec, incr, wvec, incr) / 2.0d0 + 
     $           ddot(n, c, incr, xxvec, incr)

C            write(*,399),valgo,obvalue,ob
C 399        format('valgo = ',f10.6,' obvalue = ',f10.6,
C     $           ' ob = ',f10.6)

            valgo = min(valgo, obvalue - ob + eps)
         endif
         lambda = 2.0d0 * lambda

C          write(*,400),lambda
C  400     format ('lambda = ',f12.6)

         if (lambda .ge. (1.0d0 + abs(obvalue))/toler) then
            write(*,*) 'The problem seems unbounded.'
            info2 = 999
            return
         endif
         goto 200
*     end of while (valgo .le .0.0) loop
      endif

*     y = -u(n+1:n+m)
      call dcopy(m, u(n+1), incr, y, incr)
      blasa = -1.0d0
      call dscal(m, blasa, y, incr)
*     from now on we only use first n elements of u
      nora = dvmin(n, u)
      if (nora .lt. 0.0d0) then
         nora = -alpha / nora
      elseif (nora .eq. 0.0d0) then
         nora = alpha
      else
         nora = infty
      endif

*     u = x.* u
*     w1 = u'*Q*u
*     w2 = -u'*gg
      call demvv(n, x, u)
      uplo = 'U'
      blasa = 1.0d0
      blasb = 0.0d0
*     blasb is zero so wvec can be uninitialized on entry to DSYMV
      call dsymv(uplo, n, blasa, Q, ldq, u, incr, blasb, wvec, incr)
C      call dgemv('N',n, n, blasa, Q, ldq, u, incr, blasb, wvec, incr)
      w1 = ddot(n, u, incr, wvec, incr)
      w2 = -1.0d0 * ddot(n, u, incr, gg, incr)
      if (w1 .gt. 0.0d0) then
         nora = min(w2/w1, nora)
      endif
      if (nora .eq. infty) then
         ob = neginf
      else
*        x = x + nora*u
         call daxpy(n, nora, u, incr, x, incr)
*        ob = (x' * Q * x)/2 + c'*x
         uplo = 'U'
         blasa = 1.0d0
         blasb = 0.0d0
*        blasb is zero so wvec can be uninitialized on entry to DSYMV
C         call dsymv(uplo, n, blasa, Q, ldq, x, incr, blasb, wvec, incr)
         call dgemv('N',n, n, blasa, Q, ldq, x, incr, blasb, wvec, incr)
         ob = ddot(n, x, incr, wvec, incr) / 2.0d0 +
     $        ddot(n, c, incr, x, incr)
      endif

*     end of SPPHS2 , phase 2 subroutine of SPSOQP.

      info2 = 0
      return
      end

