*=======================================================================
* File:    dlphs2.f
* Author:  Alex Stivala, based on spphase1.m by Prof. Yinyu Ye.
* Created: July 2008
*
* This is a reimplementation using FORTRAN-77 and LAPACK
* of spphase1.m from Prof. Yinyu Ye
* http://www.stanford.edu/~yyye/matlab/spphase2.m
*
* Note that the original MATLAB code uses MATLAB sparse matrix
* facilities, implemented in MATALB with the CHOLMOD routines
* from Tim Davis et al. LAPACK is for dense systems, this implementation
* uses just LAPACK so does not use sparseness (hence no 'sp' prefix
* on name, 'dl' is for 'dense LAPACK').
*
* This appears to be faster than the sparse version (spphs2.f) which
* uses UMFPACK though.
*
* $Id: dlphs2.f 2056 2009-01-21 23:54:02Z astivala $
*=======================================================================

      subroutine dlphs2(Q, ldq, A, lda, m, n, c, alpha, beta, lambda,
     $     toler, obvalue, ob,
     $     x, y, E, lde, u, ldu, xxvec, gg, XX, ldxx, 
     $     AA, ldaa, ipiv, wvec, Ewk, ldewk, work, lwork, info2)

      implicit none
*
* dlphs2 - phase 2 of dsoqp:
*          Repeatedly solve an ellipsoid constrained QP problem
*          by solving a linear system until a positive solution is
*          found.
*
*     .. Scalar Arguments ..
      integer            ldq,lda,m,n,lde,ldxx,ldaa,ldu,ldewk,lwork,info2
      double precision   alpha,beta,lambda,toler,obvalue,ob
*     ..
*     .. Array Arguments ..
      double precision   Q(ldq,*),A(lda,*),c(*),x(*),y(*),
     $     E(lde, *),u(*),xxvec(*),gg(*),XX(ldxx,*),
     $     AA(ldaa, *), wvec(*),Ewk(ldewk, *),work(*)
      integer            ipiv(*)

*     ..
*
* The phase 1 procedure called by dlsoqp.
* For more details see dlsoqp.f.
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
* E      (output) DOUBLE PRECISION, dimension (m+n, m+n)
*        workspace
*
* lde    (input) INTEGER
*        leading dimension of E. lde >= m+n
*
* u      (output) DOUBLE PREISION vector, dimension (m+n)
*        workspace
*
* ldu    (input) INTEGER
*        leading dimension of u, ldu >= m+n
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
* ipiv   (output) INTEGER, dimension (m+n)
*        workspace
*
* wvec    (output) DOUBLE PRECISION, dimension (n)
*        workspace     
*
* Ewk    (output) DOUBLE PRECISION, dimension (m+n, m+n)
*        workspace
*
* ldewk  (input) INTEGER
*        leading dimension of Ewk. lde >= m+n
*
* work   (output) DOUBLE PRECISION vector, dimension (lwork)
*        workspace for DSYSV
*
* lwork  (input) INTEGER
*        length of work vector. lwork >= (m+n)*NB where NB is optimal
*        blocksize for DSYTRF.
*
* info2  (output) INTEGER
*        = 0  : successful
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
      integer i,j,info,nrhs,incr,order,edim
      double precision w1,w2,blasa,blasb,nora,valgo
      character uplo
*     ..
*     .. Local Arrays ..

*     ..
*     .. Intrinsic Functions ..
      intrinsic min,abs
*     ..
*     .. External Subroutines and Functions ..
      external demvv
      double precision dvmin
      external dvmin
*     LAPACK subroutines
      external dlange,dpotrf,dpotrs,dpocon,dsysv
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
*     E*u = F where
*
*     E (m+n, m+n)  =    ( XXlam (n, n)   AA' (n, m) )
*                        ( AA    (m, n)    0  (m, m) )
*
*     F (m+n)       =    ( dx (n, 1) )
*                        ( 0  (m, 1) )
*
*     dx (n, 1) is elementwise -x.*gg where gg = Q*x+c
*     XXlam = XX + lambda*eye(n,n) where XX = dx * Q * dx
*                                  where dx(n) is x along diagonal
*
*     until a positive solution is found.
*
*     Inside the loop, only the top left submatrix of E is updated
*     (lambda is modified), other 3 submatrices are constant.
*     However since DSYSV stores the U*D*U**T
*     decomposition in this matrix parameter, must rebuild 
*     entire matrix each time in the loop.
*     Profiling showed doing this was the majority of time spent
*     in program, so rather than rebuilding entirely we build
*     it first outside the loop in another array,
*     then inside at each iteration
*     copy the one built outside the loop into working array E
*     using BLAS level 1 DCOPY on columns, then update the 
*     top left part by adding lambda to diagonal elements.
*     XX can be computed only once before the loop though,
*     and used in rebuilding E in the loop. This really made little
*     or no difference though, most time still spend building E
*     (but now outside loop) - need to use sparseness.
*     Part of F is constant 0 as well, but in fact we do not even
*     have F, we just use storage u for both input 'F' and output
*     u since LAPACK subroutine DSYSV stores solution in rhs
*     matrix. So we rebuid u on each iteration, on exit of
*     loop u contains the solution vector we want.      
*
*     XX (and therefore XXlam) is square, symmetric and sparse.
*     E is square, symmetric, and sparse. 
*     It may be positive definite, which can efficiently be tested
*     for with Cholesky factorization, and solved by that means
*     more efficiently than the general case using LU decomposition.
*
*     ==================================================================

      lambda = (1.0d0 - beta) * lambda

      valgo = 0.0d0

*     usually using increment of 1 in BLAS calls
      incr = 1

*     gg = Q*x + c
*     Q is symmetric  (using upper tri arbitrarily here)
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
       

*     E and its template Ewk are m+n square matrices
      edim = n + m

*     Build the Ewk matrix to be copied to E for input to DSYSV, since
*     DSYSV writes U*D*U**T decomposition over it so need top copy back
*     each iteration, then update by adding lambda. It is symmetric, we
*     will only build the upper triangle.
         do 130 j = 1, n
            call dcopy(j, XX(1,j), incr, Ewk(1,j), incr)
 130     continue
         do 170 j = 1, m
*           note use of INCX=ldaa for AA to get rows
            call dcopy(n, AA(j, 1), ldaa, Ewk(1,n+j), incr)
 170     continue
         do 190 j = n+1, edim
            do 180 i = n+1, j
               Ewk(i, j) = 0.0d0
 180        continue
 190     continue


*     Also since Q, AA, XX, XXlam all symmetric, should only store
*     one triangle. AA, XX, XXlam also are
*     all sparse so should also be more efficient to use specialized
*     routines rather than BLAS and LAPACK (as MATLAB does in the
*     MATLAB implementaton), but turns out actualy to be slower
*     when this was done using UMFPACK.

 200  if (valgo .le. 0.0d0) then

*        copy Ewk to E since previous E overwritten by DSYSV solution
*        and update E by adding lambda along diagonal elements of XX
*        submatrix in top left
*        only using upper triangle
         do 210 j = 1, n
            call dcopy(j, Ewk(1, j), incr, E(1, j), incr)
            E(j, j) = E(j, j) + lambda
 210     continue
         do 220 j = n+1, edim
            call dcopy(j, Ewk(1, j), incr, E(1, j), incr)
 220     continue
         
*        u(1:n) = -x .* gg
         call dcopy(n, x, incr, u, incr)
         call demvv(n, gg, u)
         blasa = -1.0d0
         call dscal(n, blasa, u, incr)
*        u(n+1:n+m) = 0
         do 300 i = n+1, edim
            u(i) = 0.0d0
 300     continue

C          write(*,380),( (E(i,j), j=1,m+n), i=1,m+n)
C  380     format (21(f8.4))
 
*     TODO: use DPOTRS to solve if postive definite, test with DPOTRF

         order = m + n
         nrhs = 1
         uplo = 'U'
         call dsysv(uplo, order, nrhs, E, lde, ipiv, u, ldu, 
     $        work, lwork, info)
         if (info .ne. 0) then
            write(*,*) 'dlphs2 DSYSV failed'
            info2 = info
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
         call dsymv(uplo, n, blasa, Q, ldq, x, incr, blasb, wvec, incr)
         ob = ddot(n, x, incr, wvec, incr) / 2.0d0 +
     $        ddot(n, c, incr, x, incr)
      endif

*     end of DLPHS2 , phase 2 subroutine of DLSOQP.

      info2 = 0
      return
      end

