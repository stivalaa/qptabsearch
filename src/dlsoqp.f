*=======================================================================
* File:    dlsoqp.f
* Author:  Alex Stivala, based on spsolqp.m by Prof. Yinyu Ye.
* Created: July 2008
*
* This is a reimplementation using FORTRAN-77 and LAPACK
* of spsolqp.m from Prof. Yinyu Ye
* http://www.stanford.edu/~yyye/matlab/spsolqp.m
*
* Note that the original MATLAB code uses MATLAB sparse matrix
* facilities, implemented in MATALB with the CHOLMOD routines
* from Tim Davis et al. LAPACK is for dense systems, this implementation
* uses just LAPACK so does not use sparseness (hence no 'sp' prefix
* on name, 'dl' is for 'dense LAPACK').
*
* $Id: dlsoqp.F 3240 2010-01-18 03:28:54Z alexs $
*=======================================================================

      subroutine solvqp(Q, ldq, A, lda, m, n,  b, c, x, y, info)
!      USE IFPORT
      implicit none
*
* solvqp - solve quadratic program in standard form
*
*  This program solves quadratic program in standard form:
*
*     minimize    0.5*(x'*Q*x)+c'*x
*     subject to  A*x=b, x>=0.
*
*     This program is the implementation of the interior ellipsoidal
*     trust region and barrier function algorithm with dual solution
*     updating technique in the standard QP form. Two phases are used:
*     the first uses DLPHS1 find an interior feasible point and the
*     second uses DLPHS2 to find a local optimal solution.
*
*  Technical Reference
*  
*     Y. Ye, "An extension of Karmarkar's algorithm and the trust region
*     method for convex quadratic programming," in Progress in
*     Mathematical Programming (N. Megiddo ed.), Springer-Verlag, NY
*     (1989) 49-63.
*
*     Y. Ye, "On affine-scaling algorithms for nonconvex quadratic
*     programming," Math. Programming 56 (1992) 285-300.
*
* See also
*
*      Y. Ye.  Interior Point Algorithms: Theory and Analysis.
*      Wiley-Interscience 1997.
*
*  Comment: Each iteration we solve a linear KKT system like
*
*  ( Q+mu X^{-2}   A^T )(dx) = c'
*  (     A          0  )(dy) = 0
*
*  where X = diag(x)  which is a positive diagonal matrix.

*
*     .. Scalar Arguments ..
      integer            ldq,lda,m,n,info
*     ..
*     .. Array Arguments ..
      double precision   Q(ldq,*), A(lda,*), x(*), y(*), b(*), c(*)
*     ..
*

*
* Arguments
* =========
*
* Q      (input) DOUBLE PRECISION array, dimension (n, n)
*        the symmetric objective matrix
* 
* ldq    (input) INTGER
*        leading dimension of Q. ldq >= n
*
* A      (input) DOUBLE PRECISION array, dimension (m, n)
*        the constraint left-hand matrix
*
* lda    (input) INTEGER
*        leading dimenion of A. lda >= m
*
* m      (input) INTEGER
*        number of rows in matrix A
*        m <= n
*
* n      (input) INTEGER
*        number of rows and columns in marix Q, also the
*        number of columns in matrix A
*
* b      (input) DOUBLE PRECISION vector, dimension(m)
*        constraint right-hand column vector
*
* c      (input) DOUBLE PRECISION vector, dimension(m)
*        objective column vector
*
* x      (output) DOUBLE PRECISION vector, dimension (n+1)
*        solution vector with objective value appended as last element
*
* y      (output) DOUBLE PRECISION vector, dimension (m)
*        optimal dual solution (Lagrangian multiplier)
*
* info   (output) INTEGER
*        on exit, status of the computation
*        =  0 : successful exit
*        = -1 : phase 1 failed 
*        = -2 : sizes exceed maximum dimensions 
*        = -3 : LAPACK DSYSV failed
*        = -4 : phase 2 failed
*        = -5 : max number of iterations exceeded
*        =  1 : no feasible point
*        =  2 : problem appears unbounded  
*
*=======================================================================
*

*     .. Parameters ..

*     relative stopping tolerance: the objective value close to the
*     local optimum in the range of tolerance.
      double precision toler
      parameter (toler = 1.0d-5)

*     step size 0 < beta < 1
      double precision beta
      parameter (beta = 0.8d0)

*     largest values of array dimensions m and n allowed
      integer mmax,nmax
      parameter(mmax = 220, nmax = 4000)

*     limit on number of iterations
      integer maxiter
      parameter (maxiter = 1000)

*     -inf value (may be system dependent)
      double precision neginf
      parameter (neginf=-1.79769d+308)
*     ..
*     .. Local Scalars ..
      character trans,uplo
      integer i,j,incr,order,nrhs,iter,lde,ldf,ldxx,ldaa,ldewk,lwork
      double precision alpha,gap,ob,z,blasa,blasb,lambda,obvalue,
     $     nora,lower
*     ..
*     .. Local Arrays ..
      double precision avec(mmax),zhis(maxiter)
*     workspace arrays for DLPHS1
*     E, F (each column of F separately), dx, ipiv,v reused for DLPHS2
*     v,ipiv,E and F are also used as working storage in DLSOQP itself
      double precision dx(nmax),E(mmax+nmax, mmax+nmax),
     $     F(mmax+nmax, 2), work(mmax+nmax)
      integer ipiv(mmax+nmax)
      double precision  u(nmax+2), v(nmax+2),
     $     umzv(nmax+2), obz(nmax+2)
*     additional workspace arrays for DLPHS2
      double precision XX(nmax,nmax),AA(mmax,nmax),
     $     Ewk(mmax+nmax, mmax+nmax)
*     ..
*     .. Intrinsic Functions ..
      intrinsic min,abs,max,dsqrt,dble
#if defined(__INTEL_COMPILER) || defined(__PORTLAND_COMPILER)
      intrinsic RANDOM_NUMBER
#else
      intrinsic rand,srand
#endif
*     ..
*     .. External Subroutines and Functions ..
      external dlphs1,dlphs2,drecip,demvv
      double precision dvmin,dvmax
      external dvmin,dvmax
      external rand55
*     LAPACK subroutines
      external dsysv
*     BLAS level 1 subroutines/functions
      external dcopy,dscal,daxpy
      double precision ddot
      external ddot
*     BLAS level 2 subroutines
      external dgemv,dsymv
*     ..
*     .. Executable Statements ..
*

      if  (m .gt. mmax) then
         write(*,*) 'maximum dimensions exceeded, increase mmax'
         info = -2
         return
      endif
      if  (n .gt. nmax) then
         write(*,*) 'maximum dimensions exceeded, increase nmax'
         info = -2
         return
      endif


      lde = mmax+nmax
      ldf = mmax+nmax
      ldxx = nmax
      ldaa = mmax
      ldewk = mmax+nmax
      
*     length of work vector for DSYSV >= (nmax+mmax)*NB where NB is
*     optimal blocksize for DSYSTRF. NB we are using NB=1, this
*     was optimal on x86 Linux anyway.
      lwork = mmax+nmax

C      call srand(1)

* ----------------------------------------------------------------------
*
*     phase 1: try to find an interior feasible point
*
* ----------------------------------------------------------------------

      alpha = 0.95d0
      do 10 i = 1, n+1
         x(i) = 1.0d0
 10   continue
      ob = x(n+1)
      z = 0.0d0
      
*     avec = b - A*ones(n,1) (note using x as ones here)
      trans = 'N'
      blasa = -1.0d0
      blasb = 1.0d0
      incr = 1
      call dcopy(m, b, incr, avec, incr)
      call dgemv(trans, m, n, blasa, A, lda, x, incr, blasb,
     $     avec, incr)

      gap = ob - z
 500  if (gap .ge. toler) then
         call dlphs1(A, lda, m, n, avec, b, x, y, z, ob, alpha,
     $        dx, E, lde, F, ldf, ipiv, u, v, umzv, obz, 
     $        work, lwork, info)
         if (info .ne. 0) then
            write(*,*) 'phase 1 failed'
            info = -1
            return
         endif

         ob = x(n+1)
         gap = ob - z
         if (z .gt. 0.0d0) then
            gap  = -1.0d0
            write(*,*) 'The system has no feasible solution.'
            info = 1
            return
         endif
         goto 500
      endif

* ----------------------------------------------------------------------
*
*     phase 2: search for an optimal solution
*
* ----------------------------------------------------------------------

*     only using first n elements of x from now on

*     always using increment of 1 in BLAS calls
      incr = 1

      alpha =  0.9d0

*     ==================================================================
*
*     Solve the linear system
*
*     E * comp = F where
*
*     E (m+n, m+n)  =    ( eye (n, n)   A' (n, m) )
*                        ( A  (m, n)    0  (m, m) )
*
*     F (m+n)       =    ( r (n, 1) )
*                        ( 0 (m, 1) )
*
*     where r is an n-vector of psuedorandom number uniformly
*     distributed in [0,1].
*
*     E is square, symmetric, and sparse. 
*     It may be positive definite, which can efficiently be tested
*     for with Cholesky factorization, and solved by that means
*     more efficiently than the general case using LU decomposition.
*
*
*     ==================================================================

*      symmetric matrix so store only upper triangle and use DSYSV
*     (tried using packed format and DPPSV but that is slower)
      do 730 j = 1, n
         do 720 i = 1, n
            if (i .eq. j) then
               E(i,j) = 1.0d0
            else
               E(i,j) = 0.0d0
            endif
 720     continue
 730  continue
      do 760 j = 1,m 
*        note using INCX=lda for A to get rows
         call dcopy(n, A(j, 1), lda, E(1, n+j), incr)
 760  continue
      do 790 j = n+1, n+m
         do 780 i = n+1, j
            E(i, j) = 0.0d0
 780     continue
 790  continue

*     use first column of F as RHS for linear system, also solution
      do 800 j = 1, n
#if defined(__INTEL_COMPILER) || defined(__PORTLAND_COMPILER)
         CALL RANDOM_NUMBER(F(j,1))
#else
         F(j,1) = dble(rand(0))
#endif
 800  continue
C      call rand55(n, F)
      do 810 j = n+1, n+m
         F(j, 1) = 0.0d0
 810  continue

C       write(*,888),( (E(i,j), j=1,m+n), i=1,m+n)
C  888  format (21(f8.4))
C       write(*,890),(F(i, 1), i = 1, m+n)
C  890  format(f8.4)

      order = m + n
      nrhs = 1
      uplo = 'U'
      call dsysv(uplo, order, nrhs, E, lde, ipiv, F, ldf, 
     $           work, lwork, info)
      if (info .ne. 0) then
         info = -3
         return
      endif
*     nora = min(F(1:n,1)./x)
*     DSYSV solution is in F(:,1) and using F(:,2) as temp storage
      call dcopy(n, x, incr, F(1,2), incr)
      call drecip(n, F(1,2))
      call demvv(n, F, F(1,2))
      nora = dvmin(n, F(1,2))
      
      if (nora .le. 0.0d0) then
         nora = -0.01d0 / nora
      else
         nora = dvmax(n, F(1,2))
         if (nora .eq. 0.0d0) then
            write(*,*) 'The problem has a unique feasible point.'
            x(n+1) = ob
            info = 0
            return
         endif
         nora = 0.01d0 / nora
      endif
*     x = x + nora*F(1:n,1)
      call daxpy(n, nora, F, incr, x, incr)
*     obvalue = (x'*(Q*x))/2 + c'*x
      uplo = 'U'
      blasa = 1.0d0
      blasb = 0.0d0
*     blasb is zero so v can be uninitialized on entry to DSYMV
      call dsymv(uplo, n, blasa, Q, ldq, x, incr, blasb, v, incr)
      obvalue = ddot(n, x, incr, v, incr) / 2.0d0 +
     $          ddot(n, c, incr, x, incr)

      lower = neginf
      gap = 1.0d0
      lambda = max(1.0d0, abs(obvalue)/dsqrt(dsqrt(dble(n))))
      iter = 0
 900  if (gap .ge. toler) then
         iter = iter + 1
         if (iter .ge. maxiter) then
            write(*,*) 'Maximum number of iterations exceeded.'
            x(n+1) = ob
            info = -5
            return
         endif
         
         call dlphs2(Q, ldq, A, lda, m, n, c, 
     $        alpha, beta, lambda, toler,
     $        obvalue, ob, x, y, E, lde, F(1,1), ldf, dx, F(1,2),
     $        XX, ldxx, AA, ldaa, ipiv, 
     $        v, Ewk, ldewk, work, lwork, info)

         if (info .ne. 0) then
            write(*,*) 'phase 2 failed'
            info = -4
            return
         endif

C          write(*,910),iter,ob
C  910     format('iter = ',i5,' ob = ',f8.4)

         if (ob .eq. neginf) then
            gap = 0.0d0
            write(*,*) 'The problem is unbounded.'
            info = 2
            return
         else
*           F(1:n,1) = Q*x + c - A'*y
*           first compute F(1:n,1) = Q*x + c, then -A'*y + F(1:n,1)
            uplo = 'U'
            blasa = 1.0d0
            blasb = 1.0d0
            call dcopy(n, c, incr, F, incr)
            call dsymv(uplo, n, blasa, Q, ldq, x, incr, blasb, F, incr)
            trans = 'T'
            blasa = -1.0d0
            blasb = 1.0d0
            call dgemv(trans, m, n, blasa, A, lda, y, incr, blasb, 
     $           F, incr)

            if (dvmin(n, F) .ge. 0) then
               zhis(iter+1) = ob - ddot(n, x, incr, F, incr)
               lower = zhis(iter+1)
               gap = (ob-lower)/(1.0d0+abs(ob))
               obvalue = ob
            else
               zhis(iter+1) = zhis(iter)
               lower = zhis(iter+1)
               gap = (obvalue - ob) / (1.0d0 + abs(ob))
               obvalue = ob
            endif
         endif
         goto 900
*     end of while (gap .ge. toler) loop
*     a (local) optimal solution has been found
      endif

*     end of DLSOQP
      x(n+1) = ob
      info = 0
      return
      end


