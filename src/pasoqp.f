*=======================================================================
* File:    pasoqp.f
* Author:  Alex Stivala, based on spsolqp.m by Prof. Yinyu Ye.
* Created: July 2008
*
* This is a reimplementation using FORTRAN-77 and PARDISO
* of spsolqp.m from Prof. Yinyu Ye
* http://www.stanford.edu/~yyye/matlab/spsolqp.m

* This uses sparse symmetric linear solver PARDISO
* in the Intel Math Kernel Library (MKL) (10.1.0).
*
* LAPACK and BLAS are also required.
*
* $Id: pasoqp.f 2199 2009-04-03 06:33:11Z astivala $
*=======================================================================

      subroutine solvqp(Q, ldq, A, lda, m, n,  b, c, x, y, info)
!      use IFPORT
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
*     the first uses PAPHS1 find an interior feasible point and the
*     second uses PAPHS2 to find a local optimal solution.
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
*        = -3 : a PARDISO routine failed
*        = -4 : phase 2 failed
*        = -5 : max number of iterations exceeded
*        = -6 : numer of nonzero entries exceeds maximum dimensioned for
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

*     largest number of nonzero entries in sparse E matrix
      integer nzmax
      parameter(nzmax = 10000000)

*     limit on number of iterations
      integer maxiter
      parameter (maxiter = 1000)

*     -inf value (may be system dependent)
      double precision neginf
      parameter (neginf=-1.79769d+308)
*     ..
*     .. Common block Arrays ..
*     PARDISO internal solver memory pointer
      integer*8 pt(64)
*     lsyfac2 is true after paphs2 symbolic factorization done
*     lsyfac1 is true after paphs1 symbolic factorization done
      logical lsyfac2,lsyfac1
      common /pacommon/ pt,lsyfac2,lsyfac1
*     ..
*     .. Local Scalars ..
      character trans,uplo
      integer i,j,incr,order,iter,ldf,ldxx,ldaa,
     $     ixE
*     PARDISO variables
      integer maxfct,mnum,mtype,phase,nrhs,error,msglvl,idum
      double precision ddum
*     other local scalars
      double precision alpha,gap,ob,z,blasa,blasb,lambda,obvalue,
     $     nora,lower,xdotv,cdotx
*     ..
*     .. Local Arrays ..
      double precision avec(mmax),zhis(maxiter),comp(mmax+nmax)
*     PARDISO parameter block
      integer iparm(64)
*     sparse storage for matrix E for PARDISO in compressed row format
      integer Ep(nmax+mmax+1),Ei(nzmax)
      double precision Ex(nzmax)
*     workspace arrays for PAPHS1
*     F (each column of F separately), dx, ipiv,v reused for PAPHS2
*     v,ipiv,E and F are also used as working storage in PASOQP itself
      double precision dx(nmax), F(mmax+nmax, 2)
      double precision  u(nmax+2), v(nmax+2),
     $     umzv(nmax+2), obz(nmax+2), y1(mmax+nmax), y2(mmax+nmax)
*     additional workspace arrays for PAPHS2
      double precision XX(nmax,nmax),AA(mmax,nmax)
*     ..
*     .. Intrinsic Functions ..
!      intrinsic min,abs,max,dsqrt,dble,RANDOM_NUMBER
      intrinsic rand,min,abs,max,dsqrt,dble,srand
!      intrinsic min,abs,max,dsqrt,dble
*     ..
*     .. External Subroutines and Functions ..
      external paphs1,paphs2,drecip,demvv
      double precision dvmin,dvmax
      external dvmin,dvmax
      external rand55
*     PARDISO subroutines
      external pardiso
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


      ldf = mmax+nmax
      ldxx = nmax
      ldaa = mmax
      
C      call srand(1)

*     set parameters for PARDISO
      do 1,i = 1, 64
         iparm(i) = 0
 1    continue
      iparm(1) = 1 ! no solver default
      iparm(2) = 2 ! fill-in reordering from METIS
      iparm(3) = 1 ! numbers of processors
      iparm(4) = 0 ! no iterative-direct algorithm
      iparm(5) = 0 ! no user fill-in reducing permutation
      iparm(6) = 0 ! =0 solution on the first n compoments of x
      iparm(7) = 0 ! not in use
      iparm(8) = 2 ! numbers of iterative refinement steps
      iparm(9) = 0 ! not in use
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(12) = 0 ! not in use
      iparm(13) = 0 ! maximum weighted matching algorithm is switched-off
      iparm(14) = 0 ! Output: number of perturbed pivots
      iparm(15) = 0 ! not in use
      iparm(16) = 0 ! not in use
      iparm(17) = 0 ! not in use
      iparm(18) = 0 ! Output: number of nonzeros in the factor LU
      iparm(19) = 0 ! Output: Mflops for LU factorization
      iparm(20) = 0 ! Output: Numbers of CG Iterations
      error = 0 ! initialize error flag
      msglvl = 0 ! print statistical information
      mtype = -2 ! symmetric, indefinite
      maxfct = 1
      mnum = 1
      nrhs = 1
C.. Initiliaze the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
      do 2,i = 1, 64
         pt(i) = 0
 2    continue

      lsyfac2 = .false.
      lsyfac1 = .false.

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
         call paphs1(A, lda, m, n, avec, b, x, y, z, ob, alpha,
     $        dx, Ep, Ei, Ex, nzmax, F, ldf, u, v, umzv, obz, 
     $        y1, y2, info)
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
*     E * comp = y1 where
*
*     E (m+n, m+n)  =    ( eye (n, n)   A' (n, m) )
*                        ( A  (m, n)    0  (m, m) )
*
*     y1 (m+n)      =    ( r (n, 1) )
*                        ( 0 (m, 1) )
*
*     where r is an n-vector of psuedorandom number uniformly
*     distributed in [0,1].
*
*     E is square, symmetric, and sparse. But it is not necessarily
*     positive definite so we cannot always use CHOLMOD to efficiently
*     solve the system. We therefore use PARDISO with sparse matrix
*     E stored in compressed row format.
*     In this format, the nonzero entries are stored (row major)
*     in Ex. For each row the column indices of the nonzero entries
*     are stored in Ei. For each row k, Ep(k) and Ep(k+1)-1 are
*     respecitively the first and last indices of column numbers in Ei
*     and values in Ex for that column.
*     NB PARDISO uses 1-based indexing in Ep and Ei.
*     No duplicate row indices are allowed, and row indicies for each
*     column must be in ascending order.
*     Only the upper triangle is stored. Zeros on the diagonal must
*     be explicitly stored.
* 
*     If E were positive definite, which can efficiently be tested
*     for with Cholesky factorization, it could be solved by that means
*     more efficiently than the general case using LU decomposition.
*
*     TODO: use CHOLMOD to test for postive definiteness and use 
*           Cholesky factorization if it is.
*
*     ==================================================================

*     ixE is the current index (FORTRAN 1-based) in 
*     Ei and Ex.
      ixE = 1
      Ep(1) = 1
      do 730 i = 1, n
         Ei(ixE) = i
         Ex(ixE) = 1.0
         ixE = ixE + 1
         do 720 j = 1,m
            if (A(j,i) .ne. 0.0d0) then
               Ei(ixE) = n + j
               Ex(ixE) = A(j,i)
               ixE = ixE +1
               if (ixE .gt. nzmax) then
                  write(*,*) 'Max number of nonzero values exceeded.'
                  info = -6
                  return
               endif
            endif
 720     continue
         Ep(i+1) = ixE
 730  continue
*     need to explicitly store the remaining zeros on the diagonal
      do 760 i = n+1,n+m
         Ei(ixE) = i
         Ex(ixE) = 0.0d0
         ixE = ixE + 1
         Ep(i+1) = ixE
 760     continue
      
*     use y1  as RHS for linear system
      do 800 j = 1, n
!         CALL RANDOM_NUMBER(y1(j))
          y1(j) = dble(rand(0))
C         y1(j) = 0.5d0
 800  continue
      do 810 j = n+1, n+m
         y1(j) = 0.0d0
 810  continue


      order = m + n

*     numeric factorization and
*     solve E*comp=y1 with iterative refinement
*     NB we are reusing the symbolic factorization done in paphs1
      phase = 23
      CALL pardiso (pt, maxfct, mnum, mtype, phase, order, Ex, Ep, Ei,
     1 idum, nrhs, iparm, msglvl, y1, comp, error)
      if (error .ne. 0) then
         write(*,*) 'PARDISO factorization + solve failed, info = ',
     $        error
         info = -3
         return
      endif

*     nora = min(comp./x)
*     solution is in comp and using y2 as temp storage
      call dcopy(n, x, incr, y2, incr)
      call drecip(n, y2)
      call demvv(n, comp, y2)
      nora = dvmin(n, y2)
      
      if (nora .le. 0.0d0) then
         nora = -0.01d0 / nora
      else
         nora = dvmax(n, y2)
         if (nora .eq. 0.0d0) then
            write(*,*) 'The problem has a unique feasible point.'
            x(n+1) = ob
            info = 0
            return
         endif
         nora = 0.01d0 / nora
      endif
*     x = x + nora*comp
      call daxpy(n, nora, comp, incr, x, incr)
*     obvalue = (x'*(Q*x))/2 + c'*x
      uplo = 'U'
      blasa = 1.0d0
      blasb = 0.0d0
*     blasb is zero so v can be uninitialized on entry to DSYMV
      call dsymv(uplo, n, blasa, Q, ldq, x, incr, blasb, v, incr)
C      obvalue = ddot(n, x, incr, v, incr) / 2.0d0 +
C     $          ddot(n, c, incr, x, incr)
      xdotv = ddot(n, x, incr, v, incr)
      cdotx = ddot(n, c, incr, x, incr)
      obvalue = xdotv / 2.0d0 + cdotx


*     release internal memory from phase1 and our own PARDISO calls here
      phase = -1
      CALL pardiso (pt, maxfct, mnum, mtype, phase, order, Ex, Ep, Ei,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, error)
      lsyfac1 = .false.


      lower = neginf
      zhis(1) = lower
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
         
         call paphs2(Q, ldq, A, lda, m, n, c, 
     $        alpha, beta, lambda, toler,
     $        obvalue, ob, x, y, Ep,Ei,Ex, nzmax,
     $        y2, y1, dx, F,
     $        XX, ldxx, AA, ldaa,
     $        v, info)

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
*           comp = Q*x + c - A'*y
*           first compute comp <- Q*x + c, then comp <- -A'*y + comp
            uplo = 'U'
            blasa = 1.0d0
            blasb = 1.0d0
            call dcopy(n, c, incr, comp, incr)
            call dsymv(uplo,n,blasa,Q,ldq,x,incr,blasb,comp,incr)
            trans = 'T'
            blasa = -1.0d0
            blasb = 1.0d0
            call dgemv(trans, m, n, blasa, A, lda, y, incr, blasb, 
     $           comp, incr)

            if (dvmin(n, comp) .ge. 0.0d0) then
               zhis(iter+1) = ob - ddot(n, x, incr, comp, incr)
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

*     release internal memory
      phase = -1
      CALL pardiso (pt, maxfct, mnum, mtype, phase, order, Ex, Ep, Ei,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, error)
      lsyfac2 = .false.
      lsyfac1 = .false.

      x(n+1) = ob
      info = 0

*     end of PASOQP
      return
      end


