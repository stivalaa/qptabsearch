*=======================================================================
* File:    paphs1.f
* Author:  Alex Stivala, based on spphase1.m by Prof. Yinyu Ye.
* Created: July 2008
*
* This is a reimplementation using FORTRAN-77 and PARDISO
* of spphase1.m from Prof. Yinyu Ye
* http://www.stanford.edu/~yyye/matlab/spphase1.m
*
* This uses sparse symmetric linear solver PARDISO
* in the Intel Math Kernel Library (MKL) (10.1.0).
*
* $Id: paphs1.f 2199 2009-04-03 06:33:11Z astivala $
*=======================================================================

      subroutine paphs1(A, lda, m, n, avec, b, x, y, z, ob, alpha,
     $     dx, Ep, Ei, Ex, nzmax, F, ldf, u, v, umzv, obz,
     $     y1, y2,
     $     info)

      implicit none
*
* paphs1 - phase 1 of pasoqp: find a feasible point
*
*     .. Scalar Arguments ..
      integer            lda,m,n,nzmax,ldf
      double precision   z,ob,alpha
*     ..
*     .. Array Arguments ..
      double precision   A(lda,*), x(*), y(*), avec(*), b(*)
      double precision   dx(*),Ex(*),F(ldf, *),
     $     u(*), v(*), umzv(*), obz(*), y1(*), y2(*)
      integer            Ep(*),Ei(*)

*     ..
*
* The phase 1 procedure called by pasoqp.
* For more details see pasoqp.f.
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
*        sparse storage for matrix E for PARDISO
*
* Ei     (output) INTEGER vector, dimensino (nzmax)
*        workspace
*        sparse storage for matrix E for PARDISO
*
* Ex     (output) DOUBLE PRECICISION vector, dimension (nzmax)
*        workspace
*        sparse storage for matrix E for PARDISO
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
*        solution from PARDISO with rhs column 1 of F
*
* y2     (output) DOUBLE PRECISION vector, dimension(m+n)
*        workspace
*        solution from PARDISO with rhs column 2 of F
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
* info   (output) INTEGER
*        = 0  : successful
*        = 1  : max number of nonzero values exceeded, increase nzmax
*       other  : info code from PARDISO
*
*=======================================================================
*

*     .. Parameters ..

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
      integer i,j,order,info,incr,uvlen,xlen,ixE
*     PARDISO variables
      integer maxfct,mnum,mtype,phase,nrhs,error,msglvl,idum
      double precision ddum
*     other local scalars
      double precision w1,w2,obrec,obrec2,adoty2,blasa,blasb,nora,vn2rec
      character trans
C      character uplo
*     ..
*     .. Local Arrays ..
*     PARDISO parameter block
      integer iparm(64)
*     ..
*     .. Intrinsic Functions ..
      intrinsic min,int
*     ..
*     .. External Subroutines and Functions ..
      external drecip,demvv
      double precision dvmin,dvmax
      external dvmin,dvmax
*     PARDISO subroutines
      external pardiso
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
*     solve the system. We therefore use PARDISO with sparse matrix
*     E stored in compressed row format.
*     In this format, the nonzero entries are stored (row major)
*     in Ex. For each row the column indices of the nonzero entries
*     are stored in Ei. For each row k, Ep(k) and Ep(k+1)-1 are
*     respecitively the first and last indices of column numbers in Ei
*     and values in Ex for that row.
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
*     TODO: use CHOLMOD to test for postive definiteness and use if
*           it is.
*
*     ==================================================================

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
      iparm(27 ) = 0   ! check parameters if 1
      error = 0 ! initialize error flag
      msglvl = 0 ! print statistical information
      mtype = -2 ! symmetric, indefinite
      maxfct = 1
      mnum = 1
      nrhs = 1

      incr = 1
      call dcopy(n, x, incr, dx, incr)
*     seems to be no BLAS subroutine for elementwise reciprocal, use own
      call drecip(n, dx)

*     ixE is the current index (FORTRAN 1-based) in Ei and Ex.
      ixE = 1
      Ep(1) = ixE
      do 30 i = 1, n
         Ei(ixE) = i
         Ex(ixE) = dx(i)**2
         ixE = ixE + 1
         do 20 j = 1, m
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
 20      continue
         Ep(i+1) = ixE
 30   continue
*     need to expclitly store the remaining zeros on the diagonal
      do 40 i = n+1,n+m
         Ei(ixE) = i
         Ex(ixE) = 0.0d0
         ixE = ixE + 1
         Ep(i+1) = ixE
 40   continue


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

*     reordering and symbolic factorization and
*     numeric factorization
*     We do this once on the first iteration of loop in first call
*     of this function, it is then reused on subsequent 
*     iterations + calls as the structure of the matrix does not change
*     and analysis is an expensive operation
      if (.not. lsyfac1) then
         phase = 11
         CALL pardiso (pt, maxfct, mnum, mtype, phase, order,Ex,Ep,Ei,
     1    idum, nrhs, iparm, msglvl, ddum, ddum, error)
         if (error .ne. 0) then
            write(*,*) 'PARDISO reordering+symfac failed, info = ',error
            info = error
            return
         endif
         lsyfac1 = .true.
      endif
*     numeric factorizatino and
*     solve E*y1=F(:,1) with iterative refinement
      phase = 23
      CALL pardiso (pt, maxfct, mnum, mtype, phase, order, Ex, Ep, Ei,
     1 idum, nrhs, iparm, msglvl, F, y1, error)
      if (error .ne. 0) then
         write(*,*) 'PARDISO back subst [1] failed, info = ', error
         info = error
         return
      endif
*     solve E*y2=F(:,2) with iterative refinement
      phase = 33
      CALL pardiso (pt, maxfct, mnum, mtype, phase, order, Ex, Ep, Ei,
     1 idum, nrhs, iparm, msglvl, F(1,2), y2, error)
      if (error .ne. 0) then
         write(*,*) 'PARDISO back subst [2] failed, info = ', error
         info = error
         return
      endif
*     NB we do not free PARDISO memory as we reuse the symbolic
*     factorization done here back in the caller pasoqp, where
*     the matrix has the same structure.

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

