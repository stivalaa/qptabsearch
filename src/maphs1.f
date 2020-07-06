*=======================================================================
* File:    maphs1.f
* Author:  Alex Stivala, based on spphase1.m by Prof. Yinyu Ye.
* Created: July 2008
*
* This is a reimplementation using FORTRAN-77 and MA57
* of spphase1.m from Prof. Yinyu Ye
* http://www.stanford.edu/~yyye/matlab/spphase1.m
*
* This uses sparse symmetric linear solver MA57 from HSL 2007
* http://hsl.rl.ac.uk/hsl2007/distrib/hsl2007.html
* needs registration as a registred researcher
*
* $Id: maphs1.F 3032 2009-12-10 05:22:56Z alexs $
*=======================================================================

      subroutine maphs1(A, lda, m, n, avec, b, x, y, z, ob, alpha,
     $     dx, Eirn, Ejcn, Ea, lkeep, keep, 
     $     nzmax, F, ldf, u, v, umzv, obz,
     $     y1, y2,
     $     info)

      implicit none
*
* maphs1 - phase 1 of masoqp: find a feasible point
*
*     .. Scalar Arguments ..
      integer            lda,m,n,nzmax,ldf,lkeep
      double precision   z,ob,alpha
*     ..
*     .. Array Arguments ..
      double precision   A(lda,*), x(*), y(*), avec(*), b(*)
      double precision   dx(*),Ea(*),F(ldf, *),
     $     u(*), v(*), umzv(*), obz(*), y1(*), y2(*)
      integer            Eirn(*),Ejcn(*),keep(*)

*     ..
*
* The phase 1 procedure called by masoqp.
* For more details see masoqp.f.
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
* Eirn   (output) INTEGER vector, dimension (nzmax)
*        workspace
*        sparse storage for matrix E for MA47
*
* Ejcn  (output) INTEGER vector, dimension (nzmax)
*        workspace
*        sparse storage for matrix E for MA57
*
* Ea     (output) DOUBLE PRECICISION vector, dimension (nzmax)
*        workspace
*        sparse storage for matrix E for MA57
*
* lkeep  (input) INTEGER   
*        length of keep array for MA57.
*        must be >= 5*maxorder+nzmax+max(maxorder,nzmax)+42
*         where maxorder = nmax*mmax
*
* keep   (output) INTEGER vector, dimension (lkeep)
*        workspace for MA57
*
* nzmax  (input) INTEGER
*        maximum allowed number of nonzero values in E. 
*        Dimensions Eirn and Ejcn and Ea.
*
* F      (output) DOUBLE PRECISION, dimension (m+n, 2)
*        workspace
*
* ldf    (input) INTEGER
*        leading dimension of F. ldf >= m+n
*
* y1     (output) DOUBLE PRECISION vector, dimension(m+n)
*        workspace
*        solution from MA57 with rhs column 1 of F
*
* y2     (output) DOUBLE PRECISION vector, dimension(m+n)
*        workspace
*        solution from MA57 with rhs column 2 of F
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
*       other  : info code from MA57
*
*=======================================================================
*

*     .. Parameters ..
      integer mmax,nmax,lwork,lfact,lifact
      parameter(mmax = 220, nmax = 4000)
      integer iworklen
      parameter(iworklen = 5*mmax*nmax)
      parameter(lwork = mmax*nmax)
      parameter(lfact = 5*mmax*nmax)
      parameter(lifact = 5*mmax*nmax)
*     ..
*     .. Common block Arrays ..
*     lsyfac2 is true after maphs2 symbolic factorization done
*     lsyfac1 is true after maphs1 symbolic factorization done
      logical lsyfac2,lsyfac1
      common /macommon/ lsyfac2,lsyfac1
*     ..
*     .. Local Scalars ..
      integer i,j,order,info,incr,uvlen,xlen,ixE
*     MA57 variables
      integer ne,job
*     other local scalars
      double precision w1,w2,obrec,obrec2,adoty2,blasa,blasb,nora,vn2rec
      character trans
C      character uplo
*     ..
*     .. Local Arrays ..
*     MA57 workspace
      integer iwork(iworklen)
      integer icntl(20),mainfo(40)
      double precision cntl(5),rinfo(20)
      double precision work(lwork)
      integer ifact(lifact)
      double precision fact(lfact)
*     ..
*     .. Intrinsic Functions ..
      intrinsic min,int
*     ..
*     .. External Subroutines and Functions ..
      external drecip,demvv
      double precision dvmin,dvmax
      external dvmin,dvmax
*     MA57 subroutines
      external MA57ID,MA57AD,MA57BD,MA57CD
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
*     solve the system. We therefore use MA57 with sparse matrix
*     E stored in the HSL2007 MA57 format.
*     In this format, the i and j indices of the nonzero entries of E
*     are stored in Eirn and Ejcn respectively. (Only one triangle
*     is stored). The nonzero values are stored in Ea, corresponding
*     to the i,j indices in Eirn and Ejcn.
*
*     If E were positive definite, which can efficiently be tested
*     for with Cholesky factorization, it could be solved by that means
*     more efficiently than the general case using LU decomposition.
*
*     TODO: use CHOLMOD to test for postive definiteness and use if
*           it is.
*
*     ==================================================================

*     set parameters for MA57
      call MA57ID(cntl, icntl)
*      icntl(5) = 4  ! debug messages enabled

      incr = 1
      call dcopy(n, x, incr, dx, incr)
*     seems to be no BLAS subroutine for elementwise reciprocal, use own
      call drecip(n, dx)

*     ixE is the current index in Eirn and Ejcn and Ea
      ixE = 1
      do 30 i = 1, n
         Eirn(ixE) = i
         Ejcn(ixE) = i
         Ea(ixE) = dx(i)**2
         ixE = ixE + 1
         do 20 j = 1, m
            if (A(j,i) .ne. 0.0d0) then
               Eirn(ixE) = i
               Ejcn(ixE) = n + j
               Ea(ixE) = A(j,i)
               ixE = ixE +1
               if (ixE .gt. nzmax) then
                  write(*,*) 'Max number of nonzero values exceeded.'
                  info = -6
                  return
               endif
            endif
 20      continue
 30   continue
      ne = ixE - 1

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
         CALL MA57AD (order, ne, Eirn, Ejcn, lkeep, keep, iwork, 
     $        icntl, mainfo, rinfo)
         if (mainfo(1) .lt. 0) then
            write(*,*) 'MA57 symbolic analysis failed, info=',mainfo(1)
            info = mainfo(1)
            return
         endif
         lsyfac1 = .true.
      endif
*     numeric factorizatino and
*     solve E*y1=F  without iterative refinement
      CALL MA57BD (order, ne, Ea, fact, lfact, ifact, lifact,
     $     lkeep, keep, iwork,
     $     icntl, cntl, mainfo, rinfo)
      if (mainfo(1) .lt. 0) then
         write(*,*) 'MA57 factorization failed, info = ', mainfo(1)
         info = mainfo(1)
         return
      endif
      job = 0
      CALL MA57CD (job, order, fact, lfact, ifact, lifact,
     $     2, F, ldf, work, lwork, iwork, icntl, mainfo)
      if (mainfo(1) .lt. 0) then
         write(*,*) 'MA57 solve [1] failed, info = ', mainfo(1)
         info = mainfo(1)
         return
      endif
      call dcopy(order, F, incr, y1, incr)
      call dcopy(order, F(1,2) ,incr, y2, incr)

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

