*=======================================================================
* File:    tmatd.f 
* Author:  Alex Stivala
* Created: July 2008
*
* FORTRAN-77 implementation of tableau matching (discrete)
* using the SOLVQP implementation
* of Prof. Yinyu Ye's interior point QP solver.
*
* The QP (relaxed from integer QP) we need to solve is
*
*     maximize:
*
*       \sum_{i,j,k,l} \zeta(\omega^A_{ik}, \omega^B_{jl})x_{ik}x_{jl}
*
*       where x_{ij} is a binary variable meaning SSE i in A is matched
*       with SSE j in B and \zeta is the tableau scoring function
*       (TSCORD)
*
*     subject to:
* 
*      \sum_{j=1}^{N_B} x_{ij} <= 1,    1 <= i <= N_A             (1)
*      \sum_{i=1}^{N_A} x_{ij} <= 1,    1 <= j <= N_B             (2)
*
* Constraints (1) and (2) ensure that each SSE in one tableau is 
* matched with at most one SSE in the other.
* This amounts to N_A+N_B constraints.
*
* We may optionally add the further constraints:
*
*     x_{ik} = 0,  if \omega^A_{ii} <> \omega^B_{kk}              (3)
*                  1 <= i <= N_A, 1 <= k <= N_B
*
* This amounts to (up to) a further N_A*N_B constraints.
*
*     x_{ij} + x_{kl} <= 1,  1 <= i < k < N_A, 1 <= l < j <= N_B  (4)
*
* This amounts to a further ((N_A*(N_A-1))*(N_B*(N_B-1)))/4 
* constraints. i.e. a  lot (quartic)
*
*
*     x_{ij} + x_{kl} <= 1,  if |d^A_{ik} - d^B_{jl}| > T         (5)
*                  1 <= i,k <= N_A, 1 <= j,l <= N_B
* 
*       where d^A and d^B are SSE distance matrices and T is
*       distance difference threshold.
*
* This amounts to (up to) a further N_A*N_B*N_A*N_B constraints.
*
* Constraint (3) disallows matchings of SSEs that are not of the 
* same type (i.e. strand or helix), by using the SSE type information
* stored on the main diagonal of the tableau matrix.
* Constraint (4) disallows "out of order" matchings, i.e. prevents
* non-linear matchings.
* Constraint (5) disallows matchings between SSEs where the 
* difference between the distances between the SSEs (midpoints)
* is larger than the treshold T.
*
*
* SOLVQP (dense version DLOSQP or sparse version SPSOQP)
* solves the quadratic program in standard form:
*
*     minimize    0.5*(x'*Q*x)+c'*x
*     subject to  A*x=b, x>=0.
*
* To convert the maximization problem to minimization we simply negate
* the Q matrix (c is all zero so irrelevant).
* Note also the constraints are inequality constraints, we convert to
* equality by introducing slack variables eg constraint (1) is
* expressed as
*    
*      \sum_{j=1}^{N_B} ( x_{ij} + s_i ) = 1,    1 <= i <= N_A 
* 
* where s_i are the N_A 'slack' variables. Similarly for constraints
* (2) (adding N_B) slack variables, (4) (adding a further
* ((N_A*(N_A-1))*(N_B*(N_B-1)))/4 slack variables)
* and (5) (adding up to a further N_A*N_B*N_A*N_B slack variables)
*
* Note using all of these constraints results in a very large number
* of constraints (and slack variables), which therefore also increases
* order of Q matrix)),
* so much so that the A matrix (let alone Q as well) if fully stored
* is too large for virtual memory / 32 bit addressing, so most
* constraints are actually relaxed and implemented as penalties
* in object function (Q matrix values) instead.
*
*
* $Id: tmatd.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================

      subroutine tmatd(tab1,ldt1,n1,tab2,ldt2,n2,
     $     ltype,lorder,dmat1,dmat2,score,soln,info)

      implicit none
*
* tmatd - find matching score of two tableaux using relaxed QP
*
*   The tableaux are specified here as arrays of 2 char tableau codes
*

*
*     .. Scalar Arguments ..
      integer ldt1,n1,ldt2,n2,info
      logical ltype,lorder
      double precision score
      double precision soln(*)
*     ..
*     .. Array Arguments ..
      character*2   tab1(ldt1,*),tab2(ldt2,*)
      double precision dmat1(ldt1,*),dmat2(ldt2,*)
*     ..
*

*
* Arguments
* =========
*
* tab1  (input) CHARACHTER*2 array, dimension (n1,n1)
*        Tableau for one structure. Symmetric.
*
* ldt1   (input) INTEGER
*        Leading dimension of tab1 array. ldt1 >= n1
*
* n1     (input) INTEGER
*        Dimension of tab1 array.
*
* tab2  (input) CHARACTER*2 array, dimension (n2,n2)
*        Tableau for second structure. Symmetric.
*
* ldt2   (input) INTEGER
*        Leading dimension of tab2 array. ldt2 >= n2
*
* n2     (input) INTEGER
*        Dimension of tab2 matrix.
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
*        indices in tab1 and j,l indices in tab2.
*
* dmat1  (input) DOUBLE PRECISION array, dimension (n1,n1)
*        SSE distance matrix for one structure. symmetric
*
* dmat2  (input) DOUBLE PRECISION array, dimension (n2,n2)
*        SSE distance matrix for second structure. symmetric
*
* score  (output) DOUBLE PRECISION
*        score of matching tab1 with tab2 using relaxed QP.
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

*     largest tableau  allowed
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
      character*2 tc1,tc2
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
      intrinsic dble,abs
*     ..
*     .. External Subroutines and Functions ..
      external dcopy
      external solvqp
      integer tscord
      external tscord

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

*
* 
*     Set up A matrix and b vector for constraints (1) and (2)
*     which ensure that each SSE in one tableau is 
*     matched with at most one SSE in the other.
*
*      \sum_{j=1}^{N_B} x_{ij} <= 1,    1 <= i <= N_A             (1)
*      \sum_{i=1}^{N_A} x_{ij} <= 1,    1 <= j <= N_B             (2)
*
*
*

*     set b to ones, since all contraints are <= 1
      do 30 j = 1, n1+n2
         b(j) = 1.0d0
 30   continue
*
*     setup the constraint matrix A
*
*     constraint (1)
      do 50 j = 1, n2
         do 40 i = 1, n1
            A(i, (i-1)*n2+j) = 1.0d0
 40      continue
 50   continue
*     constraint (2)
      do 70 j = 1, n1
         do 60 i = 1, n2
            A(n1+i, (j-1)*n2+i) = 1.0d0
 60      continue
 70   continue
*     slack variables for constraints (1) and (2)
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
     $                   (tab1(i,i) .ne. tab2(j,j) .or.
     $                    tab1(k,k) .ne. tab2(l,l))) then
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
                           tc1 = tab1(i,k)
                           tc2 = tab2(j,l)
                           Q((i-1)*n2 + j, (k-1)*n2 + l) = -1.0d0 *
     $                          dble(tscord(tc1, tc2))
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

