* testpardiso.f
* FORTRAN-77 test of PARDISO sparse symmetric direct solver
* interface.
* This is basically pardiso_sym_f.f example program from PARDSIO
* http://www.pardiso-project.org
* NB however that this program was tested using the PARDISO included in
* the Intel Math Kernel Library (MKL) not standalone PARDISO.
*
* Tested using Intel MKL included with Intel compilers version 11.0
* (Linux 64 bit).
*
* $Id: testpardiso.f 2198 2009-04-03 05:41:34Z astivala $
      program testpardiso
      implicit none

      external pardiso
*     internal solver memory pointer
      integer*8 pt(64)
*     other variables
      integer n,nz,maxfct,mnum,mtype,phase,nrhs,error,msglvl
      integer idum,i
      real*8  ddum
      integer iparm(64)
      parameter (n=8,nz=18,nrhs=1,maxfct=1,mnum=1)

      integer Ap(n+1), Ai(nz)
      real*8 Ax(nz), x(n), b(n)


C.. Fill all arrays containing matrix data.

      DATA Ap /1,5,8,10,12,15,17,18,19/
      DATA Ai
     1 /1,  3,    6,7,
     2    2,3,  5,
     3      3,        8,
     4        4,    7,
     5          5,6,7,
     6            6,  8,
     7              7,
     8                8/
      DATA Ax
     1 /7.d0,     1.d0,          2.d0,7.d0,
     2       -4.d0,8.d0,     2.d0,
     3            1.d0,                    5.d0,
     4                 7.d0,     9.d0,
     5                      5.d0,1.d0,5.d0,
     6                           -1.d0,     5.d0,
     7                                11.d0,
     8                                     5.d0/


C..
C.. Set up PARDISO control parameter
C..
      do i = 1, 64
         iparm(i) = 0
      end do
      iparm(1) = 1 ! no solver default
      iparm(2) = 2 ! fill-in reordering from METIS
      iparm(3) = 1 ! numbers of processors
      iparm(4) = 0 ! no iterative-direct algorithm
      iparm(5) = 0 ! no user fill-in reducing permutation
      iparm(6) = 0 ! =0 solution on the first n compoments of x
      iparm(7) = 0 ! not in use
      iparm(8) = 9 ! numbers of iterative refinement steps
      iparm(9) = 0 ! not in use
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(12) = 0 ! not in use
      iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
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
      mtype = -2 ! symmetric, indefinite is -2
C.. Initiliaze the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
      do i = 1, 64
         pt(i) = 0
      end do
C.. Reordering and Symbolic Factorization, This step also allocates
C all memory that is necessary for the factorization
      phase = 11 ! only reordering and symbolic factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, Ax, Ap, Ai,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, error)
*      WRITE(*,*) 'Reordering completed ... '
      IF (error .NE. 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
         STOP
      END IF
*      WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
*      WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
C.. Factorization.
      phase = 22 ! only factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, Ax, Ap, Ai,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, error)
*      WRITE(*,*) 'Factorization completed ... '
      IF (error .NE. 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
         STOP
      ENDIF
C.. Back substitution and iterative refinement
      iparm(8) = 2 ! max numbers of iterative refinement steps
      phase = 33 ! only factorization
      do i = 1, n
         b(i) = 1.d0
      end do
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, Ax, Ap, Ai,
     1 idum, nrhs, iparm, msglvl, b, x, error)
*      WRITE(*,*) 'Solve completed ... '
*      WRITE(*,*) 'The solution of the system is '
      WRITE(*,100) (i, x(i), i = 1, n)
 100  FORMAT('x(',i1,') = ',f8.4)
         
C.. Termination and release of memory
      phase = -1 ! release internal memory
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, error)
      END
