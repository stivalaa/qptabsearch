* umftst.f
* FORTRAN-77 version of testumfpack.c to test UMFPACK FORTRAN 
* interface.
* The C example program is from the UMFPACK (v5.2.0) User Guide 
* by Timothy A. Davis, University of Florida 
* http://www.cise.ufl.edu/research/sparse/umfpack/
*
* $Id: umftst.f 1690 2008-07-16 04:35:47Z astivala $
      program umftst
      implicit none

      external umf4def,umf4sym,umf4num,umf4solr,umf4fsym,umf4fnum
      integer n,nz
      parameter (n=5, nz=12)

      integer Ap(n+1), Ai(nz), numeric, symbolic, sys, i
      double precision Ax(nz), x(n), b(n), control(20), info(90)

      data Ap /0, 2, 5, 9, 10, 12/
      data Ai /0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4/
      data Ax /2.0d0, 3.0d0, 3.0d0, -1.0d0, 4.0d0, 4.0d0, -3.0d0, 1.0d0,
     $  2.0d0, 2.0d0, 6.0d0, 1.0d0/
      data b /8.0d0, 45.0d0, -3.0d0, 3.0d0, 19.0d0/

*     set default parameters
      call umf4def(control)

*     pre-order and symbolic analysis
      call umf4sym(n, n, Ap, Ai, Ax, symbolic, control, info)
      if (info(1) .lt. 0) then
        write(*,*) 'UMF4SYM failed, info = ', info(1)
        stop
      endif

*     numeric factorization
      call umf4num(Ap, Ai, Ax, symbolic, numeric, control, info)
      if (info(1) .lt. 0) then
        write(*,*) 'UMF4NUM failed, info = ', info(1)
        stop
      endif

*     solve Ax=b, with iterative refinement
      sys = 0
      call umf4solr(sys, Ap, Ai, Ax, x, b, numeric, control, info)
      if (info(1) .lt. 0) then
        write(*,*) 'UMFSOLR failed, info = ', info(1)
        stop
      endif

*     free the numeric factorization
      call umf4fnum(numeric)

      write(*,100) (i, x(i), i = 1, n)
 100  format('x(',i1,') = ',f8.4)

*     free the symbolic analysis
      call umf4fsym(symbolic)

      end

