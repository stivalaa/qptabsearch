*=======================================================================
* File:    testp1.f
* Author:  Alex Stivala
* Created: July 2008
*
* testp1 - test harness for dlphs1
*
* $Id: testp1.f 1689 2008-07-16 04:21:41Z astivala $
*=======================================================================

      program testp1

      implicit none

      integer           m,n,j,incx,incy,lda,info,lde,ldf,lwork
      parameter         (m = 6, n = 15)
      double precision  z,ob,blasa,blasb,alpha
      double precision  A(m, n),x(n+1),y(m),b(m),onevec(n)
      double precision  avec(m)
      character         trans

*     workspace arrays for DLPHS1
      double precision dx(n),E(m+n, m+n),F(m+n, 2)
      integer ipiv(m+n)
      double precision  u(n+2), v(n+2), umzv(n+2), obz(n+2),work(m+n)

      external dlphs1


*     BLAS level 1 routines
      external dcopy
*     BLAS level 2 routines
      external dgemv


      data A /1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
     $        1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
     $        1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
     $        0.0, 1.0, 0.0, 1.0, 0.0, 0.0,
     $        0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
     $        0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
     $        0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
     $        0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
     $        0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
     $        1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     $        0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
     $        0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
     $        0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
     $        0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
     $        0.0, 0.0, 0.0, 0.0, 0.0, 1.0/

      data b /6*1.0/

      data x /16*1.0/
      ob = x(n+1)
      z = 0.0

*     length of work vector for DSYSV >= (nmax+mmax)*NB where NB is
*     optimal blocksize for DSYSTRF. NB we are using NB=1, this
*     was optimal on x86 Linux anyway.
      lwork = m+n

*     avec = b - A*ones(n,1)
      data onevec /n*1.0/
      trans = 'N'
      blasa = -1.0
      blasb = 1.0
      incx = 1
      lda = m
      incy = 1
      call dcopy(m, b, incx, avec, incy)
      call dgemv(trans, m, n, blasa, A, lda, onevec, incx, blasb,
     $           avec, incy)
      lde = m+n
      ldf = m+n
      call dlphs1(A, lda, m, n, avec, b, x, y, z, ob, alpha, 
     $     dx, E, lde, F, ldf, ipiv, u, v, umzv, obz, 
     $     work, lwork,info)


      write(*,290),info
 290  format ('info = ',i4)
      write(*,300),z,ob
 300  format ('z = ',f8.4,' ob = ',f8.4)
      write(*,*),'x = '
      write(*,410), (x(j), j = 1, n)
 410  format (f8.4)
      end

