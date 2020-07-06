*=======================================================================
* File:    matstsp1.f
* Author:  Alex Stivala
* Created: July 2008
*
* matstsp1 - test harness for maphs1
*
* $Id: matstsp1.f 3026 2009-12-10 00:57:36Z alexs $
*=======================================================================

      program testp1

      implicit none

      integer           m,n,j,incx,incy,lda,info,lde,ldf
      parameter         (m = 6, n = 15)
      integer           nzmax,lkeep
      parameter         (nzmax=100)
      parameter         (lkeep = 5 *(m+n) + nzmax+nzmax+42)
      double precision  z,ob,blasa,blasb,alpha
      double precision  A(m, n),x(n+1),y(m),b(m),onevec(n)
      double precision  avec(m)
      character         trans

*     workspace arrays for SPPHS1
      integer           Eirn(nzmax),Ejcn(nzmax)
      double precision  Ea(nzmax)
      double precision  y1(m+n),y2(m+n)
      double precision dx(n),F(m+n, 2)
      double precision  u(n+2), v(n+2), umzv(n+2), obz(n+2)
      integer           keep(lkeep)

      external paphs1


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
      call maphs1(A, lda, m, n, avec, b, x, y, z, ob, alpha, 
     $     dx, Eirn, Ejcn, Ea, lkeep, keep, nzmax, F, ldf,
     $     u, v, umzv, obz,
     $     y1, y2, info)


      write(*,290),info
 290  format ('info = ',i4)
      write(*,300),z,ob
 300  format ('z = ',f8.4,' ob = ',f8.4)
      write(*,*),'x = '
      write(*,410), (x(j), j = 1, n)
 410  format (f8.4)
      end

