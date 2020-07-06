*=======================================================================
* File:    tstsp1.f
* Author:  Alex Stivala
* Created: July 2008
*
* testp1 - test harness for spphs1
*
* $Id: patstsp1.f 2072 2009-02-26 23:55:28Z astivala $
*=======================================================================

      program testp1

      implicit none

      integer           m,n,j,incx,incy,lda,info,lde,ldf
      parameter         (m = 6, n = 15)
      integer           nzmax
      parameter         (nzmax=100)
      double precision  z,ob,blasa,blasb,alpha
      double precision  A(m, n),x(n+1),y(m),b(m),onevec(n)
      double precision  avec(m)
      character         trans

*     workspace arrays for SPPHS1
      integer           Ep(m+n+1),Ei(nzmax)
      double precision  Ex(nzmax)
      double precision  y1(m+n),y2(m+n)
      double precision dx(n),F(m+n, 2)
      double precision  u(n+2), v(n+2), umzv(n+2), obz(n+2)

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
      call paphs1(A, lda, m, n, avec, b, x, y, z, ob, alpha, 
     $     dx, Ep, Ei, Ex, nzmax, F, ldf, u, v, umzv, obz,
     $     y1, y2, info)


      write(*,290),info
 290  format ('info = ',i4)
      write(*,300),z,ob
 300  format ('z = ',f8.4,' ob = ',f8.4)
      write(*,*),'x = '
      write(*,410), (x(j), j = 1, n)
 410  format (f8.4)
      end

