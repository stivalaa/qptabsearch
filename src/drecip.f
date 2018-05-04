      SUBROUTINE DRECIP(N,DX)
*     .. Scalar Arguments ..
      INTEGER N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
**
*     inplace elementwise reciprocal of vector
*     BLAS seems to have no subroutine for this
*     so this is a trivial implementation
*     (not using BLAS style loop unrolling etc.,
*     no increment, etc.)
*     just to keep things simple.
*     Alex Stivala 04July2008.
*     $Id: drecip.f 1565 2008-07-05 03:32:54Z astivala $
*
*
*     .. Local Scalars ..
      INTEGER I
*     ..
*
      DO 10 I = 1,N
          DX(I) = 1/DX(I)
   10 CONTINUE
      RETURN
      END

