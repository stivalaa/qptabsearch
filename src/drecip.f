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
*     $Id: drecip.F 2969 2009-11-22 01:29:39Z astivala $
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

