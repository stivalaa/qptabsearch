      SUBROUTINE DEMVV(N,X,Y)
*     .. Scalar Arguments ..
      INTEGER N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*),Y(*)
*     ..
*
*  Purpose
*  =======
**
*     inplace elementwise vector-vector multiply.
*
*     Y <- X .* Y
*
*     BLAS seems to have no subroutine for this
*     so this is a trivial implementation
*     (not using BLAS style loop unrolling etc.,
*     no increment, etc.)
*     just to keep things simple.
*     An alternative implementation would be to
*     put one vector on main diagonal of a matrix
*     and use Level 2 BLAS subroutine DGBMV 
*     Alex Stivala 04July2008.
*     $Id: demvv.F 2969 2009-11-22 01:29:39Z astivala $
*
* Arguments
* =========
*
*  X      - DOUBLE PRECISION vector of dimension N
*           input vector X
*
*  Y      - DOUBLE PRECISION vector of dimension N
*           on input, values of vector Y
*           on exit each element Y(i)
*           is set to X(i)*y(i)
*           where y(i) is input Y value.
*
*
*     .. Local Scalars ..
      INTEGER I
*     ..
*
      DO 10 I = 1,N
          Y(I) = X(I) * Y(I)
   10 CONTINUE
      RETURN
      END

