      DOUBLE PRECISION FUNCTION DVMIN(N,DX)
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
*     Return the minimum element in a vector.
*     FORTRAN-90 has an intrinsic for this,
*     but not FORTRAN-77 (and not in BLAS).
*     Alex Stivala 05July2008.
*     $Id: dvmin.F 2969 2009-11-22 01:29:39Z astivala $
*
*     .. Parameters ..
      DOUBLE PRECISION INFTY
      PARAMETER (INFTY=1.79769d+308)
*     .. Local Scalars ..
      INTEGER I
      DOUBLE PRECISION MINVAL
*     ..
*
      MINVAL = INFTY
      DO 10 I = 1,N
         IF (DX(I) .LT. MINVAL) THEN
            MINVAL = DX(I)
         ENDIF
   10 CONTINUE
      DVMIN = MINVAL
      RETURN
      END
