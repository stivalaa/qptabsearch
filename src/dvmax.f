      DOUBLE PRECISION FUNCTION DVMAX(N,DX)
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
*     Return the maximum element in a vector.
*     FORTRAN-90 has an intrinsic for this,
*     but not FORTRAN-77 (and not in BLAS).
*     Alex Stivala 05July2008.
*     $Id: dvmax.f 2056 2009-01-21 23:54:02Z astivala $
*
*     .. Parameters ..
      DOUBLE PRECISION NEGINF
      PARAMETER (NEGINF=-1.79769d+308)
*     .. Local Scalars ..
      INTEGER I
      DOUBLE PRECISION MAXVAL
*     ..
*
      MAXVAL = NEGINF
      DO 10 I = 1,N
         IF (DX(I) .GT. MAXVAL) THEN
            MAXVAL = DX(I)
         ENDIF
   10 CONTINUE
      DVMAX = MAXVAL
      RETURN
      END
