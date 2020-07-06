*=======================================================================
* File:    irrand.f
* Author:  Alex Stivala
* Created: November 2009
*
* Function to return pseudo-random integer in [1, N]
*
* $Id: irrand.F 3240 2010-01-18 03:28:54Z alexs $
*=======================================================================

      INTEGER FUNCTION IRRAND(N)
*     .. Scalar Arguments ..
      INTEGER N
*     ..
*
*  Purpose
*  =======
*
*     Return a pseudo-random integer from 1..N
*
*     .. Intrinsic Functions ..
#if defined(__INTEL_COMPILER) || defined(__PORTLAND_COMPILER)
      intrinsic RANDOM_NUMBER,int
#else
      intrinsic rand,int
#endif
*     ..
*
      DOUBLE PRECISION X
#if defined(__INTEL_COMPILER) || defined(__PORTLAND_COMPILER)
      CALL RANDOM_NUMBER(X)
      IRRAND = INT(X * N) + 1
#else
      IRRAND = INT(RAND(0) * N) + 1
#endif
      RETURN
      END
