      DOUBLE PRECISION FUNCTION norm2(score,size1,size2)
*     .. Scalar Arguments ..
      DOUBLE PRECISION score
      INTEGER size1
      INTEGER size2
*     ..
*     .. Array Arguments ..
*     ..
*
*  Purpose
*  =======
*
*     Compute size normalization (norm2 of Pelta et al. (2008)) of
*     tableau matching score
*    
*     norm2(struct1, struct2) = 2*score / (#sses(struct1) + #sses(struct2))
*
*   Parmeters:
*        score - tableau matching score
*        size1 - number of SSEs in struccture 1
*        size2 - number of SSEs in structure 2
*
*  Return value: norm2 as above.
*
*     $Id: drecip.F 2969 2009-11-22 01:29:39Z astivala $
*
      intrinsic dble
      norm2 = 2.0d0 * score / dble(size1 + size2)
      END
