*=======================================================================
* File:    gumbel.F
* Author:  Alex Stivala
* Created: July 2010
*
* Functions to compute z-score and p-value from tableau matching score,
* according to Gumbel distribution parameters.
* See thesis Ch. 6 (s6.2.2) and Ortiz et al (2002), 
* Abagyan & Batalov (1997), Versetr0m & Taylor (2006),
* Kolbeck et al (2006), Levitt & Gerstein (1998).
*
* This file contains double precision functions zgumbel and pvgumbel
* 
* $Id: rdtabd.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================


      double precision function pvgumbel(z)

      implicit none
*
* pvgumbel - compute p-value from Gumbel distribution
*
*     .. Scalar Arguments ..
      double precision z
*
*
* Arguments
* =========
*
* z (input) DOUBLE PRECISION
*        Gumbel z-score from pvgumbel to compute pvalue for
*
* Return value
* ============
*     p-value computed for x according to Gumbel(a,b) distribution
*    
*=======================================================================
*

*     .. Parameters ..
      double precision pi
      parameter(pi=3.141592653589793D0)
*     Euler-Mascheroni constant
      double precision egamma
      parameter(egamma = 0.5772156649015328606d0)
*     ..
*     .. Local Scalars ..
*     ..
*     .. Local Arrays ..
*     ..
*     .. Intrinsic Functions ..
      intrinsic sqrt,exp
*     ..
*     .. External Subroutines and Functions ..
*     ..
*     .. Executable Statements ..
*     
      pvgumbel = 1.0d0-exp(-1.0d0*exp(-1.0d0*(z*pi/sqrt(6.0d0)+egamma)))
      end

*=======================================================================



*=======================================================================

      double precision function zgumbel(x)

      implicit none
*
* zgumbel - compute Z-score from Gumbel distribution
*
*
*  The Gumbel distribution parameters a and b are estimated by MLE
*  for structures in different folds (see thesis ch.6 (s6.2.2)
*  for query200 in ASTRAL SCOP 1.75 95% nr.
*
*     .. Scalar Arguments ..
      double precision x
*
*
* Arguments
* =========
*
* x  (input) DOUBLE PRECISION
*        Tableau matching score to compute z-score for
*
* Return value
* ============
*     Z-score computed for x according to Gumbel(a,b) distribution
*    
*=======================================================================
*

*     .. Parameters ..
      double precision pi
      parameter(pi=3.141592653589793D0)
*     Euler-Mascheroni constant
      double precision egamma
      parameter(egamma = 0.5772156649015328606d0)
*     Gumbel distribution location parameter
      double precision a
      parameter(a = 0.4731609844614057d0)
*     Gumbel distribution scale parameter
      double precision b
      parameter(b = 0.44602089192192d0)
*     ..
*     .. Local Scalars ..
      double precision mu
      double precision sigma
*     ..
*     .. Local Arrays ..
*     ..
*     .. Intrinsic Functions ..
      intrinsic sqrt
*     ..
*     .. External Subroutines and Functions ..
*     ..
*     .. Executable Statements ..
*     
      mu = a + b * egamma
      sigma = (pi / sqrt(6.0d0)) * b
      zgumbel  = (x - mu) / sigma
      end

