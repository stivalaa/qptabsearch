*=======================================================================
* File:    tscorn.f
* Author:  Alex Stivala
* Created: July 2008
*
* Numeric tableau (Omega matrix) element matching score function
* 
* $Id: tscorn.f 1700 2008-07-17 09:05:38Z astivala $
*=======================================================================

      double precision function tscorn(x, y)

      implicit none
*
* tscorn - Numeric tableau (Omega matrix) matching score function
*
*    Return the tableau matching score between two Omega matrix entries
*    x and y as per Kamat et al (2008)
*
*     .. Scalar Arguments ..
      double precision x, y
*     ..
*
* Arguments
* =========
*
* x (input) DOUBLE PRECISION 
*     angle in (-pi, pi)
*
* y (input) DOUBLE PRECISION 
*     angle in (-pi, pi)
*
*=======================================================================
*

*     .. Parameters ..

      double precision pi,twopi,pion4
      parameter (pi =  3.14159265358979323846d0,
     $     twopi = 2.0d0 * pi,
     $     pion4 = pi / 4.0d0)
*     ..
*     .. Local Scalars ..
      double precision delta
*     ..
*     .. Local Arrays ..
*     ..
*     .. Intrinsic Functions ..
      intrinsic min,abs
*     ..
*     .. External Subroutines and Functions ..
*     ..
*     .. Executable Statements ..
*     
      delta = min( abs(x - y), twopi - abs(x - y) )
      tscorn = pion4 - delta
      return
      end
