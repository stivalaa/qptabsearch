*=======================================================================
* File:    icopy.f
* Author:  Alex Stivala
* Created: November 2009
*
* Copy integer vector y <- x
*
* $Id: icopy.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================
      subroutine icopy(n, x, y)
      implicit none
*
* icopy - copy integer vector: y <- x
*     
      integer       n
      integer       x(*)
      integer       y(*)

* Arguments
* =========
*
*
* n      (input) INTEGER
*        Length of vector to copy
*
* x      (input) INTEGER vector, dimension(n)
*        source vector
* 
* y      (output) INTEGER vector, dimension(n)
*        destination vector
* 
      integer i
      do 2000 i = 1, n
         y(i) = x(i)
 2000 continue
      end

