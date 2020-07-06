*=======================================================================
* File:    tscord.f
* Author:  Alex Stivala
* Created: July 2008
*
* Tableau (discrete) element matching score function
* 
* $Id: tscord.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================

      integer function tscord(x, y)

      implicit none
*
* tscord - Tableau (discrete) matching score function
*
*    Return the tableau matching score between two tableau entries
*    x and y.
*    The scores are all postiive to help ensure integer solutions
*    to relaxed QP.
*    The score is 2 if the tableau entries are equal, 1 if they are
*    equal in only one position, else 0.
*
*     .. Scalar Arguments ..
      character  x(*),y(*)
*     ..
*
* Arguments
* =========
*
* x (input) CHARACTER*2
*     two character tableau code
*
* y (input) CHARACTER*2
*     two character tableau code
*
*=======================================================================
*

*     ..
*     .. Executable Statements ..
*     
      if (x(1) .eq. y(1)) then
        if (x(2) .eq. y(2)) then
            tscord = 2
        else
            tscord = 1
        endif
      elseif (x(2) .eq. y(2)) then
        tscord = 1
      else
        tscord = -2
      endif
      return
      end 

