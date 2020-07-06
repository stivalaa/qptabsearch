*=======================================================================
* File:    tmtest2.f
* Author:  Alex Stivala
* Created: July 2008
*
* Test program for tmatn.f
*
* $Id: tmtss2.f 2968 2009-11-22 01:21:50Z astivala $
*=======================================================================

      program tmtest2

      implicit none
      
      external tmatn

      integer n1,n2,info
      parameter(n1 = 7, n2 = 7)
      double precision omega1(n1,n1),omega2(n2,n2)
      double precision dmat1(n1,n1),dmat2(n2,n2)
      double precision score
      double precision soln(n1*n2+n1+n2)


* pytableaucreate.py -t dssp -n /var/tmp/d1hhoa_.ent
      data omega1 / 0.00,  1.75, -2.51,  0.93,  1.46,  1.07, -1.67, 
     $     1.75,  0.00,  2.01, -1.06, -0.33,  2.18, -1.33, 
     $     -2.51, 2.01,  0.00,  2.66, -2.32, -1.59,  1.72, 
     $     0.93, -1.06,  2.66,  0.00,  0.73,  1.98, -0.94, 
     $     1.46, -0.33, -2.32,  0.73,  0.00,  2.15, -1.20, 
     $     1.07,  2.18, -1.59,  1.98,  2.15,  0.00,  2.60, 
     $     -1.67, -1.33,  1.72, -0.94, -1.20, 2.60,  0.00 /

* pytableaucreate.py -t dssp -n /var/tmp/d1hhoa_.ent
      data omega2 / 0.00,  1.75, -2.51,  0.93,  1.46,  1.07, -1.67, 
     $     1.75,  0.00,  2.01, -1.06, -0.33,  2.18, -1.33, 
     $     -2.51, 2.01,  0.00,  2.66, -2.32, -1.59,  1.72, 
     $     0.93, -1.06,  2.66,  0.00,  0.73,  1.98, -0.94, 
     $     1.46, -0.33, -2.32,  0.73,  0.00,  2.15, -1.20, 
     $     1.07,  2.18, -1.59,  1.98,  2.15,  0.00,  2.60, 
     $     -1.67, -1.33,  1.72, -0.94, -1.20, 2.60,  0.00 /


      data dmat1 / 0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0 /

      data dmat2 / 0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0 /

      call tmatn(omega1, n1,n1, omega2, n2 , n2, .false., .false.,
     $           dmat1, dmat2,
     $           score, soln, info)

*     score from MATLAB is -16.4933

      if (info .ne. 0) then
         write(*,*),'TABMATN FAILED'
      else
         write(*,900),score
 900     format('score = ',f6.2)
      endif
      end

