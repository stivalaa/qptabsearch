*=======================================================================
* File:    tmtest3.f
* Author:  Alex Stivala
* Created: July 2008
*
* Test program for tmatn.f
*
* $Id: tmtss3.f 2968 2009-11-22 01:21:50Z astivala $
*=======================================================================

      program tmtest3

      implicit none
      
      external tmatn

      integer n1,n2,info
      parameter(n1 = 7, n2 = 6)
      double precision omega1(n1,n1),omega2(n2,n2)
      double precision dmat1(n1,n1),dmat2(n2,n2)
      double precision score
      double precision soln(n1*n2+n1+n2)

*     MATLAB GetTableau('1ABA')
      data omega1  / 0d0,  -2.8600d0,  -0.3600d0,  -2.2100d0,
     $     2.9200d0,  -0.8100d0,  -2.7700d0,
     $     -2.8600d0,   1.0000d0,  -2.9100d0,   1.1200d0,
     $     -0.4500d0,   2.1300d0,   0.3300,
     $     -0.3600d0,  -2.9100d0,        0d0,  -2.1800d0,   2.5600d0,
     $     -1.1600d0,  -0.1200d0,
     $     -2.2100d0,   1.1200d0,  -2.1800d0,
     $     1.0000d0,   1.0200d0,  -1.8800d0,   0.8400d0,
     $     2.9200d0,  -0.4500d0,   2.5600d0, 
     $     1.0200d0,        0d0,   2.5600d0,   0.6000d0,
     $     -0.8100d0,   2.1300d0,  -1.1600d0,
     $     -1.8800d0,   2.5600d0,        0d0,   1.9700d0,
     $     -2.7700d0,   0.3300d0,  -0.1200d0,
     $     0.8400d0,   0.6000d0,   1.9700d0,   1.0000d0/


*     MATLAB GetTableau('1UBQ')
      data omega2 / 0d0, 2.6500d0,  -1.1700d0, 2.0600d0, -1.2600d0, 
     $     -0.5900d0,
     $     2.6500d0,    0d0, 2.1400d0,  -1.4300d0, 1.5700d0, 2.1000d0,
     $     -1.1700d0, 2.1400d0, 1.0000d0, 2.0900d0,  -1.0400d0,
     $     -1.2200d0,
     $     2.0600d0,  -1.4300d0, 2.0900d0,    0d0, 2.9600d0, 2.5800d0,
     $     -1.2600d0, 1.5700d0,  -1.0400d0, 2.9600d0,   
     $     0d0,  -0.7400d0,
     $     -0.5900d0, 2.1000d0,  -1.2200d0, 2.5800d0,  -0.7400d0,  0d0/

      data dmat1 / 0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0 /

      data dmat2 / 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0 /

*     should get score approx -4 (MATLAB TableauMatch -4.3232 or -4.5832)
*     -3.5124 when not using diagnoal elements

      call tmatn(omega1, n1, n1, omega2, n2, n2, .false., .false.,
     $           dmat1, dmat2,
     $           score, soln, info)
      if (info .ne. 0) then
         write(*,*),'TABMATN FAILED'
      else
         write(*,900),score
 900     format('score = ',f6.2)
      endif
      end

