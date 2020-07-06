*=======================================================================
* File:    tmtest6.f
* Author:  Alex Stivala
* Created: July 2008
*
* Test program for tmatn.f
*
* example for trying out TABMATN, like the example
* in Konagurthu et al 2008 of d1xdtt3 and d1hhoa_
* (Diptheria toxin and human alpha-haemoglobin)
*
*
* $Id: tmtst6.f 2968 2009-11-22 01:21:50Z astivala $
*=======================================================================

      program tmtest6

      implicit none
      
      external tmatn

      integer n1,n2,info
      parameter(n1 = 6, n2 = 7)
      double precision omega1(n1,n1),omega2(n2,n2)
      double precision dmat1(n1,n1),dmat2(n2,n2)
      double precision score
      double precision soln(n1*n2+n1+n2)


*     MATLAB GetTableau('1UBQ')
      data omega1 / 0 , 2.6500 ,  -1.1700 , 2.0600 , -1.2600 ,  -0.5900,
     $     2.6500 ,    0 , 2.1400 ,  -1.4300 , 1.5700 , 2.1000,
     $     -1.1700 , 2.1400 , 1.0000 , 2.0900 ,  -1.0400 ,  -1.2200,
     $     2.0600 ,  -1.4300 , 2.0900 ,    0 , 2.9600 , 2.5800,
     $     -1.2600 , 1.5700 ,  -1.0400 , 2.9600 ,    0 ,  -0.7400,
     $     -0.5900 , 2.1000 ,  -1.2200 , 2.5800 ,  -0.7400 ,    0 /


* pytableaucreate.py -t dssp -n /var/tmp/d1hhoa_.ent
      data omega2 / 0.00,  1.75, -2.51,  0.93,  1.46,  1.07, -1.67, 
     $     1.75,  0.00,  2.01, -1.06, -0.33,  2.18, -1.33, 
     $     -2.51, 2.01,  0.00,  2.66, -2.32, -1.59,  1.72, 
     $     0.93, -1.06,  2.66,  0.00,  0.73,  1.98, -0.94, 
     $     1.46, -0.33, -2.32,  0.73,  0.00,  2.15, -1.20, 
     $     1.07,  2.18, -1.59,  1.98,  2.15,  0.00,  2.60, 
     $     -1.67, -1.33,  1.72, -0.94, -1.20, 2.60,  0.00 /

      data dmat1 / 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0 /

      data dmat2 / 0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0,
     $             0, 0, 0, 0, 0, 0, 0 /

      call tmatn(omega1, n1, n1, omega2, n2, n2, .false., .false.,
     $           dmat1, dmat2,
     $           score, soln, info)

*     MATLAB TAbleauMatch gives score -8.274

      if (info .ne. 0) then
         write(*,*),'TABMATN FAILED'
      else
         write(*,900),score
 900     format('score = ',f6.2)
      endif
      end
