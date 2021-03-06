*=======================================================================
* File:    tmscord.f 
* Author:  Alex Stivala
* Created: November 2009
*
* Compute the score for a given SSE matching between two structures
* given their tableaux (discrete version), and distnace matrices.
*
* The score computed is
*
* \sum{i=1,j=1}^{N_A} \sum{j=1,k=1}^{N_B} \zeta(T_{ik},T_{kl}) x_{ik}x{jl}
*
* in the QIP formulation where x_{ij} is the binary indicator variable
* indication SSE i in A matched with SSE j in B. 
*
* But actually here we are representing the matching with the ssemap
* vector so can much more efficiently compute this in only
* O(N_A^2) with 2 nested loops over the ssemap vector rather than 
* requring O(N_A^2 N_B^2) with 4 nested loops in the naive implentation
* of the score computation using indicator variables (required only
* for using a general purpose QP solver, can do it more efficiently here).
* 
* Furthermore, we can actually halve the computation since the tableaux
* matrices are symmetric by only iterating from k = i .. N_A 
* inside the outer loop i = 1 .. N_A.
*
* $Id: tmscord.F 3265 2010-01-25 02:48:12Z alexs $
*=======================================================================

      integer function tmscord(tab1,ldt1,n1,tab2,ldt2,n2,
     $     dmat1, dmat2,
     $     ssemap)

      implicit none
*
* tmscord - Return the score for given SSE matching between two tableaux
*
*   The tableaux are specified here as arrays of 2 char tableau codes
*

*
*     .. Scalar Arguments ..
      integer ldt1,n1,ldt2,n2
*     ..
*     .. Array Arguments ..
      character*2   tab1(ldt1,*),tab2(ldt2,*)
      double precision dmat1(ldt1,*),dmat2(ldt2,*)
      integer ssemap(*)
*     ..
*

*
* Arguments
* =========
*
* tab1  (input) CHARACHTER*2 array, dimension (n1,n1)
*        Tableau for one structure. Symmetric.
*
* ldt1   (input) INTEGER
*        Leading dimension of tab1 array. ldt1 >= n1
*
* n1     (input) INTEGER
*        Dimension of tab1 array.
*
* tab2  (input) CHARACTER*2 array, dimension (n2,n2)
*        Tableau for second structure. Symmetric.
*
* ldt2   (input) INTEGER
*        Leading dimension of tab2 array. ldt2 >= n2
*
* n2     (input) INTEGER
*        Dimension of tab2 matrix.
*
* dmat1  (input) DOUBLE PRECISION array, dimension (n1,n1)
*        SSE distance matrix for one structure. symmetric
*
* dmat2  (input) DOUBLE PRECISION array, dimension (n2,n2)
*        SSE distance matrix for second structure. symmetric
*
* ssemap (input) INTEGER vector, dimension(n1)
*        SSE map vector of dimension n1. Each ssemap(i) is the SSE index
*        in tab2 that SSE i in tab1 is matched with.
*
*=======================================================================
*

*     .. Parameters ..
*     Threshold for penalizing difference in distances between SSEs
*     (Angstroms) - if difference in distances between SSEs in the two
*     structures exceeds this then penalize that match
      double precision mxssed
      parameter (mxssed = 4.0d0)
*     ..
*     .. Local Scalars ..
      integer i,j,k,l
      character*2 tc1,tc2
      integer score
*     ..
*     .. Local Arrays ..
*     ..
*     .. Intrinsic Functions ..
      intrinsic abs
*     ..
*     .. External Subroutines and Functions ..
      integer  tscord
      external tscord

*     ..
*     .. Executable Statements ..
      score = 0
      do 200 i = 1, n1
         do 180 k = i + 1, n1
            j = ssemap(i)
            l = ssemap(k)
*           only add to score if both mapped to something
            if (j .gt. 0 .and. l .gt. 0) then
*              diagonal entries are SSE type not angle
               if (i .ne. k .and. j .ne. l) then
*                 don't add score where difference between SSE
*                 distances exceeds threshold
                  if (abs(dmat1(i,k)-dmat2(j,l)) .le. mxssed) then
                     tc1 = tab1(i,k)
                     tc2 = tab2(j,l)
                     score = score + tscord(tc1, tc2)
                  endif
               endif
            endif
 180     continue
 200  continue
      tmscord = score
      return
      end

