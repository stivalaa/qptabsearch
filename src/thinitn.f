*=======================================================================
* File:    thinitn.f 
* Author:  Alex Stivala
* Created: November 2009
*
* Build the initial mapping of the two structurs for heruristic
* tableaux matching algoriths. Numeric (omega matrix) version.
*
* Here we (unlike QP tableau search (except the experimental versions
* where constraints used to enforce all SSEs in first (query) structure
* matched)) create an asymmetry between
* the first (query) and second (database) structures. We force all
* SSEs in the query structure to match to something --- if it is
* larger than the dbase structure then we return a 'bottom' score (no
* match). (We only do this if LUNMAP is not set so unmapping of a query
* SSE is not an allowed move).
*
* Otherwise we make an initial matching where we just go along
* the sequence set match of same SSEs e.g. if 1st in query is helix,
* match that to first helix in db struture, and so on.
* (Unless LTYPE flag not set, then we don't care about SSE types and 
* just go along sequence of SSEs).
* Then compute the score.
*
* $Id: thinitn.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================

      subroutine thinitn(omega1,ldt1,n1,omega2,ldt2,n2,
     $     ltype,lorder,lunmap,ssemap,revmap,info)

      implicit none
*
* thinitn - initial mapping between two structrs for herurist methods
*
*   The tableaux are specified here as omega matrices
*

*
*     .. Scalar Arguments ..
      integer ldt1,n1,ldt2,n2,info
      logical ltype,lorder,lunmap
*     ..
*     .. Array Arguments ..
      double precision  omega1(ldt1,*),omega2(ldt2,*)
      integer ssemap(*)
      integer revmap(*)
*     ..
*

*
* Arguments
* =========
*
* omega1  (input) DOUBLE PRECISION array, dimension (n1,n1)
*        Omega matrix for one structure. Symmetric.
*
* ldt1   (input) INTEGER
*        Leading dimension of omega1 array. ldt1 >= n1
*
* n1     (input) INTEGER
*        Dimension of omega1 array.
*
* omega2  (input) DOUBLE PRECISION array, dimension (n2,n2)
*        Omega matrix for second structure. Symmetric.
*
* ldt2   (input) INTEGER
*        Leading dimension of omega2 array. ldt2 >= n2
*
* n2     (input) INTEGER
*        Dimension of omega2 matrix.
*
* ltype  (input) LOGICAL
*        if true, disalow matches between SSEs that are not the same 
*        type. This uses the entry on the main diagonal as the indicator
*        of SSE type, and penalizes matches between entries with
*        different types.
*
* lorder (input) LOGICAL
*        if true, penalize matches between SSEs not maintaining sequence
*        order between the tableaux i.e. if i < k and j >= l for i,k
*        indices in omega1 and j,l indices in omega2.
*
* lunmap (input) LOGICAL
*        if true, the query struct (omega1) is allowed to have one or
*        more SSEs not mapped (0 in the ssemap). Otherwise, all of
*        its SSEs must be mapped and so omega1 must be smaller than omega2.
*        
* ssemap (output) INTEGER vector, dimension (n1)
*        solution SSE map vector of dimension n1. 
*        Each ssemap(i) is the SSE index
*        in omega2 that SSE i in omega1 is matched with.
*
* revmap (output) INTEGER vector, dimension(n2)
*     reverse ssemap: revmap(j) for j index in omega2 is the index i
*     in omega1 that matches that sse i.e. if ssemap(i) = j then
*     revmap(j) = ssemap(i) and vice versa, for quick lookup of what
*     is matched so we can easily check that one-to-one mapping maintained
*
* info   (output) INTEGER
*        on exit, status of the computation
*        =  0 : successful exit
*        =  1 : cannot setup intial ssemap with both lorder and ltype
*
*=======================================================================
*

*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      integer i,j
*     ..
*     .. Local Arrays ..
*     ..
*     .. Intrinsic Functions ..
      intrinsic min
*     ..
*     .. External Subroutines and Functions ..
*     ..
*     .. Executable Statements ..
*     
*     initialize ssemap to all zero meaning no match for each sse
      do 2 i = 1, n1
         ssemap(i) = 0
 2    continue
      do 3 j = 1, n2
         revmap(j) = 0
 3    continue

      if (.not. lunmap .and. (n1 .gt. n2)) then
         write(*,*) 'ERROR: LUNMAP NOT SET BUT N1 GT N2'
         return
      endif      

*     Initial SSE map set by matching along sequence, only matching
*     SSEs of same type if LTYPE flag is set.
*     TODO: maybe do something else if LORDER flag is not set
      if (ltype) then
         j = 1
         do 10 i = 1, n1
 5          if (j .lt. n2 .and. (omega1(i,i) .ne. omega2(j,j))) then
               j = j + 1
               goto 5
            endif
            if (j .ge. n2) then
C               write (*,*) 'CANNOT INITIALIZE WITH ORDER CONSTRAINT'
*              TODO: do something about this - if not using LORDER
*               then could assign in some other order
               if (.not. lunmap) then
*                 couldn't map all SSEs in omega1, but were required to
                  info = 1
               else
*                 not all SSEs in omega1 are mapped, but that's OK
                  info = 0
               endif
               return
            else
               ssemap(i) = j
               revmap(j) = i
               j = j + 1
            endif
 10      continue
      else
         do 20 i = 1, min(n1,n2)
            ssemap(i) = i
            revmap(i) = i
 20      continue
      endif
      info = 0
      return
      end
