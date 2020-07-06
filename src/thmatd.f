*=======================================================================
* File:    thmatd.f 
* Author:  Alex Stivala
* Created: November 2009
*
* FORTRAN-77 implementation of tableau matching (discrete)
* using local search heuristics.
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
* Then we use local search to try to improve the score.
* We work backwards from last SSE in query SSE sequence to first.
* For each one, we find the max score for matching to any of the
* other SSEs in other structure (obeying the constraints in use
* ie.. LORDER, LTYPE), and also if LUNMAP is set of mapping the SSE
* to no SSE in other structure (ie removing it from mapping between
* query and db structure). So we make a greedy choice of what to 
* map it to. Note that this algorithm (working from last SSE in
* sequence) is chosen specifically for  LORDER constraint,
* since we only then are allowed to map to SSEs in the other
* sequence in such a way as to have a 'non-crossing' mapping;
* if we started from the beginning then we would most likely
* have few or no choices without crossing mapping of second in
* sequence etc. -- only by starting from end do we open up
* more choices to try.
*
* TODO: always assuming lorder is set at the moment
*
* $Id: thmatd.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================

      subroutine thmatd(tab1,ldt1,n1,tab2,ldt2,n2,
     $     ltype,lorder,dmat1,dmat2,score,ssemap,info)

      implicit none
*
* thmatd - find match score of two tableaux using local search heuristic
*
*   The tableaux are specified here as arrays of 2 char tableau codes
*

*
*     .. Scalar Arguments ..
      integer ldt1,n1,ldt2,n2,info
      logical ltype,lorder
      integer score
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
* ltype  (input) LOGICAL
*        if true, disalow matches between SSEs that are not the same 
*        type. This uses the entry on the main diagonal as the indicator
*        of SSE type, and penalizes matches between entries with
*        different types.
*
* lorder (input) LOGICAL
*        if true, penalize matches between SSEs not maintaining sequence
*        order between the tableaux i.e. if i < k and j >= l for i,k
*        indices in tab1 and j,l indices in tab2.
*
* dmat1  (input) DOUBLE PRECISION array, dimension (n1,n1)
*        SSE distance matrix for one structure. symmetric
*
* dmat2  (input) DOUBLE PRECISION array, dimension (n2,n2)
*        SSE distance matrix for second structure. symmetric
*
* score  (output) INTEGER
*        score of matching tab1 with tab2 using the heuristic methods
*
* ssemap (output) INTEGER vector, dimension (n1)
*        solution SSE map vector of dimension n1. 
*        Each ssemap(i) is the SSE index
*        in tab2 that SSE i in tab1 is matched with.
*
* info   (output) INTEGER
*        on exit, status of the computation
*        =  0 : successful exit
*        =  1 : cannot setup intial ssemap with both lorder and ltype
*
*=======================================================================
*

*     .. Parameters ..

*     Threshold for penalizing difference in distances between SSEs
*     (Angstroms) - if difference in distances between SSEs in the two
*     structures exceeds this then penalize that match
      double precision mxssed
      parameter (mxssed = 4.0d0)
*
*     The bottom score is used for setting 'no match' when e.g. the
*     query structure is larger than the db structure
      integer bottom
      parameter (bottom = -99999)
*
*     limit on number of iterations
      integer maxiter
      parameter (maxiter = 1000)
*     largest tableau (omega matrix) dimension allowed
      integer maxdim
      parameter(maxdim = 110)
*
*     ..
*     .. Local Scalars ..
      integer k,iter,ssei,oldj,newj,oldscore
      integer startj,maxj,maxscore,endj
      logical lunmap
*     ..
*     .. Local Arrays ..
*     reverse ssemap: revmap(j) for j index in tab2 is the index i
*     in tab1 that matches that sse i.e. if ssemap(i) = j then
*     revmap(j) = ssemap(i) and vice versa, for quick lookup of what
*     is matched so we can easily check that one-to-one mapping maintained
*     revmap has dimension (n2)
      integer revmap(maxdim)
*     ..
*     .. Intrinsic Functions ..
      intrinsic min
*     ..
*     .. External Subroutines and Functions ..
      integer  tmscord
      external tmscord
      integer  irrand
      external irrand
      external firsttypeind
      integer  firsttypeind
      external icopy
      external thinit
*     ..
*     .. Executable Statements ..
*     

*     TODO make this an option 
      lunmap = .true.

      if (.not. lunmap .and. (n1 .gt. n2)) then
*        query structure is larger than db structure: no match
         score = bottom
         return
      endif

*     setup initial mapping
      call thinit(tab1, ldt1, n1, tab2, ldt2, n2, ltype, lorder,
     $     lunmap, ssemap, revmap, info)
      if (info .ne. 0) then
*       constraints (order,type,unmap) don't allow any matching
         score = bottom
         return
      endif

C      write(*,*) (ssemap(i),i=1,n1)
C      write(*,*) (revmap(i),i=1,n2)
      
      score = tmscord(tab1, ldt1, n1, tab2, ldt2, n2, 
     $     dmat1, dmat2, ssemap)
      oldscore = score
      iter = 0
      ssei = n1
 100  if (ssei .gt. 0 .and. iter .lt. maxiter) then
*        work backwards from last SSE in query trying to improve score
         startj = ssemap(ssei)
         k = ssei
 102     if (startj .eq. 0 .and. k .gt. 0) then
            startj = ssemap(k)
            k = k - 1
            goto 102
         endif
         if (startj .eq. 0) then
            startj = n2+1
         endif

*        'hill climbing': find max score for sse in tab2 (db structure)
*        TODO at the moment, ONLY going forwards i.e. assuming LORDER
         oldscore = score
         maxscore = score
         maxj = startj
         if (ssei .eq. n1) then
            endj = n2+1
         elseif (ssemap(ssei+1) .eq. 0) then
            endj = 0
            k = 1
 105        if (endj .eq. 0 .and. ssei + k .le. n1) then
               endj = ssemap(ssei + k)
               k = k + 1
               goto 105
            endif
         else
            endj = ssemap(ssei+1)
         endif
         newj = firsttypeind(tab2, ldt2, n2, startj, 
     $        tab1(ssei,ssei), revmap, endj)
C         write (*,*) (ssemap(k), k= 1, n1)
C         write (*,*) iter, ssei, startj, endj, newj
         oldj = ssemap(ssei)
 110     if (newj .ne. 0) then
            ssemap(ssei) = newj
            revmap(oldj) = 0
            revmap(newj) = ssei
            score = tmscord(tab1, ldt1, n1, tab2, ldt2, n2,
     $           dmat1, dmat2, ssemap)
            if (score .gt. maxscore) then
               maxscore = score
               maxj = newj
C               write(*,*) 'x',ssei,maxj
            endif
            ssemap(ssei) = oldj
            revmap(oldj) = ssei
            revmap(newj) = 0
            if (ssei .eq. n1) then
               endj = n2+1
            elseif (ssemap(ssei+1) .eq. 0) then
               endj = 0
               k = 1
 115           if (endj .eq. 0 .and. ssei + k .le. n1) then
                  endj = ssemap(ssei + k)
                  k = k + 1
                  goto 115
               endif
            else
               endj = ssemap(ssei+1)
            endif
            newj = firsttypeind(tab2, ldt2, n2, newj+1, 
     $           tab1(ssei,ssei), revmap, endj)
C            write (*,*) (ssemap(k), k= 1, n1)
C            write (*,*) iter, ssei, startj, endj, newj
            goto 110
*         end of while (newj .ne. 0)
         endif

         if (lunmap) then
*        also try removing the SSE from matching
            revmap(ssemap(ssei)) = 0
            ssemap(ssei) = 0
            revmap(oldj) = 0
            score = tmscord(tab1, ldt1, n1, tab2, ldt2, n2,
     $           dmat1, dmat2, ssemap)
            if (score .gt. maxscore) then
               maxscore = score
               maxj = 0
            endif
            ssemap(ssei) = oldj
            revmap(oldj) = ssei
         endif

         if (maxscore .gt. oldscore) then
            if (ssemap(ssei) .gt. 0) then
               revmap(ssemap(ssei)) = 0
            endif
            ssemap(ssei) = maxj
            if (maxj .gt. 0) then
               revmap(maxj) = ssei
            endif
            score = maxscore
         endif
         ssei = ssei - 1
         iter = iter + 1
         goto 100
*     end of while (iter .lt. maxiter) 
      endif
      
      info = 0
      return
      end
