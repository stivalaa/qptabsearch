*=======================================================================
* File:    tsamtn.f 
* Author:  Alex Stivala
* Created: November 2009
*
* FORTRAN-77 implementation of tableau matching (numeric)
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
* Then we use simulated annealing to improve the score. At each 
* iteration a random SSE is chosen to be remapped to a random
* other SSE (obeying constraints that are set) or (if LUNMAP is set)
* mapped to no SSE in the other structure.
*
* $Id: tsamtn.F 3242 2010-01-18 05:13:04Z alexs $
*=======================================================================

      subroutine tsamtn(omega1,ldt1,n1,omega2,ldt2,n2,
     $     ltype,lorder,dmat1,dmat2,score,ssemap,info)

      implicit none
*
* tsamtn - find match score of two tableaux using simulated annealing
*
*   The tableaux are specified here as omega matrices
*

*
*     .. Scalar Arguments ..
      integer ldt1,n1,ldt2,n2,info
      logical ltype,lorder
      double precision score
*     ..
*     .. Array Arguments ..
      double precision omega1(ldt1,*),omega2(ldt2,*)
      double precision dmat1(ldt1,*),dmat2(ldt2,*)
      integer ssemap(*)
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
* dmat1  (input) DOUBLE PRECISION array, dimension (n1,n1)
*        SSE distance matrix for one structure. symmetric
*
* dmat2  (input) DOUBLE PRECISION array, dimension (n2,n2)
*        SSE distance matrix for second structure. symmetric
*
* score  (output) DOUBLE PRECISION
*        score of matching omega1 with omega2 using the heuristic methods
*
* ssemap (output) INTEGER vector, dimension (n1)
*        solution SSE map vector of dimension n1. 
*        Each ssemap(i) is the SSE index
*        in omega2 that SSE i in omega1 is matched with.
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
      double precision bottom
      parameter (bottom = -99999.0d0)
*
*     limit on number of restarts
      integer maxstart
      parameter (maxstart = 100)
*
*     limit on number of iterations
      integer maxiter
      parameter (maxiter = 100)
*
*     largest tableau (omega matrix) dimension allowed
      integer maxdim
      parameter(maxdim = 110)
*     initial temperature
      double precision temp0
      parameter(temp0 = 10.0d0)
*     factor to multiply temperature by at each step
      double precision alpha
      parameter(alpha = 0.95d0)
*
*     ..
*     .. Local Scalars ..
      integer k,iter,ssei,oldj,newj
      integer startj,endj,restart
      logical lunmap
      double precision newscore,maxscore
      double precision temp
      double precision x
*     ..
*     .. Local Arrays ..
*     reverse ssemap: revmap(j) for j index in omega2 is the index i
*     in omega1 that matches that sse i.e. if ssemap(i) = j then
*     revmap(j) = ssemap(i) and vice versa, for quick lookup of what
*     is matched so we can easily check that one-to-one mapping maintained
*     revmap has dimension (n2)
      integer revmap(maxdim)
*     best ssemap found. this has dimenion (n1)
      integer bestmap(maxdim)
*     ..
*     .. Intrinsic Functions ..
      intrinsic exp,dble,min
#if defined(__INTEL_COMPILER) || defined(__PORTLAND_COMPILER)
      intrinsic RANDOM_NUMBER
#else
      intrinsic rand
#endif
*     ..
*     .. External Subroutines and Functions ..
      double precision  tmscorn
      external tmscorn
      integer  irrand
      external irrand
      external randtypeindn
      integer  randtypeindn
      external icopy
      external thinitn
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

      maxscore = bottom
      restart = 0
 10   if (restart .lt. maxstart) then
*        setup initial mapping
         call thinitn(omega1, ldt1, n1, omega2, ldt2, n2, ltype, lorder,
     $        lunmap, ssemap, revmap, info)
         if (info .ne. 0) then
*           constraints (order,type,unmap) don't allow any matching
            score = bottom
            return
         endif

         score = tmscorn(omega1, ldt1, n1, omega2, ldt2, n2, 
     $        dmat1, dmat2, ssemap)
         if (score .gt. maxscore) then
            maxscore = score
            call icopy(n1, ssemap, bestmap)
         endif
         iter = 0
         temp = temp0
 100     if (iter .lt. maxiter) then
*        generate neighbour state by picking random SSE in omega1
*        and moving its mapping to a random SSE in omega2, mainting
*        constraints
            ssei = irrand(n1)
*           TODO for now always obeying type constraint (LTYPE)
            if (lorder) then
               startj = ssemap(ssei)
               k = ssei
 110           if (startj .eq. 0 .and. k .gt. 0) then
                  startj = ssemap(k)
                  k = k - 1
                  goto 110
               endif
               if (startj .eq. 0) then
                  startj = n2+1
               endif
               if (ssei .eq. n1) then
                  endj = n2+1
               elseif (ssemap(ssei+1) .eq. 0) then
                  endj = 0
                  k = 1
 120              if (endj .eq. 0 .and. ssei + k .le. n1) then
                     endj = ssemap(ssei + k)
                     k = k + 1
                     goto 120
                  endif
               else
                  endj = ssemap(ssei+1)
               endif
            else
               startj = 1
               endj = n2
            endif
            newj = randtypeindn(omega2, ldt2, n2, startj, 
     $           omega1(ssei,ssei), revmap, endj)
C     write (*,*) (ssemap(k), k= 1, n1)
C     write (*,*) iter, ssei, startj, endj, newj
            oldj = ssemap(ssei)
            if (newj .ne. 0) then
               ssemap(ssei) = newj
               if (oldj .ne. 0) then
                  revmap(oldj) = 0
               endif
               revmap(newj) = ssei
            elseif (lunmap) then
*              also try removing the SSE from matching
               if (oldj .ne. 0) then
                  revmap(ssemap(ssei)) = 0
                  revmap(oldj) = 0
               endif
               ssemap(ssei) = 0
            endif

            newscore = tmscorn(omega1, ldt1, n1, omega2, ldt2, n2, 
     $           dmat1, dmat2, ssemap)
            if (newscore .gt. maxscore) then
               maxscore = newscore
               call icopy(n1, ssemap, bestmap)
            endif
C     write (*,*) iter, temp, score,newscore,
C     $        exp(dble(newscore - score) / temp)
#if defined(__INTEL_COMPILER) || defined(__PORTLAND_COMPILER)
            CALL RANDOM_NUMBER(x)
#else
            x = dble(rand(0))
#endif
            if (exp(dble(newscore - score) / temp).gt. x) then
               score = newscore
            else
               ssemap(ssei) = oldj
               if (oldj .ne. 0) then
                  revmap(oldj) = ssei
               endif
            endif
            
            temp = temp * alpha
            iter = iter + 1
*           end of while (iter .lt. maxiter) 
            goto 100
         endif
         restart = restart + 1
*        end of while (restart .lt. maxstart)
         goto 10
      endif

      score = maxscore
      call icopy(n1, bestmap, ssemap)
      info = 0
      return
      end
