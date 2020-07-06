*=======================================================================
* File:    tftind.f
* Author:  Alex Stivala
* Created: November 2009
*
*  find the index of the first SSE of same type in tableau
*  that is not already mapped, or 0 if not found
*
* $Id: tftind.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================

      integer function firsttypeind(tab, ldim, n1, startind, ssetype,
     $     smap, endind)
      implicit none
*
* firsttypeind - find the index of the first SSE of same type in tableau
*             that is not already mapped, or 0 if not found
*     
      integer       ldim
      integer       n1
      integer       startind
      integer       endind
      character*2   ssetype
      character*2   tab(ldim,*)
      integer       smap(*)

* Arguments
* =========
*
* tab   (input) CHARACHTER*2 array, dimension (n1,n1)
*        Tableau for one structure. Symmetric.
*
* ldt   (input) INTEGER
*        Leading dimension of tab array. ldt >= n1
*
* n1     (input) INTEGER
*        Dimension of tab array.
* 
* startind (input) INTEGER
*        SSE index to start at in tab
*
* ssetype (input) CHARACTER*2
*        SSE type as two charcter string 'xa' etc.
*
* smap    (input) INTEGER vector, dimension(n1)
*         each smap(i) is index in other tableau it is already mapped
*         to, or 0 for not mapped.
*
* endind (input) INTEGER
*         last SSE index to consider in tab 
*
      integer i
      logical found

      i = startind
      found = .false.
      firsttypeind = 0
 1000 if (.not. found .and. (i .lt. endind)) then
         if (tab(i,i) .eq. ssetype .and. smap(i) .eq. 0) then
            found = .true.
            firsttypeind = i
         endif
         i = i + 1
         goto 1000
      endif
      return
      end
