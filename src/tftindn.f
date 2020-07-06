*=======================================================================
* File:    tftindn.f
* Author:  Alex Stivala
* Created: November 2009
*
*  find the index of the first SSE of same type in omega matrix
*  that is not already mapped, or 0 if not found
*
* $Id: tftindn.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================

      integer function firsttypeindn(omega, ldim, n1, startind, ssetype,
     $     smap, endind)
      implicit none
*
* firsttypeindn - find the index of the first SSE of same type in omega matrix
*             that is not already mapped, or 0 if not found
*     
      integer       ldim
      integer       n1
      integer       startind
      integer       endind
      double precision   ssetype
      double precision  omega(ldim,*)
      integer       smap(*)

* Arguments
* =========
*
* omega   (input) DOUBLE PRECISION array, dimension (n1,n1)
*        Omega matrix for one structure. Symmetric.
*
* ldt   (input) INTEGER
*        Leading dimension of omega array. ldt >= n1
*
* n1     (input) INTEGER
*        Dimension of omega array.
* 
* startind (input) INTEGER
*        SSE index to start at in omega
*
* ssetype (input) DOUBLE PRECISION
*        SSE type as 1.0, 2.0, etc.
*
* smap    (input) INTEGER vector, dimension(n1)
*         each smap(i) is index in other omega matrix it is already mapped
*         to, or 0 for not mapped.
*
* endind (input) INTEGER
*         last SSE index to consider in omega 
*
      integer i
      logical found

      i = startind
      found = .false.
      firsttypeindn = 0
 1000 if (.not. found .and. (i .lt. endind)) then
         if (omega(i,i) .eq. ssetype .and. smap(i) .eq. 0) then
            found = .true.
            firsttypeindn = i
         endif
         i = i + 1
         goto 1000
      endif
      return
      end
