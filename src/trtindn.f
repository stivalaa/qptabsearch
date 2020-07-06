*=======================================================================
* File:    trtindn.f
* Author:  Alex Stivala
* Created: November 2009
*
*  find the index of a random SSE of same type in tableau
*  that is not already mapped, or 0 if not found
*
* $Id: trtindn.F 2969 2009-11-22 01:29:39Z astivala $
*=======================================================================

      integer function randtypeindn(omega, ldim, n1, startind, ssetype,
     $     smap, endind)
      implicit none
*
* randtypeindn - find the index of a random SSE of same type in tableau
*             that is not already mapped, or 0 if not found
*     
      integer       ldim
      integer       n1
      integer       startind
      integer       endind
      character*2   ssetype
      character*2   omega(ldim,*)
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
*        Dimension of tab array.
* 
* startind (input) INTEGER
*        SSE index to start at in omega
*
* ssetype (input) DOUBLE PRECISION
*        SSE type as 1.0, 2.0, etc.
*
* smap    (input) INTEGER vector, dimension(n1)
*         each smap(i) is index in other tableau it is already mapped
*         to, or 0 for not mapped.
*
* endind (input) INTEGER
*         last SSE index to consider in omega 
*
      integer  irrand
      external irrand
      
      integer maxdim
      parameter(maxdim=110)

      integer indlist(maxdim)
      integer i,indi

      i = startind
      indi = 1
      randtypeindn = 0
 1000 if (i .lt. endind) then
         if (omega(i,i) .eq. ssetype .and. smap(i) .eq. 0) then
            indlist(indi) = i
            indi = indi + 1
         endif
         i = i + 1
         goto 1000
      endif

      if (indi .eq. 2) then
         randtypeindn = indlist(1)
      elseif (indi .gt. 2) then
         randtypeindn = indlist(irrand(indi-1))
      endif
      end
