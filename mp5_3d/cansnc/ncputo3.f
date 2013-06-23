c======================================================================|
      subroutine ncputo3(idf,idi,idj,idk,name,data_dbl,ix,jx,kx)
c======================================================================|
c
c NAME  ncputo3
c
c PURPOSE
c    output 3-D double-floating-point variable
c
c OUTPUTS
c    None 
c    
c INPUTS
c    idf: [integer] ID number for the output file
c    idi,idj,idk: [integer] ID number for the coordinate dimension
c    name: [character] variable's name
c    data_dbl(ix,jx,kx): [double] value for variable
c    ix,jx,kx: [integer] dimension size
c    NOTE: idf,idi,idj,idk are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c     
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(mrank=3)
      character*(*) name
      dimension data_dbl(ix,jx,kx)
      integer mdims(mrank)
      integer mstart(mrank),mcount(mrank)
      real data(ix,jx,kx)
c----------------------------------------------------------------------|      

      do k=1,kx
      do j=1,jx
      do i=1,ix
         data(i,j,k)=data_dbl(i,j,k)
      enddo
      enddo
      enddo

      mdims(1)=idi
      mdims(2)=idj
      mdims(3)=idk

      mstatus = nf_redef(idf)
      mstatus = nf_def_var(idf,name,5,mrank,mdims,idd)
      mstatus = nf_enddef(idf)

      mcount(1)=ix
      mcount(2)=jx
      mcount(3)=kx

      mstart(1)=1
      mstart(2)=1
      mstart(3)=1

      mstatus=nf_put_vara_real(idf,idd,mstart,mcount,data)

      return
      end
