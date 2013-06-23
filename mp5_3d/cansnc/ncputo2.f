c======================================================================|
      subroutine ncputo2(idf,idi,idj,name,data_dbl,ix,jx)
c======================================================================|
c
c NAME  ncputo2
c
c PURPOSE
c    output 2-D double-floating-point variable
c
c OUTPUTS
c    None 
c    
c INPUTS
c    idf: [integer] ID number for the output file
c    idi,idj: [integer] ID number for the coordinate dimension
c    name: [character] variable's name
c    data_dbl(ix): [double] value for variable
c    ix,jx: [integer] dimension size
c    NOTE: idf,idi,idj are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c     
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(mrank=2)
      character*(*) name
      dimension data_dbl(ix,jx)
      integer mdims(mrank)
      integer mstart(mrank),mcount(mrank)
      real data(ix,jx)
c----------------------------------------------------------------------|      

      do j=1,jx
      do i=1,ix
         data(i,j)=data_dbl(i,j)
      enddo
      enddo

      mdims(1)=idi
      mdims(2)=idj

      mstatus = nf_redef(idf)
      mstatus = nf_def_var(idf,name,5,mrank,mdims,idd)
      mstatus = nf_enddef(idf)

      mcount(1)=ix
      mcount(2)=jx

      mstart(1)=1
      mstart(2)=1

      mstatus=nf_put_vara_real(idf,idd,mstart,mcount,data)

      return
      end
