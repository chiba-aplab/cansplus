c======================================================================|
      subroutine ncputs3(idf,idd,nd,data_dbl,ix,jx,kx)
c======================================================================|
c
c NAME  ncputs3
c
c PURPOSE
c    output sequential 3-D double-floating-point variable
c
c OUTPUTS
c    None 
c    
c INPUTS
c    idf: [integer] ID number for the output file
c    idd: [integer] ID number for the sequential variable
c    nd: [integer] 
c    data_dbl(ix,jx,kx): [double] value for variable
c    ix,jx,kx: [integer] dimension size
c    NOTE: idf are defined by subroutine 'ncopndm'
c          idd are defined by subroutine 'ncdefs3'
c 
c HISTORY
c    written 2002-3-1 T. Yokoyama
c     
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(mrank=4)
      dimension data_dbl(ix,jx,kx)
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

      mcount(1)=ix
      mcount(2)=jx
      mcount(3)=kx
      mcount(4)=1

      mstart(1)=1
      mstart(2)=1
      mstart(3)=1
      mstart(4)=nd

      mstatus=nf_put_vara_real(idf,idd,mstart,mcount,data)

      return
      end
