c======================================================================|
      subroutine ncgeto2(idf,name,data_dbl,ix,jx)
c======================================================================|
c
c NAME  ncgeto2
c
c PURPOSE
c    input 2-D double-floating-type variable
c
c OUTPUTS
c    data_dbl(ix,jx): [double] value for variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    ix,jx: [integer] dimension size
c    NOTE: idf are defined by subroutine 'ncopndm'
c 
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(mrank=2)
      character*(*) name
      dimension data_dbl(ix,jx)
      integer mstart(mrank),mcount(mrank)
      real data(ix,jx)
c----------------------------------------------------------------------|      

      mstatus=nf_inq_varid(idf,name,idd)

      mcount(1)=ix
      mcount(2)=jx

      mstart(1)=1
      mstart(2)=1

      mstatus=nf_get_vara_real(idf,idd,mstart,mcount,data)

      do j=1,jx
      do i=1,ix
         data_dbl(i,j)=data(i,j)
      enddo
      enddo

      return
      end
