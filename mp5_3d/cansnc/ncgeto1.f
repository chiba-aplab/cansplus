c======================================================================|
      subroutine ncgeto1(idf,name,data_dbl,ix)
c======================================================================|
c
c NAME  ncgeto1
c
c PURPOSE
c    input 1-D double-floating-type variable
c
c OUTPUTS
c    data_dbl(ix): [double] value for variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    ix: [integer] dimension size
c    NOTE: idf are defined by subroutine 'ncopndm'
c 
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(mrank=1)
      character*(*) name
      dimension data_dbl(ix)
      integer mstart(mrank),mcount(mrank)
      real data(ix)
c----------------------------------------------------------------------|      

      mstatus=nf_inq_varid(idf,name,idd)

      mcount(1)=ix
      mstart(1)=1

      mstatus=nf_get_vara_real(idf,idd,mstart,mcount,data)

      do i=1,ix
         data_dbl(i)=data(i)
      enddo

      return
      end
