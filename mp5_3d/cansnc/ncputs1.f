c======================================================================|
      subroutine ncputs1(idf,idd,nd,data_dbl,ix)
c======================================================================|
c
c NAME  ncputs1
c
c PURPOSE
c    output sequential 1-D double-floating-point variable
c
c OUTPUTS
c    None 
c    
c INPUTS
c    idf: [integer] ID number for the output file
c    idd: [integer] ID number for the sequential variable
c    nd: [integer] 
c    data_dbl(ix): [double] value for variable
c    ix: [integer] dimension size
c    NOTE: idf are defined by subroutine 'ncopndm'
c          idd are defined by subroutine 'ncdefs1'
c 
c HISTORY
c    written 2002-3-1 T. Yokoyama
c     
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(mrank=2)
      dimension data_dbl(ix)
      integer mstart(mrank),mcount(mrank)
      real data(ix)
c----------------------------------------------------------------------|      

      do i=1,ix
         data(i)=data_dbl(i)
      enddo

      mcount(1)=ix
      mcount(2)=1

      mstart(1)=1
      mstart(2)=nd

      mstatus=nf_put_vara_real(idf,idd,mstart,mcount,data)

      return
      end
