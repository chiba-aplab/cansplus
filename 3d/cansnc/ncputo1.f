c======================================================================|
      subroutine ncputo1(idf,idi,name,data_dbl,ix)
c======================================================================|
c
c NAME  ncputo1
c
c PURPOSE
c    output 1-D double-floating-point variable
c
c OUTPUTS
c    None 
c    
c INPUTS
c    idf: [integer] ID number for the output file
c    idi: [integer] ID number for the coordinate dimension
c    name: [character] variable's name
c    data_dbl(ix): [double] value for variable
c    ix: [integer] dimension size
c    NOTE: idf,idi are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c     
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) name
      dimension data_dbl(ix)
      real data(ix)
c----------------------------------------------------------------------|      

      do i=1,ix
        data(i)=data_dbl(i)
      enddo

      mstatus = nf_redef(idf)
      mstatus = nf_def_var(idf,name,5,1,idi,idd)
      mstatus = nf_enddef(idf)

      mstatus=nf_put_var_real(idf,idd,data)

      return
      end
