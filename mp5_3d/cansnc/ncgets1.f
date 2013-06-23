c======================================================================|
      subroutine ncgets1(idf,name,nd,data_dbl,ix,ndrd)
c======================================================================|
c
c NAME  ncgets1
c
c PURPOSE
c    input 1-D sequentially-output double-floating-type variable
c
c OUTPUTS
c    data_dbl(ix): [double] value for variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    idd: [integer] ID number for the defined variable
c    nd: [integer] 
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
      character*(*) name
c----------------------------------------------------------------------|      

      mstatus=nf_inq_varid(idf,name,idd)

      mstatus=nf_inq_unlimdim(idf,idnx)
      mstatus=nf_inq_dimlen(idf,idnx,nx)

      if (nd.eq.-1) then
        ndrd=nx
      else
        if (nd.le.nx) then
          ndrd=nd
        else
          ndrd=nx
        endif
      endif

      mcount(1)=ix
      mcount(2)=1

      mstart(1)=1
      mstart(2)=ndrd

      mstatus=nf_get_vara_real(idf,idd,mstart,mcount,data)

      do i=1,ix
         data_dbl(i)=data(i)
      enddo

      return
      end
