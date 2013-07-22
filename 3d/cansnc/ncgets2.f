c======================================================================|
      subroutine ncgets2(idf,name,nd,data_dbl,ix,jx,ndrd)
c======================================================================|
c
c NAME  ncgets2
c
c PURPOSE
c    input 2-D sequentially-output double-floating-type variable
c
c OUTPUTS
c    data_dbl(ix,jx): [double] value for variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    idd: [integer] ID number for the defined variable
c    nd: [integer] 
c    ix,jx: [integer] dimension size
c    NOTE: idf are defined by subroutine 'ncopndm'
c          idd are defined by subroutine 'ncdefs2'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c     
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(mrank=3)
      dimension data_dbl(ix,jx)
      integer mstart(mrank),mcount(mrank)
      real data(ix,jx)
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
      mcount(2)=jx
      mcount(3)=1

      mstart(1)=1
      mstart(2)=1
      mstart(3)=ndrd

      mstatus=nf_get_vara_real(idf,idd,mstart,mcount,data)

      do j=1,jx
      do i=1,ix
         data_dbl(i,j)=data(i,j)
      enddo
      enddo

      return
      end
