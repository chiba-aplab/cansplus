c======================================================================|
      subroutine ncgetssi(idf,name,nd,mdata_int,ndrd)
c======================================================================|
c
c NAME  ncgetssi
c
c PURPOSE
c    input scalar sequentially-output integer variable
c
c OUTPUTS
c    data_dbl: [double] value for variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    idd: [integer] ID number for the defined variable
c    nd: [integer] 
c    NOTE: idf are defined by subroutine 'ncopndm'
c          idd are defined by subroutine 'ncdefs1'
c          
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
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

      mstatus=nf_get_var1_int(idf,idd,ndrd,mdata_int)

      return
      end
