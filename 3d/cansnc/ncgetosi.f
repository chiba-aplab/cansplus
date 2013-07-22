c======================================================================|
      subroutine ncgetosi(idf,name,mdata_int)
c======================================================================|
c
c NAME  ncgetosi
c
c PURPOSE
c    input scalar integer-type variable
c
c OUTPUTS
c    mdata_int: [integer] value for variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    NOTE: idf are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) name
c----------------------------------------------------------------------|      

      mstatus=nf_inq_varid(idf,name,idd)

      mstatus=nf_get_var_int(idf,idd,mdata_int)

      return
      end
