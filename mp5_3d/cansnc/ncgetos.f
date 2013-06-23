c======================================================================|
      subroutine ncgetos(idf,name,data_dbl)
c======================================================================|
c
c NAME  ncgetos
c
c PURPOSE
c    input scalar double-floating-type variable
c
c OUTPUTS
c    data_dbl: [double] value for variable
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
      real data
c----------------------------------------------------------------------|      

      mstatus=nf_inq_varid(idf,name,idd)

      mstatus=nf_get_var_real(idf,idd,data)

      data_dbl=data

      return
      end
