c======================================================================|
      subroutine ncputos(idf,name,data_dbl)
c======================================================================|
c
c NAME  ncputos
c
c PURPOSE
c    output scalar double-floating-point variable
c
c OUTPUTS
c    None 
c    
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    data_dbl: [double] value for variable
c    NOTE: idf are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c 
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) name
c----------------------------------------------------------------------|      

      mstatus = nf_redef(idf)
      mstatus = nf_def_var(idf,name,5,0,0,idd)
      mstatus = nf_enddef(idf)

      mstatus=nf_put_var_real(idf,idd,real(data_dbl))

      return
      end
