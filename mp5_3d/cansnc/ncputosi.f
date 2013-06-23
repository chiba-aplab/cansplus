c======================================================================|
      subroutine ncputosi(idf,name,mdata_int)
c======================================================================|
c
c NAME  ncputosci
c
c PURPOSE
c    output scalar integer variable
c
c OUTPUTS
c    None 
c    
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    mdata_int: [double] value for variable
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
      mstatus = nf_def_var(idf,name,4,0,0,idd)
      mstatus = nf_enddef(idf)

      mstatus=nf_put_var_int(idf,idd,mdata_int)

      return
      end
