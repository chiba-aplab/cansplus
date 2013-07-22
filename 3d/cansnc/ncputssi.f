c======================================================================|
      subroutine ncputssi(idf,idd,nd,mdata_int)
c======================================================================|
c
c NAME  ncputssi
c
c PURPOSE
c    output sequential scalar integer variable
c
c OUTPUTS
c    None 
c    
c INPUTS
c    idf: [integer] ID number for the output file
c    idd: [integer] ID number for the sequential variable
c    nd: [integer] 
c    mdata_int: [integer] value for variable
c    NOTE: idf are defined by subroutine 'ncopndm'
c          idd are defined by subroutine 'ncdefs1'
c          
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------|      

      mstatus=nf_put_var1_int(idf,idd,nd,mdata_int)

      return
      end
