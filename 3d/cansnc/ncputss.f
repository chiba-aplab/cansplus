c======================================================================|
      subroutine ncputss(idf,idd,nd,data_dbl)
c======================================================================|
c
c NAME  ncputss
c
c PURPOSE
c    output sequential scalar double-floating-point variable
c
c OUTPUTS
c    None 
c    
c INPUTS
c    idf: [integer] ID number for the output file
c    idd: [integer] ID number for the sequential variable
c    nd: [integer] 
c    data_dbl: [double] value for variable
c    NOTE: idf are defined by subroutine 'ncopndm'
c          idd are defined by subroutine 'ncdefs1'
c          
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      real data
c----------------------------------------------------------------------|      

      data=data_dbl
      mstatus=nf_put_var1_real(idf,idd,nd,data)

      return
      end
