c======================================================================|
      subroutine ncputosc(idf,name,data_char)
c======================================================================|
c
c NAME  ncputosc
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
c    data_char: [double] value for variable
c    NOTE: idf are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c 
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) name
      character*(*) data_char
c----------------------------------------------------------------------|      

      mstatus = nf_redef(idf)
c     mstatus = nf_put_att_text(idf,0,name,len_trim(data_char),data_char)
      mstatus = nf_put_att_text(idf,0,name,len(data_char),data_char)
      mstatus = nf_enddef(idf)


      return
      end
