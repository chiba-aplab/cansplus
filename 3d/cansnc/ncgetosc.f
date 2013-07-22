c======================================================================|
      subroutine ncgetosc(idf,name,data_char)
c======================================================================|
c
c NAME  ncgetosc
c
c PURPOSE
c    input scalar character-type variable
c
c OUTPUTS
c    data_char: [double] value for variable
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
      character*(*) data_char
c----------------------------------------------------------------------|      

      mstatus = nf_get_att_text(idf,0,name,data_char)


      return
      end
