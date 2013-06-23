c======================================================================|
      subroutine ncdefssi(idf,name,idd,idn)
c======================================================================|
c
c NAME  ncdefssi
c
c PURPOSE
c    define scalar sequentially-output variable
c
c OUTPUTS
c    idd: [integer] ID number for the defined variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    idn: [integer] ID number for the temporal dimension
c    NOTE: idf,idn are defined by subroutine 'ncopndm'
c    
c HISTORY
c    written 2002-3-1 T. Yokoyama
c    
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      character*(*) name
c----------------------------------------------------------------------|

      mstatus = nf_redef(idf)
      mstatus = nf_def_var(idf,name,4,1,idn,idd)
      mstatus = nf_enddef(idf)


      return
      end
