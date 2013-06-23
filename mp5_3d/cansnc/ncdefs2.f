c======================================================================|
      subroutine ncdefs2(idf,name,idd,idn,idi,idj)
c======================================================================|
c
c NAME  ncdefs2
c
c PURPOSE
c    define 2-D sequentially-output variable
c
c OUTPUTS
c    idd: [integer] ID number for the defined variable
c
c INPUTS
c    idf: [integer] ID number for the output file
c    name: [character] variable's name
c    idn: [integer] ID number for the temporal dimension
c    idi,idj: [integer] ID number for the coordinate dimension
c    NOTE: idf,idn,idi,idj are defined by subroutine 'ncopndm'
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|
      implicit double precision (a-h,o-z)
      parameter(mrank=3)
      integer mdims(mrank)
      character*(*) name
c----------------------------------------------------------------------|

      mdims(1)=idi
      mdims(2)=idj
      mdims(3)=idn

      mstatus = nf_redef(idf)
      mstatus = nf_def_var(idf,name,5,mrank,mdims,idd)
      mstatus = nf_enddef(idf)


      return
      end
