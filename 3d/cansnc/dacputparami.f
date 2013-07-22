c======================================================================|
      subroutine dacputparami(mf_params,name,mvalue)
c======================================================================|
c
c NAME  ncdefss
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

      write(mf_params,9010) 4,name,mvalue
9010  format('#dac ',i4,a32,' ',i20)

      return
      end
