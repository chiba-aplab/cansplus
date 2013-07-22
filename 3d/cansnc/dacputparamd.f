c======================================================================|
      subroutine dacputparamd(mf_params,name,value)
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
      double precision value
c----------------------------------------------------------------------|

      write(mf_params,9020) 5,name,value
9020  format('#dac ',i4,a32,' ',e20.13)

      return
      end
