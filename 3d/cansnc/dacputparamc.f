c======================================================================|
      subroutine dacputparamc(idf,name,value)
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
      character*(*) value
c----------------------------------------------------------------------|

      write(idf,9030) 2,name,value
9030  format('#dac ',i4,a32,' ',a20)

      return
      end
