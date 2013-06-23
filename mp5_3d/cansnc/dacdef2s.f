c======================================================================|
      subroutine dacdef2s(idf,file,mtype,ix,jx)
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
      character*(*) file
c----------------------------------------------------------------------|

      open(idf,file=file,form='unformatted')

c   for byteorder
      write(idf) 1

c   version
      write(idf) 0

c   for data type: integer=4, real=5, double=6
      write(idf) mtype

c   number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
c   note that the time-sequence should be included
      write(idf) 3

c   number of gridsize: 
c     unlimited=-1 (dependent on how many steps calculated)
      write(idf) ix,jx,-1

      return
      end
