c======================================================================|
      subroutine dacopnr0s(idf,file,mtype,nx)
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
      read(idf) munity
      if (munity.ne.1) then
        write(6,*) 'Error:: endian different'
        stop
      endif

c   version
      read(idf) mversion
      if (mversion.ne.0) then
        write(6,*) 'Error:: version different. version=', mversion
        stop
      endif

c   for data type: integer=4, real=5, double=6
      read(idf) mtype

c   number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
c   note that the time-sequence should be included
      read(idf) mdim
      if (mdim.ne.1) then
        write(6,*) 'Error:: dimension different. mdim=', mdim
        stop
      endif

c   number of gridsize: 
c     unlimited=-1 (dependent on how many steps calculated)
      read(idf) nx

      return
      end
