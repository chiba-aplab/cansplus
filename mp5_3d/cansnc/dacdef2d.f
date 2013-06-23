c======================================================================|
      subroutine dacdef2d(idf,file,mtype,ix,jx)
c======================================================================|
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
      write(idf) 2

c   number of gridsize: 
c     unlimited=-1 (dependent on how many steps calculated)
      write(idf) ix,jx

      return
      end
