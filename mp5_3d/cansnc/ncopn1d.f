c======================================================================|
      subroutine ncopn1d(idf,file,idn,idi,ix)
c======================================================================|
c
c NAME  ncopndm
c
c PURPOSE
c    open netCDF file and define dimension
c
c OUTPUTS
c    idf: [integer] ID number for the output file
c    idi: [integer] ID number for the coordinate dimension
c    idn: [integer] ID number for the temporal dimension
c
c INPUTS
c    file: [character] file name
c    ix: [integer] dimension size
c
c HISTORY
c    written 2002-3-1 T. Yokoyama
c
c----------------------------------------------------------------------|

      implicit double precision (a-h,o-z)
      character*(*) file
c----------------------------------------------------------------------|      

      mstatus = nf_create(file,0,idf)

      mstatus = nf_def_dim(idf,'ix',ix,idi)
      mstatus = nf_def_dim(idf,'nx', 0,idn)

      mstatus = nf_enddef(idf)


      return
      end
