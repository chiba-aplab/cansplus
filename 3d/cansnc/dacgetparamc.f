c======================================================================|
      subroutine dacgetparamc(mf_params,mline,value)
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
      character*32 name
      character*60 cdummy
      character*20 value
c----------------------------------------------------------------------|

      rewind(mf_params)
      do m=1,mline-1
        read(mf_params,9000) cdummy
      enddo
      read(mf_params,9030) mtype,name,mvalue
9000  format(a60)
9030  format('#dac ',i4,a32,' ',a20)

      return
      end
