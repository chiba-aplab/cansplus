module dac_header
  implicit none
  private :: dacdefprep,dacopnrprep
  
  public :: dacdef0s,dacdef1s,dacdef1d,dacdef2s,dacdef2d,dacdef3s,dacdef3d &
       ,dacopnr0s,dacopnr1s,dacopnr1d,dacopnr2s,dacopnr3s &
       ,dacputparamc,dacputparami,dacputparamd

contains


  subroutine dacdefprep(idf,file,mtype)
    implicit none
    integer, intent(in) :: idf,mtype
    character(*),intent(in) :: file

    open(idf,file=file,form='unformatted')
! for byteorder
    write(idf) 1
! version
    write(idf) 0
! for data type: integer=4, real=5, double=6
    write(idf) mtype

  end subroutine dacdefprep


  subroutine dacopnrprep(idf,file,mtype)
    implicit none
    integer,intent(in) :: idf
    integer,intent(inout) :: mtype
    character(*),intent(in) :: file
    integer :: munity,mversion
    
    open(idf,file=file,form='unformatted')
! for byteorder                                                               
    read(idf) munity
    if (munity /= 1) then
       write(6,*) 'Error:: endian different'
       stop
    endif
! version                                                                     
    read(idf) mversion
    if (mversion /=0) then
       write(6,*) 'Error:: version different. version=', mversion
       stop
    endif
! for data type: integer=4, real=5, double=6                                  
    read(idf) mtype
    
  end subroutine dacopnrprep


  subroutine dacdef0s(idf,file,mtype)
    implicit none
    integer,intent(in) :: idf,mtype
    character(*),intent(in) :: file

    call dacdefprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
! note that the time-sequence should be included
    write(idf) 1
! number of gridsize: 
! unlimited=-1 (dependent on how many steps calculated)
    write(idf) -1

    return
  end subroutine dacdef0s


  subroutine dacdef1d(idf,file,mtype,ix)
    implicit none
    integer,intent(in) :: idf,mtype,ix
    character(*),intent(in) :: file
    
    call dacdefprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
! note that the time-sequence should be included
    write(idf) 1
! number of gridsize: 
! unlimited=-1 (dependent on how many steps calculated)
    write(idf) ix

    return
  end subroutine dacdef1d


  subroutine dacdef1s(idf,file,mtype,ix)
    implicit none
    integer,intent(in) :: idf,mtype,ix
    character(*),intent(in) :: file
    
    call dacdefprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
! note that the time-sequence should be included
    write(idf) 2
! number of gridsize: 
! unlimited=-1 (dependent on how many steps calculated)
    write(idf) ix,-1

    return
  end subroutine dacdef1s


  subroutine dacdef2d(idf,file,mtype,ix,jx)
    implicit none
    integer,intent(in) :: idf,mtype,ix,jx
    character(*),intent(in) :: file
    
    call dacdefprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
! note that the time-sequence should be included
    write(idf) 2
! number of gridsize: 
! unlimited=-1 (dependent on how many steps calculated)
    write(idf) ix,jx

    return
  end subroutine dacdef2d


  subroutine dacdef2s(idf,file,mtype,ix,jx)
    implicit none
    integer,intent(in) :: idf,mtype,ix,jx
    character(*),intent(in) :: file
    
    call dacdefprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
! note that the time-sequence should be included
    write(idf) 3
! number of gridsize: 
! unlimited=-1 (dependent on how many steps calculated)
    write(idf) ix,jx,-1

    return
  end subroutine dacdef2s


  subroutine dacdef3d(idf,file,mtype,ix,jx,kx)
    implicit none
    integer,intent(in) :: idf,mtype,ix,jx,kx
    character(*),intent(in) :: file

    call dacdefprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
! note that the time-sequence should be included
    write(idf) 3
! number of gridsize: 
! unlimited=-1 (dependent on how many steps calculated)
    write(idf) ix,jx,kx

    return
  end subroutine dacdef3d


  subroutine dacdef3s(idf,file,mtype,ix,jx,kx)
    implicit none
    integer,intent(in) :: idf,mtype,ix,jx,kx
    character(*),intent(in) :: file
    
    call dacdefprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...
! note that the time-sequence should be included
    write(idf) 4
! number of gridsize: 
! unlimited=-1 (dependent on how many steps calculated)
    write(idf) ix,jx,kx,-1

    return
  end subroutine dacdef3s


  subroutine dacopnr0s(idf,file,mtype,nx)
    implicit none
    integer,intent(in) :: idf
    integer,intent(inout) :: mtype,nx
    character(*),intent(in) :: file
    integer :: mdim

    call dacopnrprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...            
! note that the time-sequence should be included                              
    read(idf) mdim
    if (mdim /= 1) then
       write(6,*) 'Error:: dimension different. mdim=', mdim
       stop
    endif
! number of gridsize:                                                         
! unlimited=-1 (dependent on how many steps calculated)                     
    read(idf) nx

    return
  end subroutine dacopnr0s


  subroutine dacopnr1s(idf,file,mtype,ix,nx)
    implicit none
    integer,intent(in) :: idf
    integer,intent(inout) :: mtype,ix,nx
    character(*),intent(in) :: file
    integer :: mdim

    call dacopnrprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...            
! note that the time-sequence should be included                              
    read(idf) mdim
    if (mdim /= 2) then
       write(6,*) 'Error:: dimension different. mdim=', mdim
       stop
    endif
! number of gridsize:                                                         
! unlimited=-1 (dependent on how many steps calculated)                     
    read(idf) ix,nx

    return
  end subroutine dacopnr1s


  subroutine dacopnr1d(idf,file,mtype,ix)
    implicit none
    integer,intent(in) :: idf
    integer,intent(inout) :: mtype,ix
    character(*),intent(in) :: file
    integer :: mdim

    call dacopnrprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...            
! note that the time-sequence should be included                              
    read(idf) mdim
    if (mdim /= 1) then
       write(6,*) 'Error:: dimension different. mdim=', mdim
       stop
    endif
! number of gridsize:                                                         
! unlimited=-1 (dependent on how many steps calculated)                     
    read(idf) ix

    return
  end subroutine dacopnr1d


  subroutine dacopnr2s(idf,file,mtype,ix,jx,nx)
    implicit none
    integer,intent(in) :: idf
    integer,intent(inout) :: mtype,ix,jx,nx
    character(*),intent(in) :: file
    integer :: mdim

    call dacopnrprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...            
! note that the time-sequence should be included                              
    read(idf) mdim
    if (mdim /= 3) then
       write(6,*) 'Error:: dimension different. mdim=', mdim
       stop
    endif
! number of gridsize:                                                         
! unlimited=-1 (dependent on how many steps calculated)                     
    read(idf) ix,jx,nx

    return
  end subroutine dacopnr2s


  subroutine dacopnr3s(idf,file,mtype,ix,jx,kx,nx)
    implicit none
    integer,intent(in) :: idf
    integer,intent(inout) :: mtype,ix,jx,kx,nx
    character(*),intent(in) :: file
    integer :: mdim

    call dacopnrprep(idf,file,mtype)
! number of dimension: scalar=0, 1-dimension=1, 2-dimension=2, ...            
! note that the time-sequence should be included                              
    read(idf) mdim
    if (mdim /= 4) then
       write(6,*) 'Error:: dimension different. mdim=', mdim
       stop
    endif
! number of gridsize:                                                         
! unlimited=-1 (dependent on how many steps calculated)                     
    read(idf) ix,jx,kx,nx

    return
  end subroutine dacopnr3s


  subroutine dacputparamc(idf,name,put_value)
    implicit none
    integer,intent(in) :: idf
    character(*),intent(in) :: name,put_value

    write(idf,9030) 2,name,put_value
9030 format('#dac ',i4,a32,' ',a20)

    return
  end subroutine dacputparamc


  subroutine dacputparamd(idf,name,put_value)
    implicit none
    integer,intent(in) :: idf
    character(*),intent(in) :: name
    real(8),intent(in) :: put_value

    write(idf,9020) 5,name,put_value
9020 format('#dac ',i4,a32,' ',e20.13)

    return
  end subroutine dacputparamd


  subroutine dacputparami(idf,name,put_value)
    implicit none
    integer,intent(in) :: idf,put_value
    character(*),intent(in) :: name

    write(idf,9010) 4,name,put_value
9010 format('#dac ',i4,a32,' ',i20)

    return
  end subroutine dacputparami


end module dac_header
