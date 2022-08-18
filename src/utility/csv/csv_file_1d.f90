! cvs_file_1d.f90 --
!     Include file for csv_file.f90:
!     contains the body of the one-dimensional version of the
!     writing routines.
!
!     $Id: csv_file_1d.f90,v 1.2 2006/03/26 19:03:53 arjenmarkus Exp $
!
    integer, intent(in)                 :: lun
    logical, intent(in), optional       :: advance
    character(*), intent(in), optional   :: sep
    

    character(:), allocatable     :: sp
    logical                             :: adv
    integer                             :: i

    adv = .true.
    sp = ','

    if ( present(advance) ) adv = advance

    if (present(sep) ) sp = sep

    do i = 1,size(array)-1
        call csv_write( lun, array(i), sp, .false. )
    enddo

    call csv_write( lun, array(size(array)), sp, adv )
! 
! end of body
