program test_loaddata
    use,intrinsic :: iso_fortran_env, only : int16, int32, real32, real64

    implicit none
   
    character(len=*), parameter  :: mdl_name = 'test_loaddata'
    integer, parameter :: ip = int32
    integer, parameter :: sp = real32
    integer, parameter :: dp = real64
    

    call test_load


    contains

    subroutine  test_load
        character(len=*), parameter ::  subrtn_name = 'test_load', &
                            fullPathSubrtn = mdl_name//'.'//subrtn_name

        character(len=:), allocatable :: fname;
        real(dp), dimension(:), allocatable :: TpTab, RhoTab, lnTimeTab
        real(dp), dimension(:,:,:), allocatable :: arr_dump
        integer :: n_tp, n_rho, n_times
        integer :: ierr, ui
        integer :: it

        fname = 'data/neM20Ni01Z002.dump'

        write(*,"(A,A)") 'Loading from file: ', trim(fname);

        open(newunit=ui, iostat=ierr, file=fname, status='unknown', form='formatted');
        read(ui,*) n_tp, n_rho, n_times
        write(*,*) ' n_tp, n_rho, n_times: ', n_tp, n_rho, n_times

        call dalloc1d('TpTab',  n_tp, TpTab, initial=0._dp, path=fullPathSubrtn)
        call dalloc1d('RhoTab',  n_rho, RhoTab, initial=0._dp, path=fullPathSubrtn)
        call dalloc1d('lnTimeTab',  n_times, lnTimeTab, initial=0._dp, path=fullPathSubrtn)

        allocate(arr_dump(n_tp,n_rho,n_times), STAT=ierr);

        read(ui,*) TpTab;
        read(ui,*) RhoTab;
        read(ui,*) lnTimeTab
        read(ui,*) arr_dump

        close(ui)

        write(*,"(A,A/)") ' Show data ';

        
        write(*,*) 'TpTab: ', TpTab
        write(*,*) 'RhoTab: ', RhoTab
        write(*,*) 'lnTimeTab: ', lnTimeTab
        write(*,*) 'arr_dump: ', arr_dump

        ! do it = 1, n_times
        !     write(*,*) it, exp(lnTimeTab(it))

        !     ! write(*,*) arr_dump
        ! enddo

    end subroutine  test_load
    

    subroutine dalloc1d(name, n, a, initial, path)
      character(len=*), intent(in) :: name
      integer, intent(in) :: n
      real(dp), intent(in), optional :: initial
      character(len=*), optional ::  path;
      real(dp), dimension(:), allocatable,  intent(out) :: a

      character(len=*), parameter ::  subrtn_name = 'dalloc1d', &
                     fullPathSubrtn = mdl_name//'.'//subrtn_name
    !   character(len=999) ::  str;
      character(len=99) ::  l_path;
      integer :: ierr;

      l_path = fullPathSubrtn
      if( present(path)) l_path = path;

      allocate (a(n), STAT=ierr)
      if (ierr /= 0) then
            write(*, '(4a, i4)') fullPathSubrtn, ' Not enough memory for ', name, ' N=',n;
            error stop 666;
      endif

      a(:) = 0._dp
      if ( present(initial)) then
            a(:) = initial
      endif
   end subroutine dalloc1d

end program test_loaddata