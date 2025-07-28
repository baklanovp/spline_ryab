module datas
 real(8), dimension(:), allocatable :: TpTab, RhoTab, lnTimeTab
 real(8), dimension(:,:,:), allocatable :: arr_dump
 integer :: n_tp, n_rho, n_times
end module

!========================================================
!********************************************************

    subroutine  load
	    use datas
        implicit none
        character(30) fname
        integer :: ierr, ui
        integer :: it

        fname = 'data/neM20Ni01Z002.dump'

        write(*,"(A,A)") 'Loading from file: ', trim(fname);

        open(10, file=trim(fname), status='unknown', form='formatted');
        read(10,*) n_tp, n_rho, n_times
        write(*,*) ' n_tp, n_rho, n_times: ', n_tp, n_rho, n_times

        allocate(TpTab(n_tp),RhoTab(n_rho),lnTimeTab(n_times))
        allocate(arr_dump(n_tp,n_rho,n_times))

        read(10,*) TpTab
        read(10,*) RhoTab
        read(10,*) lnTimeTab
        read(10,*) arr_dump

        close(10)

        write(*,"(A,A/)") ' Show data ';

        
        write(*,*) 'TpTab: ', TpTab
!		read*
        write(*,*) 'RhoTab: ', RhoTab
!		read*
        write(*,*) 'lnTimeTab: ', lnTimeTab
!		read*
        write(*,*) 'arr_dump: ', arr_dump
        ! read*

    end subroutine  load
    