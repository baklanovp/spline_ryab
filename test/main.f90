program MAIN
    use datas
    use RYABMOD
    implicit none
    real(8) FUNCINTERP
    real(8) s,per(NDIM)
    real :: t1, t2

    integer K, c

    call load

    open(11,file='result.dat')

    ! do K=1,n_tp
    !   write(11,'(10es23.15)') TpTab(K),arr_dump(K,5,3)
    ! end do

    ! close(11)
    ! stop


    c = 0
    s=0.d0
    CALL CPU_TIME(t1)

    do while(s.le.1.d0)
      c = c+1
      per(1)=TpTab(2)+s*(TpTab(n_tp-1)-TpTab(2))
      per(2)=RhoTab(5) !RhoTab(2)+s*(RhoTab(n_rho-1)-RhoTab(2))
      per(3)=lnTimeTab(3) !lnTimeTab(2)+s*(lnTimeTab(n_times-1)-lnTimeTab(2))

      write(11,'(10es23.15)') per(1),FUNCINTERP(per)
      ! print*,s
      !  read*
      s=s+1.d-4
    end do

    close(11)

    ! Code segment to be timed
    CALL CPU_TIME(t2)

    PRINT *, 'Time taken for ', c, ' calls of FUNCINTERP: ', t2 - t1, ' seconds.'

    stop
end program
