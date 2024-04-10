module mod_TimeSchemes

    use mod_constantes
    use mod_flux
    

    implicit none

    contains


    function RK2SSP(tn, Un, dt, imax, bi, dx)result(Unp1)

        integer, intent(in)                         :: imax
        real(pr), intent(in)                        :: dt, tn, dx
        real(pr), dimension(:), intent(in)          :: bi
        real(pr), dimension(0:imax+1,2), intent(in) :: Un

        real(pr), dimension(1:imax,2)               :: Unp1


        ! -- Local
        real(pr), dimension(0:imax+1, 2) :: k1
        real(pr), dimension(1:imax, 2)   :: k2

        k1(1:imax,:) = Application(3, imax, tn, Un, bi, dx)
        k1(0,:) = k1(1,:)
        k1(imax+1,:) = k1(imax,:)
        ! k1(0,:) = 0.
        ! k1(imax+1,:) = 0.
        k2 = Application(3, imax, dt + tn, Un + dt*k1, bi, dx)

        Unp1(1:imax,:) = Un(1:imax,:) + dt/2*(k1(1:imax,:)+k2(1:imax,:))

    end function



end module