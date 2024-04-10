module mod_flux

    use mod_constantes

    implicit none


contains


    function FluxNum(imax, Ug, Ud, bi)result(F)

        integer, intent(in)                          :: imax
        real(pr), dimension(1:imax+1, 2), intent(in) :: Ud
        real(pr), dimension(0:imax, 2), intent(in)   :: Ug
        real(pr), dimension(imax+1), intent(in)      :: bi

        real(pr), dimension(1:imax+1, 2)             :: F

        F(1:imax+1,1) = 0.5*(Ud(1:imax+1,2) + Ug(0:imax,2)) - 0.5*bi(1:imax+1)*(Ud(1:imax+1,1) - Ug(0:imax,1))
        F(1:imax+1,2) = (Ug(0:imax,2)*Ug(0:imax,2)/Ug(0:imax,1) + &
        g*Ug(0:imax,1)*Ug(0:imax,1)/2. + Ud(1:imax+1,2)*Ud(1:imax+1,2)/Ud(1:imax+1,1)&
        + g*Ud(1:imax+1,1)*Ud(1:imax+1,1)/2.)*0.5 - 0.5*bi(1:imax+1)*(Ud(1:imax+1,2) - Ug(0:imax,2))

    end function


    ! -- Clacul du coefficient bi
    function f_bi(hg, hd, ug, ud)result(bi)
        real(pr), intent(in)                :: hg, hd
        real(pr), intent(in)                :: ud, ug
        real(pr)                            :: bi

        bi = max(abs(ug + sqrt(g*hg)), abs(ud + sqrt(g*hd)), abs(ug - sqrt(g*hg)), abs(ud - sqrt(g*hd)))

    end function f_bi


    subroutine RightLeftValues(Un, imax, Ug, Ud)

        integer, intent(in)                         :: imax
        real(pr), dimension(0:imax+1,2), intent(in) :: Un !Besoin de revoir les dimension
        real(pr), dimension(0:imax,2), intent(out)  :: Ug
        real(pr), dimension(1:imax+1,2), intent(out):: Ud
        ! -- Local
        !real(pr), dimension(2) :: U_right
        integer :: i

        ! -- On pose que Un(i+2,:)=Un(i+1,:)
        ! Ud(imax+1,:) = Un(imax+1,:)
        ! do i = 1,imax
        !     Ug(i,1) =  1/12.*Un(i-1,1) + 1/3.*Un(i,1) + 7/12.*Un(i+1,1)
        !     Ug(i,2) =  1/12.*Un(i-1,2) + 1/3.*Un(i,2) + 7/12.*Un(i+1,2)
        ! end do
        ! do i = 1,imax
        !     Ud(i,1) = 7/12.*Un(i-1,1) + 1/3.*Un(i,1) + 1/12.*Un(i+1,1)
        !     Ud(i,2) = 7/12.*Un(i-1,2) + 1/3.*Un(i,2) + 1/12.*Un(i+1,2)
        ! end do
        
        do i = 0,imax
            Ug(i,1) = 1/2.*(Un(i+1,1) + Un(i,1))
            Ug(i,2) = 1/2.*(Un(i+1,2) + Un(i,2))
        end do

        do i = 1,imax
            Ud(i,1) = 3/2.*Un(i,1) - 1/2.*Un(i+1,1)
            Ud(i,2) = 3/2.*Un(i,2) - 1/2.*Un(i+1,2)
        end do

    end subroutine


    function Application(sys, imax, tn, Un, bi, dx)result(res)

        integer, intent(in)                  :: imax, sys
        real(pr), intent(in)                 :: tn, dx
        real(pr), dimension(:), intent(in)   :: bi
        real(pr), dimension(0:imax+1,2), intent(in) :: Un

        real(pr), dimension(1:imax,2) :: res

        ! -- Local
        real(pr), dimension(1:imax+1, 2) :: Ud
        real(pr), dimension(0:imax, 2)   :: Ug
        real(pr), dimension(1:imax+1, 2) :: F

        !! -- sys = 1 (test function : x'=-y ; y'=x)
        !! -- sys = 2 (test function : x'=t*x**2)
        !! -- sys = 3 (seconde order discretisation for Saint-Venant)

        Ug(0:imax,:) = Un(0:imax,:)
        Ud(1:imax+1,:) = Un(1:imax+1,:)

        if (sys == 1) then

            res(1,1) = -Un(1,1)
            res(2,1) = Un(0,1)

        else if (sys == 2) then

            res(1,1) = tn*Un(0,1)**2
        
        else if (sys == 3) then

            call RightLeftValues(Un, imax, Ug, Ud)
            F = FluxNum(imax, Ug, Ud, bi)
            res(1:imax,:) = - 1/dx*(F(2:imax+1,:) - F(1:imax,:))

        end if




    end function

end module