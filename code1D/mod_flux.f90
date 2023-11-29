module mod_flux

    use mod_constantes
    use mod_vague_rupture

    implicit none

    type flux_type
        !choix du type d'approximation de flux
        integer      :: choix_approx_flux
        !les vecteurs d'entrées h et u (hauteur et vitesse init)
        real(pr), dimension(:), allocatable :: hnp1, unp1 ! //Mieux vaux appeler un veteur Un avec hn et un
        real(pr), dimension(:), allocatable :: f_h, f_q

    end type

contains


    subroutine sol_approx_tn(flux, dt, cfl, dx, tn)

        ! - Externe
        type(flux_type), intent(inout)       :: flux
        real(pr), intent(in)                 :: dx, cfl
        real(pr), intent(inout)              :: dt, tn
        

        ! - Interne
        integer                                     :: i, imax
        real(pr), dimension(1:size(flux%hnp1)-1, 2) :: Ug, Ud, F  ! - Taille imax+1
        real(pr), dimension(1:size(flux%hnp1)-2, 2) :: Un, Unp1   ! - Taille imax
        real(pr)                                    :: unsurdx
        real(pr), dimension(1:size(flux%hnp1)-1)    :: b_i


        unsurdx = 1./dx
        imax = size(flux%hnp1)-2

        ! -- Tableau des bi
        do i=0, imax
            b_i(i+1) = f_bi(flux%hnp1(i), flux%hnp1(i+1), flux%unp1(i), flux%unp1(i+1))
        end do


        ! -- Calcule du pas de temps respectant la cfl
        if (dt > dx/(2*MAXVAL(b_i))) then 
            dt = cfl * dx/(2 * MAXVAL(b_i))
        end if

        select case(flux%choix_approx_flux)

        case(1)

            Ug(1:imax+1,1) = flux%hnp1(0:imax)
            Ug(1:imax+1,2) = flux%hnp1(0:imax)*flux%unp1(0:imax)
            Ud(1:imax+1,1) = flux%hnp1(1:imax+1)
            Ud(1:imax+1,2) = flux%hnp1(1:imax+1)*flux%unp1(1:imax+1)

            F = flux_num_2(imax, Ug, Ud, b_i)

            Unp1(1:imax,1) = flux%hnp1(1:imax) - dt*unsurdx*(F(2:imax+1,1) - F(1:imax,1))
            Unp1(1:imax,2) = flux%hnp1(1:imax)*flux%unp1(1:imax) - dt*unsurdx*(F(2:imax+1,2) - F(1:imax, 2))

        case(2) ! - Ordre 2 stencile (Ui-1, Ui, Ui+1)
            

            Ug(1:imax,1) = 1./3*(flux%hnp1(0:imax-1) + flux%hnp1(1:imax) + flux%hnp1(2:imax+1))
            Ug(1:imax,2) = 1./2*(flux%hnp1(0:imax-1)*flux%unp1(0:imax-1) + flux%hnp1(1:imax)*flux%unp1(1:imax) + &
            flux%hnp1(2:imax+1)*flux%unp1(2:imax+1))
            Ud(1:imax,1) = 5./6*flux%hnp1(2:imax+1) + 1./3*flux%hnp1(1:imax) - 1./6*flux%hnp1(0:imax-1)
            Ud(1:imax,2) = 5./6*flux%hnp1(2:imax+1)*flux%unp1(2:imax+1) + 1./3*flux%hnp1(1:imax)*flux%unp1(1:imax) - &
            1./6*flux%hnp1(0:imax-1)*flux%unp1(0:imax-1)

            Un(1:imax,1) = flux%hnp1(1:imax)
            Un(1:imax,2) = flux%hnp1(1:imax)*flux%unp1(1:imax)

            Unp1 = RK2_SSP(flux, Un, Ug, Ud, dt, unsurdx, b_i, imax)

        end select


        flux%hnp1(1:imax) = Unp1(1:imax,1)
        flux%unp1(1:imax) = Unp1(1:imax,2)/Unp1(1:imax,1)

        ! -- Conditions de Neumann au bord
        flux%hnp1(0) = flux%hnp1(1)
        flux%unp1(0) = flux%unp1(1)
        flux%hnp1(imax+1) = flux%hnp1(imax)
        flux%unp1(imax+1) = flux%unp1(imax)

        tn = tn + dt


    end subroutine sol_approx_tn

! ###########################################################################################

    function flux_num_2(imax, Ug, Ud, bi)result(F)

        integer, intent(in)                          :: imax
        real(pr), dimension(1:imax+1, 2), intent(in) :: Ug, Ud
        real(pr), dimension(imax+1), intent(in)      :: bi

        real(pr), dimension(1:imax+1, 2)             :: F

        ! -- Local
        integer :: i

        do i=1, imax+1

            F(i,1) = 0.5*(Ud(i,2) + Ug(i,2)) - 0.5*bi(i)*(Ud(i,1) - Ug(i,1))
            F(i,2) = (Ug(i,2)*Ug(i,2)/Ug(i,1) + g*Ug(i,1)*Ug(i,1)/2. + Ud(i,2)*Ud(i,2)/Ud(i,1)&
            + g*Ud(i,1)*Ud(i,1)/2.)*0.5 - 0.5*bi(i)*(Ud(i,2) - Ug(i,2))

        end do

    end function

! #############################################################################################

    ! -- Clacul du coefficient bi
    function f_bi(hg, hd, ug, ud)result(bi)
        real(pr), intent(in)                :: hg, hd
        real(pr), intent(in)                :: ud, ug
        real(pr)                            :: bi

        bi = max(abs(ug + sqrt(g*hg)), abs(ud + sqrt(g*hd)), abs(ug - sqrt(g*hg)), abs(ud - sqrt(g*hd)))

    end function f_bi


! #############################################################################################

    ! -- RK2 - Heun - SSP
    function RK2_SSP(flux, Un, Ug, Ud, dt, unsurdx, bi, imax)result(Unp1)

        type(flux_type), intent(in)          :: flux
        real(pr), intent(in)                 :: dt, unsurdx
        real(pr), dimension(:,:), intent(in) :: Ug, Ud
        real(pr), dimension(:,:), intent(in) :: Un
        real(pr), dimension(:), intent(in)   :: bi
        integer, intent(in)                  :: imax

        real(pr), dimension(1:imax, 2)   :: Unp1

        ! -- Local
        real(pr), dimension(1:imax, 2)   :: k1, k2
        real(pr), dimension(1:imax+1, 2) :: F
        
        F = flux_num_2(imax, Ug, Ud, bi)
        k1(1:imax, :) = -dt*unsurdx*(F(2:imax+1,:) - F(1:imax,:))
        F = flux_num_2(imax, Un+dt*k1, Un+dt*k1, bi)
        k2(1:imax, :) = -dt*unsurdx*(F(2:imax+1,:) - F(1:imax,:))

        Unp1(1:imax,1) = flux%hnp1(1:imax) + dt*(1./2*k1(1:imax, 1) + 1./2*k2(1:imax, 1))
        Unp1(1:imax,2) = flux%hnp1(1:imax)*flux%unp1(1:imax) + dt*(1./2*k1(1:imax, 2) + 1./2*k2(1:imax, 2)) 


    end function




end module