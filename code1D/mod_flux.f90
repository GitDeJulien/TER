module mod_flux

    use mod_constantes
    use mod_vague_rupture

    implicit none

    type flux_type
        !choix du type d'approximation de flux
        integer      :: choix_approx_flux
        !les vecteurs d'entrées h et u (hauteur et vitesse init)
        real(pr), dimension(:), allocatable :: hnp1, unp1
        real(pr), dimension(:), allocatable :: f_h, f_q

        integer      :: n
    end type

contains


    subroutine sol_approx_tn(flux, dt, dx)

        ! - Externe
        type(flux_type), intent(inout)       :: flux
        real(pr), intent(in)                 :: dt, dx
        

        ! - Interne
        integer                :: i, imax
        real(pr)               :: h, q
        real(pr),dimension(2)  :: F
        real(pr)               :: unsurdx

        unsurdx = 1./dx

        imax = size(flux%hnp1)-2

        do i = 1, imax+1

            F = flux_num(flux%choix_approx_flux, flux%hnp1(i-1), flux%hnp1(i), flux%unp1(i-1), flux%unp1(i))

            flux%f_h(i) = F(1)
            flux%f_q(i) = F(2)

        end do

        do i = 1, imax


            h = flux%hnp1(i) - dt*unsurdx*(flux%f_h(i+1) - flux%f_h(i))
            print*, "h=", h
            q = flux%hnp1(i)*flux%unp1(i) - dt*unsurdx*(flux%f_q(i+1) - flux%f_q(i-1))

            flux%hnp1(i) = h
            flux%unp1(i) = q/h

        end do
        print*, ""

        ! condition de Neumann sur les bords
        flux%hnp1(0) = flux%hnp1(1)
        flux%unp1(0) = flux%unp1(1)
        flux%hnp1(imax+1) = flux%hnp1(imax)
        flux%unp1(imax+1) = flux%unp1(imax)


    end subroutine sol_approx_tn


    ! -- Calcul du flux numérique
    function flux_num(choix_flux, hg, hd, &
        ug, ud)result(F)

        real(pr), intent(in)         :: hg, hd
        real(pr), intent(in)         :: ud, ug
        integer, intent(in)          :: choix_flux
        real(pr), dimension(2)       :: F

        real(pr) :: b


        select case(choix_flux)

        case(1)

            b = max(ug + sqrt(g*hg), ud + sqrt(g*hd), ug - sqrt(g*hg), ud - sqrt(g*hd))


            F(1) = 0.5*(ud*hd + ug*hg) - 0.5*b*(hd - hg)
            F(2) = (ug*ug*hg + hg*hg/2. + ud*ud*hd + hd*hd/2.)*0.5 - 0.5*b*(hd*ud - hg*ug)


        end select


    end function flux_num



! ------------------------------------------------------





end module