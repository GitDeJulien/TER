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


    subroutine sol_approx_tn(flux, dt, cfl, dx)

        ! - Externe
        type(flux_type), intent(inout)       :: flux
        real(pr), intent(in)                 :: dx, cfl 
        real(pr), intent(inout)              :: dt
        

        ! - Interne
        integer                :: i, imax
        real(pr)               :: h, q
        real(pr),dimension(2)  :: F
        real(pr)               :: unsurdx
        real(pr), dimension(size(flux%hnp1)) :: b_i

        unsurdx = 1./dx
        imax = size(flux%hnp1)-2


        do i = 0, imax

            call flux_num(flux%choix_approx_flux, flux%hnp1(i), flux%hnp1(i+1), flux%unp1(i), flux%unp1(i+1), b_i(i+1), F)

            flux%f_h(i+1) = F(1)
            flux%f_q(i+1) = F(2)

        end do

        ! -- Calcule du pas de temps respectant la cfl
        if (dt > dx/(2*MAXVAL(b_i))) then 
            dt = cfl * dx/(2 * MAXVAL(b_i))
        end if

        ! -- Mise à jour des h et q au nouveau pas de temps sans les bords
        do i = 1, imax

            h = flux%hnp1(i) - dt*unsurdx*(flux%f_h(i+1) - flux%f_h(i))
            q = flux%hnp1(i)*flux%unp1(i) - dt*unsurdx*(flux%f_q(i+1) - flux%f_q(i))

            flux%hnp1(i) = h
            flux%unp1(i) = q/h

        end do

        ! condition de Neumann sur les bords
        flux%hnp1(0) = flux%hnp1(1)
        flux%unp1(0) = flux%unp1(1)
        flux%hnp1(imax+1) = flux%hnp1(imax)
        flux%unp1(imax+1) = flux%unp1(imax)

        ! -- Mise à jour du pas de temps
        dt = 0.05


    end subroutine sol_approx_tn


    ! -- Calcul du flux numérique
    subroutine flux_num(choix_flux, hg, hd, &
        ug, ud, bi, F)

        real(pr), intent(in)                :: hg, hd
        real(pr), intent(in)                :: ud, ug
        integer, intent(in)                 :: choix_flux
        real(pr), intent(out)               :: bi
        real(pr), dimension(2), intent(out) :: F
        

        select case(choix_flux)

        case(1)

            bi = max(abs(ug + sqrt(g*hg)), abs(ud + sqrt(g*hd)), abs(ug - sqrt(g*hg)), abs(ud - sqrt(g*hd)))


            F(1) = 0.5*(ud*hd + ug*hg) - 0.5*bi*(hd - hg)
            F(2) = (ug*ug*hg + g*hg*hg/2. + ud*ud*hd + g*hd*hd/2.)*0.5 - 0.5*bi*(hd*ud - hg*ug)


        end select


    end subroutine flux_num



! ------------------------------------------------------





end module