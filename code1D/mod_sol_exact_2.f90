module mod_sol_exact_2

    use mod_constantes
    use mod_flux
    use mod_vague_rupture

    implicit none



contains

    subroutine sol_exact_tn(flux, tn, x_i, sol_h, sol_u)

        type(flux_type), intent(in)         :: flux
        real(pr), intent(in)                :: tn
        real(pr), dimension(:), intent(in)  :: x_i
        real(pr), dimension(:), intent(out) :: sol_h, sol_u

        ! -- Variable interne
        integer  :: i, imax
        real(pr) :: sigma, lambda_g, lambda_etoile
        real(pr) :: h, u
        real(pr) :: h_init_am, h_init_av

        call hauteur(h_init_am, h_init_av)

        imax = size(x_i)-2
        call f_sigma_vp(flux, sigma, lambda_g, lambda_etoile )

        do i=1, imax+2

            if (x_i(i) < lambda_g*tn) then
                
                sol_h(i) = h_init_am
                sol_u(i) = 0.
                
            else if (x_i(i) >= lambda_g*tn .AND. x_i(i) < lambda_etoile*tn) then

                call detente_h_u(x_i(i), tn, h, u)
                sol_h(i) = h
                sol_u(i) = u

            else if (x_i(i) >= lambda_etoile*tn .AND. x_i(i) < sigma*tn) then

                sol_h(i) = h
                sol_u(i) = u

            else

                sol_h(i) = h_init_av
                sol_u(i) = u

            end if

        end do
    end subroutine

    ! ##################################################################

    subroutine f_sigma_vp(flux, sigma, lambda_g, lambda_etoile)

        type(flux_type), intent(in)  :: flux
        real(pr), intent(out)        :: sigma
        real(pr), intent(out)        :: lambda_g, lambda_etoile

        ! -- Variable local
        !real(pr), dimension(size(flux%hnp1)) :: xn
        real(pr) :: fn, fb, h_etoile, u_etoile, hnp1, eps
        integer  :: imax
        real(pr) :: h_init_am, h_init_av

        imax = size(flux%hnp1)-2
        eps = 10.**(-5)

        call hauteur(h_init_am, h_init_av)

        h_etoile = 0.
        lambda_g = -sqrt(g*h_init_am)
        hnp1 = h_init_am

        do while (abs(hnp1-h_etoile) > eps)
            
            h_etoile = hnp1

            fn = 2*(sqrt(g*h_init_am)-sqrt(g*h_etoile))+(h_init_av-h_etoile)*sqrt(g*(h_etoile+h_init_av)/(2*h_etoile*h_init_av))
            fb = 2*(sqrt(g*h_init_am)-sqrt(g*h_init_av))

            hnp1 = h_etoile - (h_init_av-h_etoile)/(fb-fn)*fn

        end do
        
        h_etoile = hnp1

        u_etoile = 2*(sqrt(g*h_init_am)-sqrt(g*h_etoile))
        sigma = h_etoile*u_etoile/(h_etoile-h_init_av)
        lambda_etoile = u_etoile - sqrt(g*h_etoile)

    end subroutine

    ! ###################################################

    subroutine detente_h_u(x, t, h, u)

        real(pr), intent(in)  :: x, t
        real(pr), intent(out) :: h, u

        ! -- Variables locales
        real(pr) :: ksi, h_init_am, h_init_av

        call hauteur(h_init_am, h_init_av)

        ksi = x/t
        h = 1./(9*g)*(2*sqrt(g*h_init_am)-ksi)**2
        u = 2*(sqrt(g*h_init_am)-sqrt(g*h))

    end subroutine

    ! ####################################################

    subroutine hauteur(h_init_am, h_init_av)

        real(pr), intent(out) :: h_init_am, h_init_av

        open(unit = 2, file = "params.dat", action = "read")

        read(2,*)
        read(2,*)
        read(2,*) 
        read(2,*) h_init_am
        read(2,*) h_init_av

        close(2)

    end subroutine


end module