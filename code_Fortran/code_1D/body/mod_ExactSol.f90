module mod_ExactSol

    use mod_constantes
    use mod_initialize
    use mod_ApproxSol
    

    implicit none

contains

    subroutine ExactSolFct(app, tn, x_i, sol_h, sol_u, sigma, params)

        type(approx_type), intent(in)       :: app
        character(len=40),intent(in)        :: params
        real(pr), intent(in)                :: tn
        real(pr), dimension(:), intent(in)  :: x_i
        real(pr), dimension(:), intent(out) :: sol_h, sol_u
        real(pr), intent(out)               :: sigma

        ! -- Variable interne
        integer  :: i, imax
        real(pr) :: lambda_g, lambda_etoile
        real(pr) :: h, u, h_etoile, u_etoile
        real(pr) :: h_init_am, h_init_av

        call hauteur(h_init_am, h_init_av, params)

        imax = size(x_i)-2
        call f_sigma_vp(app, sigma, lambda_g, lambda_etoile, h_etoile, u_etoile, params)

        do i=1, imax+2

            if (x_i(i) < lambda_g*tn) then
                
                sol_h(i) = h_init_am
                sol_u(i) = 0.
                
            else if (x_i(i) >= lambda_g*tn .AND. x_i(i) < lambda_etoile*tn) then

                call detente_h_u(x_i(i), tn, h, u, params)
                sol_h(i) = h
                sol_u(i) = u


            else if (x_i(i) >= lambda_etoile*tn .AND. x_i(i) < sigma*tn) then

                sol_h(i) = h_etoile
                sol_u(i) = u_etoile

            else

                sol_h(i) = h_init_av
                sol_u(i) = 0

            end if

        end do
    end subroutine

    ! ##################################################################

    subroutine f_sigma_vp(app, sigma, lambda_g, lambda_etoile, h_etoile, u_etoile, params)

        type(approx_type), intent(in)  :: app
        character(len=40),intent(in) :: params
        real(pr), intent(out)        :: sigma
        real(pr), intent(out)        :: lambda_g, lambda_etoile
        real(pr), intent (out)       :: h_etoile, u_etoile

        ! -- Variable local
        real(pr) :: fn, fb, hnp1, eps
        integer  :: imax
        real(pr) :: h_init_am, h_init_av

        !imax = size(flux%hnp1)-2
        imax = size(app%Unp1)-2
        eps = 10.**(-5)

        call hauteur(h_init_am, h_init_av, params)

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

    subroutine detente_h_u(x, t, h, u, params)

        real(pr), intent(in)         :: x, t
        character(len=40),intent(in) :: params
        real(pr), intent(out)        :: h, u

        ! -- Variables locales
        real(pr) :: ksi, h_init_am, h_init_av

        call hauteur(h_init_am, h_init_av, params)

        ksi = x/t
        h = 1./(9*g)*(2*sqrt(g*h_init_am)-ksi)**2
        u = 2*(sqrt(g*h_init_am)-sqrt(g*h))

    end subroutine

    ! ####################################################

    subroutine hauteur(h_init_am, h_init_av, params)

        character(len=40),intent(in) :: params
        real(pr), intent(out)        :: h_init_am, h_init_av
        

        open(unit = 2, file = params, action = "read")

        read(2,*)
        read(2,*)
        read(2,*) 
        read(2,*) h_init_am
        read(2,*) h_init_av

        close(2)

    end subroutine


end module