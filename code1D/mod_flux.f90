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


    ! -- Initialisation des paramètres

    function flux_num(choix_flux, hg, hd, &
        ug, ud)result(F)

        real(pr), intent(in)         :: hg, hd
        real(pr), intent(in)         :: ud, ug
        integer, intent(in)          :: choix_flux
        real(pr), dimension(2)       :: F

        real(pr) :: b


        select case(choix_flux)

        case(1)

            ! F(1) = (u_d*h_d - u_g*h_g)*0.5
            ! F(2) = (u_d*u_d*h_d + g*(h_d*h_d)*0.5 - u_g*u_g*h_g + g*(h_g*h_g)*0.5)*0.5
            b = max(ug + sqrt(g*hg), ud + sqrt(g*hd), ug - sqrt(g*hg), ud - sqrt(g*hd))


            F(1) = 0.5*(ud*hd + ug*hg) - 0.5*b*(hd - hg)
            F(2) = (ug*ug*hg + hg*hg/2. + ud*ud*hd + hd*hd/2.)*0.5 - 0.5*b*(hd*ud - hg*ug)


        end select


    end function flux_num

! ------------------------------------------------------





end module