module mod_error


    use mod_constantes
    use mod_initialize

    implicit none

    contains

    function Error_fct(approx , exa, dx, norme)result(err)

        ! -- Code pour les erreur 
        ! - norme = 1 (L1)
        ! - norme = 2 (L2)
        ! - norme = 3 (Linfini)

        real(pr), dimension(:), intent(in)   :: exa, approx
        integer, intent(in)                  :: norme
        real(pr), intent(in)                 :: dx
        real(pr)                             :: err

        !Interne
        real(pr), dimension(size(exa))  :: err_vect

        err = 0.
        err_vect = abs(exa - approx)

        select case(norme)
        case(1) ! Norme L1

            err = dx*sum(err_vect)
            print*, "Error L1 : ", err

        case(2) ! Norme L2

            err = sqrt(dx*sum(err_vect**2))
            print*, "Error L2 : ", err

        case(3) ! Norme L_infty

            err = maxval(err_vect)

            print*, "Error L_inifnit : ", err

        end select

    end function

end module