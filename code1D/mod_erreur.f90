module mod_erreur


    use mod_constantes
    use mod_vague_rupture

    implicit none

contains


    function Error_fct(approx, exa, norme, dx)result(err)

        ! -- Code pour les erreur 
        ! - norme = 1 (L1)
        ! - norme = 2 (L2)
        ! - norme = 3 (Linfini)

        real(pr), dimension(:), intent(in) :: approx, exa
        integer, intent(in)                :: norme
        real(pr), intent(in)               :: dx
        real(pr)                           :: err

        err = 0.
        select case(norme)
        case(1) ! Norme L1

            err = sum(dx*abs(exa(:)-approx(:)))
            print*, "Erreur L1 = ", err

        case(2) ! Norme L2

            err = sqrt(sum(dx*abs(exa(:)-approx(:))))
            print*, "Error L2 = ", err

        case(3) ! Norme L_infty
            err = maxval(abs(exa(:)-approx(:)))
            print*, "Error L_inifnit = ", err
        end select

    end function




end module