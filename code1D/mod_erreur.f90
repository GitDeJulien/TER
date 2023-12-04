module mod_erreur


    use mod_constantes
    use mod_vague_rupture

    implicit none

contains


    function Error_fct(approx, exa, dx, norme)result(err)

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

            err = dx*sum(abs(exa(:)-approx(:)))
            print*, "Erreur L1 = ", err

        case(2) ! Norme L2

            err = sqrt(dx*sum((exa-approx)**2))
            print*, "Error L2 = ", err

        case(3) ! Norme L_infty
            err = maxval(abs(exa(:)-approx(:)))
            print*, "Error L_inifnit = ", err
        end select

    end function


    subroutine OUT(vect_imax, vect_error)

        integer, dimension(:), intent(in)  :: vect_imax
        real(pr), dimension(:), intent(in) :: vect_error

        ! Local variable
        integer :: i

        open(unit = 1, file = 'OUT/error.dat', action = "write")

        do i = 1, size(vect_imax)

            write(1,*) 1./vect_imax(i), vect_error(i)
        end do

        close(1)

    end subroutine
        



end module