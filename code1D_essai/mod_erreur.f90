module mod_erreur


    use mod_constantes
    use mod_vague_rupture

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
        err_vect(:) = abs(exa - approx)

        select case(norme)
        case(1) ! Norme L1

            err = dx*sum(err_vect)
            print*, "Erreur L1 = ", err

        case(2) ! Norme L2

            err = sqrt(dx*sum(err_vect(:)**2))
            print*, "Error L2 = ", err

        case(3) ! Norme L_infty

            err = maxval(err_vect)

            print*, "Error L_inifnit = ", err

        case(4) ! Norme H1

            err = sqrt(sum((err_vect(2:size(exa))-err_vect(1:size(exa)-1))**2)/dx)
            print*, "Erreur H1 =", err

        end select

    end function


    subroutine OUT(vect_imax, vect_error)

        integer, dimension(:), intent(in)  :: vect_imax
        real(pr), dimension(:), intent(in) :: vect_error

        ! Local variable
        integer :: i

        open(unit = 1, file = 'OUT/error.dat', action = "write")

        do i = 1, size(vect_imax)

            write(1,*) 1./(vect_imax(i)), vect_error(i)
        end do

        close(1)

    end subroutine
        



end module