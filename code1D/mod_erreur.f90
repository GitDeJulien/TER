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
        real(pr), intent(in)               :: dx
        integer, intent(in)                :: norme

        real(pr)                           :: err

        ! -- Locale
        integer :: i

        err = 0.
        select case(norme)
        case(1) ! Norme L1
            do i=1,size(approx)
                err = err + dx*abs(exa(i)-approx(i))
            end do

        case(2)
            do i=1,size(approx)
                err = err + dx*abs(exa(i)-approx(i))**2
            end do
            err = sqrt(err)

        case(3)
            err = maxval(abs(exa(:)-approx(:)))
        end select

    end function


end module