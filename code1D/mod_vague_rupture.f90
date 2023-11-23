module mod_vague_rupture

    use mod_constantes

    implicit none


contains


    ! -- Initialisation des paramètres
    subroutine initialisation(imax, x_i, h_i, u_i)

        ! -- Variable externe
        integer, intent(in)                     :: imax
        real(pr), dimension(0:imax+1) , intent(in)     :: x_i
        real(pr), dimension(0:imax+1) , intent(inout)  :: h_i, u_i

        ! -- Variable interne
        real(pr)  :: h_init_am, h_init_av
        integer   :: i, choix_init

        open(unit = 2, file = "params.dat", action = "read")

        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*) h_init_am
        read(2,*) h_init_av
        read(2,*)
        read(2,*) choix_init

        close(2)

        select case(choix_init)

        case(1)

            do i = 0, imax+1

                u_i(i) = 0._pr

                if (x_i(i) < 0) then
                    h_i(i) = h_init_am
                else
                    h_i(i) = h_init_av
                end if

            end do

        case(2)

        end select

    end subroutine initialisation

    subroutine maillage(imax, tmax, dx, cfl, x_i)
        ! -- Variable externe
        integer, intent(in)                    :: imax
        real(pr), intent(inout)                :: dx, tmax
        real(pr), intent(inout)                :: cfl
        real(pr), dimension(0:imax+1) , intent(inout) :: x_i

        ! -- Variable interne
        integer   :: i
        real(pr)  :: L, x_min, x_max

        open(unit = 2, file = "params.dat", action = "read")

        read(2,*)
        read(2,*) x_min
        read(2,*) x_max
        read(2,*) 
        read(2,*) 
        read(2,*) tmax
        read(2,*) 
        read(2,*) cfl

        close(2)


        L = abs(x_max-x_min)
        dx = L/(imax)

        do i = 0, imax+1

            x_i(i) = (i-imax/2)*dx

        end do

    end subroutine



end module