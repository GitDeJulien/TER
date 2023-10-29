module mod_vague_rupture

    use mod_constantes

    implicit none


contains


    ! -- Initialisation des paramètres
    subroutine initialisation(imax, tmax, dx, dt, x_i, h_i, u_i)


        ! -- Variable externe
        integer, intent(inout)                             :: imax
        real(pr), intent(inout)                            :: dx, dt, tmax
        real(pr), dimension(:), allocatable, intent(inout) :: h_i, u_i, x_i

        ! -- Variable interne
        real(pr)  :: h_init_am, h_init_av, L, x_min, x_max
        integer   :: i, choix_init
        real(pr)  :: cfl

        open(unit = 2, file = "params.dat", action = "read")

        read(2,*) imax
        read(2,*) x_min
        read(2,*) x_max
        read(2,*) h_init_am
        read(2,*) h_init_av
        read(2,*) tmax
        read(2,*) choix_init
        read(2,*) cfl

        close(2)

        ! -- Allocation des tableaux
        allocate(h_i(0:imax+2), u_i(0:imax+2))
        allocate(x_i(0:imax+2))

        L = abs(x_max-x_min)
        dx = L/(imax)

        ! -- Calcule de dt avec la condition cfl
        dt = cfl * dx

        select case(choix_init)

        case(1)

            do i = 0, imax+2

                x_i(i) = (i-imax/2)*dx
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





end module