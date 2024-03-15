module mod_initialize

    use mod_constantes

    implicit none

contains

    ! -- Initialisation des paramètres
    subroutine initialize(imax, tmax, dx, x_i, U0, cfl, params)


        ! -- Variable externe
        integer, intent(inout)                          :: imax
        real(pr), intent(inout)                         :: cfl, tmax
        real(pr), dimension(0:imax+1), intent(out)      :: x_i
        real(pr), dimension(0:imax+1,2), intent(out)    :: U0
        real(pr), intent(out)                           :: dx
        character(len=30), intent(in)                   :: params

        ! -- Variable interne
        real(pr)  :: HInitLeft, HInitRight
        integer   :: i, ChoixInit
        real(pr)  :: L, x_min, x_max


        open(unit = 2, file = params, action = "read")
        

        read(2,*)
        read(2,*) x_min
        read(2,*) x_max
        read(2,*) HInitLeft
        read(2,*) HInitRight
        read(2,*) tmax
        read(2,*) ChoixInit
        read(2,*) cfl

        close(2)

        print*, "-- Parameters readed --"
        
        ! -- Création du maillage 1D
        L = abs(x_max-x_min)
        dx = L/(imax)

        do i = 0, imax+1

            x_i(i) = (i-imax/2)*dx

        end do

        select case(ChoixInit)

        case(1)

            do i = 0, imax+1
                U0(i,2) = 0._pr
                if (x_i(i) < 0) then
                    U0(i,1) = HInitLeft
                else
                    U0(i,1) = HInitRight
                end if
            end do

        end select

    end subroutine initialize


end module