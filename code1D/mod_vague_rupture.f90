module mod_vague_rupture

    use mod_constantes

    implicit none


contains


    ! -- Initialisation des paramètres
    subroutine initialisation(imax, x_i, h_i, u_i, cfl,nb_capteurs,Pos_capteurs,params)


        ! -- Variable externe
        integer, intent(inout)                             :: imax, nb_capteurs
        real(pr), dimension(0:imax+1) , intent(in)     :: x_i
        real(pr), dimension(0:imax+1) , intent(inout)  :: h_i, u_i
        real(pr), dimension(:), allocatable, intent(inout) :: Pos_capteurs
        real(pr), intent(inout)                            :: cfl
        character(len=30),intent(in)                       :: params

        ! -- Variable interne
        real(pr)  :: h_init_am, h_init_av
        integer   :: i, choix_init

        print*,"test 1"        
        open(unit = 2, file = params, action = "read")
        print*,"test 2"        

        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*) h_init_am
        read(2,*) h_init_av
        read(2,*)
        read(2,*) choix_init
        read(2,*) cfl
        read(2,*) nb_capteurs
        print*, " nb capteurs = ",nb_capteurs

        print*, "A la position :"
        if ( nb_capteurs /= 0 ) then
            do i = 0, nb_capteurs-1
                read(2,*) Pos_capteurs(i) 
                print*, "x = ",Pos_capteurs(i) 
            end do
        end if

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

    subroutine maillage(imax, tmax, dx, cfl, x_i, params)
        ! -- Variable externe
        integer, intent(in)                             :: imax
        real(pr), intent(inout)                         :: dx, tmax
        real(pr), intent(inout)                         :: cfl
        real(pr), dimension(0:imax+1) , intent(inout)   :: x_i
        character(len=30),intent(inout)                 :: params


        ! -- Variable interne
        integer   :: i
        real(pr)  :: L, x_min, x_max
    

        
        open(unit = 2, file = params, action = "read")

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