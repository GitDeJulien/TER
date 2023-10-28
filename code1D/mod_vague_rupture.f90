module mod_vague_rupture

    use mod_constantes

    implicit none

    type vague_type
        !choix du type de condition initial
        integer      :: choix_cond_init
        integer      :: n
        real(pr)     :: x_max, x_min !abscisse min et max
        real(pr)     :: d !profondeur eau
        real(pr)     :: Delta
        real(pr)     :: Hauteur_init
        !les vecteurs d'entrées h0 et u0 (hauteur et débit init)
        real(pr), dimension(:), allocatable :: hi_0, ui_0
        real(pr), dimension(:), allocatable :: x_i

    end type

contains


    ! -- Initialisation des paramètres

    subroutine init_vague_rupture(vague, L, dx)

        type(vague_type), intent(out) :: vague
        real(pr), intent(out)        :: L, dx

        integer :: i
        open(unit = 1, file = "params_vague_rupture.dat", action = "read")


        read(1,*) vague%n
        read(1,*) vague%x_min
        read(1,*) vague%x_max
        read(1,*) vague%d
        read(1,*) vague%Delta
        read(1,*) vague%Hauteur_init

        close(1)

        ! -- Demander à l'utilisateur de choisir la cond_init
        print*, "Entrer  votre choix de condition initiale"
        print*, "1) Front en crénau"
        read*, vague%choix_cond_init

        L = abs(vague%x_max-vague%x_min)
        dx = L/(vague%n)


        allocate(vague%hi_0(0:vague%n+2), vague%ui_0(0:vague%n+2))
        allocate(vague%x_i(0:vague%n+2))

        select case(vague%choix_cond_init)

        case(1)

            do i = 0, vague%n+1

                vague%x_i(i) = (i-vague%n/2)*dx

                vague%ui_0(i) = 0._pr
                if (vague%x_i(i)<0) then
                    vague%hi_0(i) = 1._pr
                else
                    vague%hi_0(i) = 0._pr
                end if

            end do


        case(2)

        do i = 1, vague%n+1

            vague%x_i(i) = (i-vague%n/2)*dx
            vague%hi_0(i) = vague%d + (16._pr/6) * (vague%d**3/vague%Delta**2) * 1/COSH(vague%x_i(i)/vague%Delta)**2
            vague%ui_0(i) = 0._pr

        end do

        end select



    end subroutine init_vague_rupture





end module