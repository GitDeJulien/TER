module mod_flux

    use mod_constantes
    use mod_vague_rupture

    implicit none

    type flux_type
        !choix du type d'approximation de flux
        integer      :: choix_approx_flux
        !les vecteurs d'entrées h et u (hauteur et vitesse init)
        real(pr), dimension(:), allocatable :: hnp1, unp1 ! //Mieux vaux appeler un veteur Un avec hn et un
        real(pr), dimension(:), allocatable :: f_h, f_q

    end type

contains


    subroutine sol_approx_tn(flux, dt, cfl, dx, tn)

        ! - Externe
        type(flux_type), intent(inout)       :: flux
        real(pr), intent(in)                 :: dx, cfl
        real(pr), intent(inout)              :: dt, tn
        

        ! - Interne
        integer                                     :: i, imax
        real(pr), dimension(1:size(flux%hnp1)-1, 2) :: Ud  ! - Taille imax+1 (1:imax+1)
        real(pr), dimension(0:size(flux%hnp1)-2, 2) :: Ug      ! - Taille imax+1  (0:imax)
        real(pr), dimension(1:size(flux%hnp1)-1, 2) :: F          ! - Taille imax
        real(pr), dimension(0:size(flux%hnp1)-1, 2) :: Un, Unp1        ! - Taille imax+2 (0:imax+1)
        real(pr)                                    :: unsurdx
        real(pr), dimension(1:size(flux%hnp1)-1)    :: b_i


        unsurdx = 1./dx
        imax = size(flux%hnp1)-2

        ! -- Tableau des bi
        do i=0, imax
            b_i(i+1) = f_bi(flux%hnp1(i), flux%hnp1(i+1), flux%unp1(i), flux%unp1(i+1))
        end do


        ! -- Calcule du pas de temps respectant la cfl
        if (dt > dx/(2*MAXVAL(b_i))) then 
            dt = cfl * dx/(2 * MAXVAL(b_i))
        end if

        Un(0:imax+1,1) = flux%hnp1(0:imax+1)
        Un(0:imax+1,2) = flux%hnp1(0:imax+1)*flux%unp1(0:imax+1)

        select case(flux%choix_approx_flux)

        case(1)

            Ug(0:imax,:) = Un(0:imax,:)
            Ud(1:imax+1,:) = Un(1:imax+1,:)

            !call discretisation_second_order(imax, Un, Ug, Ud)

            F = flux_num_2(imax, Ug, Ud, b_i)

            ! -- Euler sur le temps
            Unp1(1:imax,1) = flux%hnp1(1:imax) - dt*unsurdx*(F(2:imax+1,1) - F(1:imax,1))
            Unp1(1:imax,2) = flux%hnp1(1:imax)*flux%unp1(1:imax) - dt*unsurdx*(F(2:imax+1,2) - F(1:imax, 2))

            ! -- Condition de bord - Neumann
            Unp1(0,:) = Unp1(1,:)
            Unp1(imax+1,:) = Unp1(imax,:)

        case(2) ! - Ordre 2 stencile (Ui-1, Ui, Ui+1)

            Unp1 = RK2_SSP(Un, Ug, Ud, dt, unsurdx, b_i, imax)
            call test(flux, imax, unsurdx, dt, b_i, Ug, Ud, Un, Unp1)

        end select


        flux%hnp1(0:imax+1) = Unp1(0:imax+1,1)
        flux%unp1(0:imax+1) = Unp1(0:imax+1,2)/Unp1(0:imax+1,1)

        tn = tn + dt

    end subroutine sol_approx_tn

! ###########################################################################################

    function flux_num_2(imax, Ug, Ud, bi)result(F)

        integer, intent(in)                          :: imax
        real(pr), dimension(1:imax+1, 2), intent(in) :: Ud
        real(pr), dimension(0:imax, 2), intent(in)   :: Ug
        real(pr), dimension(imax+1), intent(in)      :: bi

        real(pr), dimension(1:imax+1, 2)             :: F


        F(1:imax+1,1) = 0.5*(Ud(1:imax+1,2) + Ug(0:imax,2)) - 0.5*bi(1:imax+1)*(Ud(1:imax+1,1) - Ug(0:imax,1))
        F(1:imax+1,2) = (Ug(0:imax,2)*Ug(0:imax,2)/Ug(0:imax,1) + &
        g*Ug(0:imax,1)*Ug(0:imax,1)/2. + Ud(1:imax+1,2)*Ud(1:imax+1,2)/Ud(1:imax+1,1)&
        + g*Ud(1:imax+1,1)*Ud(1:imax+1,1)/2.)*0.5 - 0.5*bi(1:imax+1)*(Ud(1:imax+1,2) - Ug(0:imax,2))

    end function

! #############################################################################################

    ! -- Clacul du coefficient bi
    function f_bi(hg, hd, ug, ud)result(bi)
        real(pr), intent(in)                :: hg, hd
        real(pr), intent(in)                :: ud, ug
        real(pr)                            :: bi

        bi = max(abs(ug + sqrt(g*hg)), abs(ud + sqrt(g*hd)), abs(ug - sqrt(g*hg)), abs(ud - sqrt(g*hd)))

    end function f_bi


! #############################################################################################

    ! -- RK2 - Heun - SSP
    function RK2_SSP(Un, Ug, Ud, dt, unsurdx, bi, imax)result(Unp1)

        real(pr), intent(in)                           :: dt, unsurdx
        integer, intent(in)                            :: imax
        real(pr), dimension(0:imax,2), intent(inout)   :: Ug
        real(pr), dimension(1:imax+1,2), intent(inout) :: Ud
        real(pr), dimension(0:imax+1,2), intent(in)    :: Un
        real(pr), dimension(:), intent(in)             :: bi

        real(pr), dimension(0:imax+1, 2)   :: Unp1

        ! -- Local
        real(pr), dimension(0:imax+1, 2) :: k1, k2, UnStar
        real(pr), dimension(1:imax+1, 2) :: F


        call discretisation_second_order(imax, Un, Ug, Ud)
        F = flux_num_2(imax, Ug, Ud, bi)
        !k1(1:imax, :) = -dt*unsurdx*(F(2:imax+1,:) - F(1:imax,:))
        k1(1:imax, :) = -unsurdx*(F(2:imax+1,:) - F(1:imax,:))

        Unstar(1:imax,:) = Un(1:imax,:) + dt*k1(1:imax,:)
        UnStar(0,:) = Un(0,:)
        UnStar(imax+1,:) = Un(imax+1,:)

        call discretisation_second_order(imax, UnStar, Ug, Ud)
        F = flux_num_2(imax, Ug, Ud, bi)
        !k2(1:imax, :) = -dt*unsurdx*(F(2:imax+1,:) - F(1:imax,:))
        k2(1:imax, :) = -unsurdx*(F(2:imax+1,:) - F(1:imax,:))

        Unp1(1:imax, :) = Un(1:imax, :) + dt*(1./2*k1(1:imax, :) + 1./2*k2(1:imax, :))

        Unp1(0,:) = Unp1(1,:)
        Unp1(imax+1,:) = Unp1(imax,:)

    end function

    subroutine discretisation_second_order(imax, Un, Ug, Ud)
        
        integer, intent(in)                   :: imax
        real(pr), dimension(0:imax+1,2), intent(in)  :: Un
        real(pr), dimension(0:imax,2), intent(out)   :: Ug
        real(pr), dimension(1:imax+1,2), intent(out) :: Ud

        !Ug(1:imax,:) = 1./3*(Un(0:imax-1,:) + Un(1:imax,:) + Un(2:imax+1,:))
        Ug(1:imax,:) = 5./6*Un(0:imax-1,:) + 1./3*Un(1:imax,:) - 1./6*Un(2:imax+1,:)
        !Ud(1:imax,:) = 5./6*Un(2:imax+1,:) + 1./3*Un(1:imax,:) - 1./6*Un(0:imax-1,:)
        Ud(1:imax,:) = 1./3*(Un(0:imax-1,:) + Un(1:imax,:) + Un(2:imax+1,:))

        Ug(0,:) = Un(0,:)
        Ud(imax+1,:) = Un(imax+1, :)

    end subroutine

    subroutine test(flux, imax, unsurdx, dt, bi, Ug, Ud, Un, Unp1)

        type(flux_type)     :: flux
        integer, intent(in) :: imax
        real(pr), intent(in):: unsurdx, dt
        real(pr), dimension(:), intent(in)             :: bi
        real(pr), dimension(0:imax,2), intent(inout)   :: Ug
        real(pr), dimension(1:imax+1,2), intent(inout) :: Ud
        real(pr), dimension(0:imax+1,2), intent(inout) :: Un
        real(pr), dimension(0:imax+1, 2), intent(inout) :: Unp1

        ! -- Local
        integer :: i
        real(pr), dimension(1:size(flux%hnp1)-1, 2) :: F
        integer :: valid
        integer :: cpt

        valid = 1
        cpt = 0

        do i=1,imax+1
            if (Unp1(i, 1)<1) then

                valid = 0
                
            end if
        end do

        if (valid == 0) then
            Un(0:imax+1,1) = flux%hnp1(0:imax+1)
            Un(0:imax+1,2) = flux%hnp1(0:imax+1)*flux%unp1(0:imax+1)

            Ug(0:imax,:) = Un(0:imax,:)
            Ud(1:imax+1,:) = Un(1:imax+1,:)

            F = flux_num_2(imax, Ug, Ud, bi)

            ! -- Euler sur le temps
            Unp1(1:imax,1) = flux%hnp1(1:imax) - dt*unsurdx*(F(2:imax+1,1) - F(1:imax,1))
            Unp1(1:imax,2) = flux%hnp1(1:imax)*flux%unp1(1:imax) - dt*unsurdx*(F(2:imax+1,2) - F(1:imax, 2))

            ! -- Condition de bord - Neumann
            Unp1(0,:) = Unp1(1,:)
            Unp1(imax+1,:) = Unp1(imax,:)

        else
            cpt = cpt + 1
            print*, "cpt =", cpt

        end if

    end subroutine

    subroutine evolution_capteur(name_file,imax,iter,dt,nb_capteurs,Pos_capteurs)
        
        !-----Variables d'entrée-----!
        integer, intent(in)             :: imax, iter, nb_capteurs
        real(pr),dimension(0:nb_capteurs-1),intent(in) :: Pos_capteurs
        character(len=30),intent(in)    :: name_file
        real(PR),intent(in)             :: dt

        !-----Variables de sortie-----!
        character(len=30)            :: name_file_out
        !character(len=10)            :: xch
        character(len=3)             :: ich

        !-----Variables locales-----!
        integer,dimension(:),allocatable  :: Pos_maillage
        real,dimension(:),allocatable     :: H
        real(pr)                          :: a, xim1, t!, x
        integer                           :: k, km1, Nmax, j, l, iostatus, i

        allocate(Pos_maillage(0:nb_capteurs-1),H(0:nb_capteurs-1))

        name_file_out = "OUT/capteur_"//trim(ich)//".txt"

        open(10,file=name_file,action="read")

        Nmax=10000
        a=-40000             !Par sécurité
        xim1=0
        k=0

        !---------- Recherche du point du maillage le plus proche pour chaque capteurs ------!

        do i = 0,nb_capteurs-1
        !--- On considère que Pos_capteurs est classé de manière croissante ---!
            print*, "position du capteur x = ", Pos_capteurs(i)
            do while (Pos_capteurs(i)>a ) 
                xim1=a                          ! On stocke le point précédent !
                read(10,*) a                    ! On lit le point suivant !
                k = k + 1                       ! k compte la position du  point du maillage !
            end do
            Pos_maillage(i) = k
            !print*, "Pour le capteur numéro ", i, ", le point de maillagee le plus proche est le numéro ", k
            !print*, "en x = ", a
            
        end do

        close(10)

        ! ---- Choix arbitraire de prendre le point juste d'après ----!

        ! ---- Ecriture de l'évolution des différents capteurs en fonction  

        
        !!name_file_out = "OUT/capteur_"//trim(ich)//"_x_"//trim(xch)//".txt"
        name_file_out = "OUT/capteur.txt"
        open(11,file=name_file_out,action="write")

        open(10,file=name_file,action="read")

        t = 0._PR
        
        do l = 1,iter
            km1 = 0
            do i = 0, nb_capteurs -1
                k = Pos_maillage(i)

                do j = 1, k-1-km1
                    read(10,*,iostat=iostatus)
                end do
                km1 = k
                read(10,*,iostat=iostatus) a, H(i)
                
                
            end do

            write(11, '(F0.3,F12.6,F12.6,F12.6)') t, (H(j), j=0, nb_capteurs-1)

            do j = k+1, imax + 3
                read(10,*,iostat=iostatus)
                if ( iostatus/=0 ) then
                    print*, "Fin du fichier"
                    exit
                end if
            end do

            t=t+dt
            
        end do
        
        close(11)
        close(10)

        print*, "Le nom du fichier du capteur est ",name_file_out
    end subroutine evolution_capteur


end module