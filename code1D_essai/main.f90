program main

    use mod_vague_rupture
    use mod_constantes
    use mod_erreur
    use mod_flux
    use mod_sol_exact_2

    implicit none


    
    type(flux_type)                     :: flux
    real(pr)                            :: tn, dx, tmax, dt, cfl
    real(pr), dimension(:), allocatable :: x_i, h_i, u_i
    real(pr), dimension(:), allocatable :: sol_exa_h, sol_exa_u, Pos_capteurs
    real(pr)                            :: errorL1, errorL2, errorLi, errorH1
    real(pr), dimension(2)              :: Vect_err
    real(pr)                            :: sigma, lambda_g, lambda_etoile, h_etoile, u_etoile

    character(len=30)                   :: name_file,params,date,num

    integer, dimension(2)              :: Vect_imax
    integer :: i, iter, Nmax, imax, cpt, nb_capteurs

    ! -- Sécurité
    Nmax = 10**6
    ! -- Initialisation d'un compteur
    cpt = 0

    print*, "Quel est la date de l'expérience ? (JJ_MM)"
    read(*,*) date
    print*, "Quel est le numéro de l'expérience ? "
    read(*,*) num

    params = "exp_"//trim(date)//"/param/exp_"//trim(num)//".dat"
    !params = "params.dat"

    open(unit = 2, file = params, action = "read")
        read(2,*) imax
        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*) 
        read(2,*) 
        read(2,*)
        read(2,*) 
        read(2,*) nb_capteurs
    close(2)

    allocate(Pos_capteurs(0:nb_capteurs-1))


    print*, "Entrer votre choix d'approximation du flux : "
    print*, " 1) Approximation de Rusanov d'ordre 1"
    print*, " 2) Approximation de Van Leer d'ordre 2"
    read*, flux%choix_approx_flux

    10 continue

    ! -- Ouverture du fichier d'écriture des résultats d'approximation
    if (flux%choix_approx_flux == 1) then
        name_file = "OUT/approx_RUSANOV.dat"
        open(unit = 10, file = "OUT/approx_RUSANOV.dat", action = "write")
    else if (flux%choix_approx_flux == 2) then
        name_file = "OUT/approx_Van_Leer.dat"
        open(unit = 10, file = "OUT/approx_Van_Leer.dat", action = "write")
    else
        print*,"Le choix de flux que vous avez fait n'est pas valide, veuillez recommencer"
        stop
    end if

    print*, "Initialisation ..."
    ! -- Allocation des tableaux de flux de hauteur et de vitesse
    allocate(x_i(0:imax+1))
    allocate(h_i(0:imax+1), u_i(0:imax+1))
    !allocate(flux%f_h(0:imax+1), flux%f_q(0:imax+1))
    !allocate(flux%hnp1(0:imax+1), flux%unp1(0:imax+1))
    allocate(flux%Unp1(0:imax+1, 2))
    allocate(sol_exa_h(0:imax+1), sol_exa_u(0:imax+1))


    print*,"imax = ",imax

    ! -- Initialisation maillage
    print*,"Calcul maillage..."
    call maillage(imax, tmax, dx, cfl, x_i, params)
    print*,"Maillage calculé."
    ! initialisation maillage + condition initiale
    call initialisation(imax, x_i, h_i, u_i, cfl,nb_capteurs,Pos_capteurs, params)

    ! -- Initialisation des flux 
    ! flux%hnp1 = h_i
    ! flux%unp1 = u_i
    flux%Unp1(:, 1) = h_i(:)
    flux%Unp1(:, 2) = u_i(:)*h_i(:) 

    ! -- Écriture de la condition initiale
    do i = 0, imax+1
        write(10,*) x_i(i), h_i(i)
    end do
    write(10,*)
    write(10,*)

    dt = 0.01
    tn = 0.
    iter = 0
    sol_exa_h = 0.
    sol_exa_u = 0.

    ! -- Ouverture du fichier d'écriture des résultats
    open(unit = 11, file = "OUT/sol_exact.dat", action = "write")
    open(unit = 12, file = "OUT/sol_h_u.dat", action = "write")

    print*, "Passage n°", cpt+1
    print*, "-----------------------------------------"
    print*, "dx =", dx
    print*, "tmax =", tmax

    ! -- Boucle en temps
    print*, "-----------------------------------------"
    print*,"Calcul de la solution..."
    ! -- Boucle en temps
    do while (tn <= tmax .AND. iter < Nmax)

        iter = iter + 1
        ! -- Solution exacte
        call sol_exact_tn(flux, tn, x_i, sol_exa_h, sol_exa_u, params)

        ! -- Solution approchée
        call sol_approx_tn(flux, dt, cfl, dx, tn, imax, iter)

        ! -- Mise à jour des flux
        ! h_i = flux%hnp1
        ! u_i = flux%unp1
        h_i = flux%Unp1(:,1)
        u_i = flux%Unp1(:,2)/h_i(:)

        ! -- Écriture dans les fichiers un fichier .dat

        do i = 0, imax+1
            write(10,*) x_i(i), h_i(i)
        end do

        write(10,*)
        write(10,*)

        if (iter == 100) then
            print*, "dt =", dt
            print*,"tn =" , tn
            do i = 0, imax+1
                write(12,*) x_i(i), h_i(i), u_i(i), sol_exa_h(i), sol_exa_u(i)
            end do
        end if

        do i = 0, imax+1
            write(11,*) x_i(i), sol_exa_h(i)
        end do

        write(11,*)
        write(11,*)

    end do

    ! -- Calcul de la vitesse de l'onde de choc
    print*, "----------------------------------"
    call f_sigma_vp(flux, sigma, lambda_g, lambda_etoile, h_etoile, u_etoile, params)
    print*, "La vitesse de l'onde de choc est de  ", sigma, "m/s"
    
    
    ! -- Calcul de l'erreur
    print*, "----------------------------------"
    errorL1 = Error_fct(flux%Unp1(:,1), sol_exa_h, dx, 1)

    errorL2 = Error_fct(flux%Unp1(:,1), sol_exa_h, dx, 2)

    errorLi = Error_fct(flux%Unp1(:,1), sol_exa_h, dx, 3)

    errorH1 = Error_fct(flux%Unp1(:,1), sol_exa_h, dx, 4)

    print*, "----------------------------------"

    ! -- Incrément du compteur
    cpt = cpt + 1
    Vect_err(cpt) = errorL2
    Vect_imax(cpt) = imax
    ! Désalocation des tableaux
    deallocate(x_i, h_i, u_i)
    !deallocate(flux%f_h, flux%f_q, flux%hnp1, flux%unp1)
    deallocate(sol_exa_h, sol_exa_u)
    deallocate(flux%Unp1)
    ! Fermeture des fichiers
    close(10)
    close(11)
    close(12)

    ! ######################## ORDRE #################### !
    if (flux%choix_approx_flux == 1 .AND. cpt < 2) then
        imax = imax *2
        goto 10
    end if

    print*, "#############################################"
    print*, "Schéma d'ordre = ", log(Vect_err(1)/Vect_err(2))/log(2.)
    print*, "#############################################"

    !call OUT(Vect_imax, Vect_err)




    ! -- Expérience -- !
    !!-- Cette partie va permettre de visualiser l'évolution -- !!
    !! -- de la hauteur d'eau sur un capteur -- !
    print*,"Calcul évolution des capteurs ..."

    !call evolution_capteur(name_file,imax,iter,dt,nb_capteurs,Pos_capteurs)

    print*,"Fin."





end program