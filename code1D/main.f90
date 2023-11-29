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
    real(pr), dimension(:), allocatable :: sol_exa_h, sol_exa_u
    real(pr)                            :: errorL1, errorL2, errorLi
    real(pr), dimension(6)             :: Vect_err 
    integer, dimension(6)              :: Vect_imax
    integer :: i, iter, Nmax, imax, cpt

    ! -- Sécurité
    Nmax = 10**7
    ! -- Initialisation d'un compteur
    cpt = 0

    open(unit = 2, file = "params.dat", action = "read")
        read(2,*) imax
    close(2)

    print*, "Entrer votre choix d'approximation du flux : "
    print*, " 1) Approximation de Rusanov d'ordre 1"
    print*, " 2) Approximation de Van Leer d'ordre 2"
    read*, flux%choix_approx_flux

    10 continue

    ! -- Ouverture du fichier d'écriture des résultats d'approximation
    if (flux%choix_approx_flux == 1) then
        open(unit = 10, file = "OUT/approx_RUSANOV.dat", action = "write")
    else if (flux%choix_approx_flux == 2) then
        open(unit = 10, file = "OUT/approx_Van_Leer.dat", action = "write")
    else
        print*,"Le choix de flux que vous avez fait n'est pas valide, veuillez recommencer"
        stop
    end if

    ! -- Allocation des tableaux de flux de hauteur et de vitesse
    allocate(x_i(0:imax+1))
    allocate(h_i(0:imax+1), u_i(0:imax+1))
    allocate(flux%f_h(0:imax+1), flux%f_q(0:imax+1))
    allocate(flux%hnp1(0:imax+1), flux%unp1(0:imax+1))
    allocate(sol_exa_h(0:imax+1), sol_exa_u(0:imax+1))

    ! -- Initialisation maillage
    call maillage(imax, tmax, dx, cfl, x_i)

    ! initialisation maillage + condition initiale
    call initialisation(imax, x_i, h_i, u_i)

    ! -- Initialisation des flux 
    flux%hnp1 = h_i
    flux%unp1 = u_i

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

    print*, "Passage n°", cpt+1
    print*, "-----------------------------------------"
    print*, "dx =", dx
    print*, "tmax =", tmax

    ! -- Boucle en temps
    do while (tn <= tmax .AND. iter < Nmax)

        iter = iter + 1
        ! -- Solution exacte
        call sol_exact_tn(flux, tn, x_i, sol_exa_h, sol_exa_u)

        ! -- Solution approchée
        call sol_approx_tn(flux, dt, cfl, dx, tn)

        ! -- Mise à jour des flux
        h_i = flux%hnp1
        u_i = flux%unp1

        ! -- Écriture dans les fichiers un fichier .dat

        do i = 0, imax+1
            write(10,*) x_i(i), h_i(i)
        end do

        write(10,*)
        write(10,*)

        do i = 0, imax+1
            write(11,*) x_i(i), sol_exa_h(i)
        end do

        write(11,*)
        write(11,*)

    end do

    ! -- Calcule de l'erreur
    print*, "-------------------------------------------"
    errorL1 = Error_fct(h_i, sol_exa_h, 1, dx)
    errorL2 = Error_fct(h_i, sol_exa_h, 2, dx)
    errorLi = Error_fct(h_i, sol_exa_h, 3, dx)
    print*, "-------------------------------------------"

    ! -- Incrément du compteur
    cpt = cpt + 1
    Vect_err(cpt) = errorL1
    Vect_imax(cpt) = imax
    ! Désalocation des tableaux
    deallocate(x_i, h_i, u_i)
    deallocate(flux%f_h, flux%f_q, flux%hnp1, flux%unp1)
    deallocate(sol_exa_h, sol_exa_u)
    ! Fermeture des fichiers
    close(10)
    close(11)

    ! ######################## ORDRE #################### !
    if (flux%choix_approx_flux == 1 .AND. cpt < 2) then
        imax = imax *2
        goto 10
    end if

    !call OUT(Vect_imax, Vect_err)
    print*, "#############################################"
    print*, "Schéma d'ordre = ", log(Vect_err(1)/Vect_err(2))/log(2.)
    print*, "#############################################"


end program