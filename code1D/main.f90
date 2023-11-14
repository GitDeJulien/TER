program main

    use mod_vague_rupture
    use mod_constantes
    use mod_flux
    use mod_sol_exact_2

    implicit none


    
    type(flux_type)                     :: flux
    real(pr)                            :: tn, dx, tmax, dt, cfl
    real(pr), dimension(:), allocatable :: x_i, h_i, u_i
    real(pr), dimension(:), allocatable :: sol_exa_h, sol_exa_u
    real(pr)                            :: error

    integer :: i, iter, Nmax, imax

    ! -- Sécurité
    Nmax = 10**6

    ! initialisation maillage + condition initiale
    call initialisation(imax, tmax, dx, x_i, h_i, u_i, cfl)

    ! -- Allocation des tableaux de flux de hauteur et de vitesse
    allocate(flux%f_h(0:imax+1), flux%f_q(0:imax+1))
    allocate(flux%hnp1(0:imax+1), flux%unp1(0:imax+1))
    allocate(sol_exa_h(0:imax+1), sol_exa_u(0:imax+1))

    ! -- Initialisation des flux 
    flux%hnp1 = h_i
    flux%unp1 = u_i

    print*, "Entrer votre choix d'approximation du flux : "
    print*, " 1) Approximation de Roe d'ordre 1"
    print*, " 2) Approximation de Roe d'ordre 2"
    print*, " 3) Approximation de Lax Vendroff d'ordre 2"
    read*, flux%choix_approx_flux

    ! -- Ouverture du fichier d'écriture des résultats d'approximation
    if (flux%choix_approx_flux == 1) then
        open(unit = 10, file = "OUT/approx_ROE_O(1).dat", action = "write")
    else if (flux%choix_approx_flux == 2) then
        open(unit = 10, file = "OUT/approx_ROE_O(2).dat", action = "write")
    else if (flux%choix_approx_flux == 3) then
        open(unit = 10, file = "OUT/approx_LW_O(2).dat", action = "write")
    else
        print*,"Le choix de flux que vous avez fait n'est pas valide, veuillez recommencer"
        stop
    end if
    
    ! -- Écriture de la condition initiale
    do i = 0, imax+1
        write(10,*) x_i(i), h_i(i)
    end do
    write(10,*)
    write(10,*)

    dt = 0.1
    tn = 0.
    iter = 0
    sol_exa_h = 0.
    sol_exa_u = 0.

    ! -- Calcule de la solution exacte au temps t = 0
    call sol_exact_tn(flux, tn, x_i, sol_exa_h, sol_exa_u)

    ! -- Ouverture du fichier d'écriture des résultats
    open(unit = 11, file = "OUT/sol_exact.dat", action = "write")


    do i = 0, imax+1
        write(11,*) x_i(i), h_i(i)
    end do
    write(11,*)
    write(11,*)

    print*, "-----------------------------------------"
    print*, "dx =", dx
    print*, "-----------------------------------------"

    tn = dt
    do while (tn <= tmax .AND. iter < Nmax)

        iter = iter + 1

        call sol_approx_tn(flux, dt, cfl, dx)
        h_i = flux%hnp1
        u_i = flux%unp1

        tn = tn + dt
        call sol_exact_tn(flux, tn, x_i, sol_exa_h, sol_exa_u)

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

        ! -- Mise à jour du pas de temps
        dt = 0.05

    end do

    ! -- Calcule de l'erreur
    print*, "----------------------------------"
    error = Error_fct(h_i, sol_exa_h, dx, 1)
    print*, "Error L1 = ", error
    error = Error_fct(h_i, sol_exa_h, dx, 2)
    print*, "Error L2 = ", error
    error = Error_fct(h_i, sol_exa_h, dx, 3)
    print*, "Error L_inifnit = ", error
    print*, "----------------------------------"


    deallocate(x_i, h_i, u_i)
    deallocate(flux%f_h, flux%f_q, flux%hnp1, flux%unp1)
    deallocate(sol_exa_h, sol_exa_u)

    close(10)
    close(11)


end program