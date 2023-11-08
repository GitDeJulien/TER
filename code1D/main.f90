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

    integer :: i, iter, Nmax, imax

    ! -- Sécurité
    Nmax = 10**6

    ! initialisation maillage + condition initiale
    call initialisation(imax, tmax, dx, x_i, h_i, u_i, cfl)

    ! -- Ouverture du fichier d'écriture des résultats
    open(unit = 10, file = "OUT/premiere_sol.dat", action = "write")
    ! -- Écriture de la condition initiale
    do i = 0, imax+1
        write(10,*) x_i(i), h_i(i)
    end do
    write(10,*)
    write(10,*)


    ! -- Allocation des tableaux de flux de hauteur et de vitesse
    allocate(flux%f_h(0:imax+1), flux%f_q(0:imax+1))
    allocate(flux%hnp1(0:imax+1), flux%unp1(0:imax+1))
    allocate(sol_exa_h(0:imax+1), sol_exa_u(0:imax+1))

    ! -- Initialisation des flux 
    flux%hnp1 = h_i
    flux%unp1 = u_i

    print*, "Entrer votre choix d'approximation du flux : "
    print*, " 1) Approximation de Lax-Friedrichs"
    read*, flux%choix_approx_flux

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

    do while (tn <= tmax .AND. iter < Nmax)

        iter = iter + 1
        tn = tn + dt

        call sol_approx_tn(flux, dt, cfl, dx)
        call sol_exact_tn(flux, tn, x_i, sol_exa_h, sol_exa_u)

        h_i = flux%hnp1
        u_i = flux%unp1

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


    deallocate(x_i, h_i, u_i)
    deallocate(flux%f_h, flux%f_q, flux%hnp1, flux%unp1)
    deallocate(sol_exa_h, sol_exa_u)

    close(10)
    close(11)


end program