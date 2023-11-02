program main

    use mod_vague_rupture
    use mod_constantes
    use mod_flux

    implicit none


    
    type(flux_type)                     :: flux
    real(pr)                            :: tn, dx, tmax, dt, cfl
    real(pr), dimension(:), allocatable :: x_i, h_i, u_i

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

    ! -- Initialisation des flux 
    flux%hnp1 = h_i
    flux%unp1 = u_i

    print*, "Entrer votre choix d'approximation du flux : "
    print*, " 1) Approximation de Lax-Friedrichs"
    read*, flux%choix_approx_flux

    dt = 0.01

    print*, "-----------------------------------------"
    print*, "dx =", dx
    print*, 'dt =', dt
    print*, "-----------------------------------------"

    tn = 0.
    iter = 0

    do while (tn <= tmax .AND. iter < Nmax)

        iter = iter + 1
        tn = tn + dt

        call sol_approx_tn(flux, dt, cfl, dx)

        h_i = flux%hnp1
        u_i = flux%unp1

        do i = 0, imax+1
            write(10,*) x_i(i), h_i(i)
        end do

        write(10,*)
        write(10,*)


    end do


    deallocate(x_i, h_i, u_i)
    deallocate(flux%f_h, flux%f_q, flux%hnp1, flux%unp1)

    close(10)


end program