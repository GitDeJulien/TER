program main

    use mod_constantes
    use mod_flux

    implicit none

    ! initialisation variables

    ! Demander à l'utilisateur de choisir sont sa condition init 
    ! et sont flux

    ! initialisation maillage + vague (init_vague)

    type(vague_type)                    :: vague
    type(flux_type)                     :: flux
    real(pr)                            :: t, L, dx, tmax, dt
    real(pr)                            :: cfl
    real(pr), dimension(2)              :: F
    real(pr)                            :: h, q
    real(pr)                            :: unsurdx

    integer :: i, iter, ntmax



    print*, "Entrer un temps de simulation en seconde :"
    read*, tmax
    print*, "------------------------------------------"

    call init_vague_rupture(vague, L, dx)


    t = 0.



    !Ouverture du fichier d'écriture des résultats
    open(unit = 1, file = "OUT/premiere_sol.dat", action = "write")
    do i = 0, vague%n+1
        write(1,*) vague%x_i(i), vague%hi_0(i)!, vague%qi_0(i)
    end do
    write(1,*)
    write(1,*)


    allocate(flux%f_h(0:vague%n+2), flux%f_q(0:vague%n+2))
    allocate(flux%hnp1(0:vague%n+2), flux%qnp1(0:vague%n+2))

    flux%hnp1 = vague%hi_0
    flux%qnp1 = vague%qi_0

    print*, "hnp1_0 = ", flux%hnp1(0)
    print*, "------------------------------------------"

    print*, "Entrer  votre choix d'approximation du flux : "
    print*, "-> 1) Moyenne aux bords de mailles"
    read*, flux%choix_approx_flux


    cfl = 10._pr**(-2)
    dt = cfl * dx
    unsurdx = 1/dx
    print*, "dx =", dx
    print*, "dt =", dt

    ntmax = int(tmax/dt)

    do iter = 1, ntmax

        t = iter * dt


        do i = 1, vague%n+1

            F = flux_num(flux%choix_approx_flux, flux%hnp1(i-1), flux%hnp1(i+1), flux%qnp1(i-1), flux%qnp1(i+1))

            flux%f_h(i) = F(1)
            flux%f_q(i) = F(2)

        end do


        do i = 1, vague%n+1


            h = flux%hnp1(i) - dt*unsurdx*(flux%f_h(i))
            !print*, "h =", h

            if (flux%hnp1(i)>0) then

                q = flux%hnp1(i)*flux%qnp1(i) - dt*unsurdx*(flux%f_q(i))

                flux%qnp1(i) = q/h
                flux%hnp1(i) = h
            else
                flux%qnp1(i) = 0.
                flux%hnp1(i) = h
            end if

        end do

        ! condition de Neumann sur les bords
        flux%hnp1(0) = flux%hnp1(1)
        flux%qnp1(0) = flux%qnp1(1)
        flux%hnp1(flux%n+1) = flux%hnp1(flux%n)
        flux%qnp1(flux%n+1) = flux%qnp1(flux%n)

        do i = 0, vague%n+1
            write(1,*) vague%x_i(i), flux%hnp1(i)!, flux%qnp1(i)
        end do

        write(1,*)
        write(1,*)


    end do



    deallocate(vague%hi_0, vague%qi_0, vague%x_i)
    deallocate(flux%f_h, flux%f_q, flux%hnp1, flux%qnp1)

    close(1)


end program