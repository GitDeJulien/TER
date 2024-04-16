program main

    use mod_constantes
    use mod_initialize
    use mod_flux
    use mod_ExactSol
    use mod_TimeSchemes
    use mod_ApproxSol
    use mod_error
    use mod_analyse_exp

    implicit none


    
    type(approx_type)                     :: app
    real(pr)                              :: tn, dx, tmax, dt, cfl, sigma
    real(pr), dimension(:), allocatable   :: x_i
    real(pr), dimension(:,:), allocatable :: ExactSol
    real(pr)                              :: errL2, errL1, errLi
    !real(pr)                            :: sigma, lambda_g, lambda_etoile, h_etoile, u_etoile
    character(len=40)                     :: params, nameOUT, num, date

    integer :: i, iter, imax

    ! params = "entry/parameters.txt"
    print*, "Entrer la date de l'expérience ? (JJ_MM)"
    read(*,*) date
    print*, "Entrer le numéro de l'expérience ? "
    read(*,*) num

    params = "entry/exp_"//trim(date)//"/param/exp_"//trim(num)//".dat"

    open(unit = 10, file = params, action = "read")
        read(10,*) imax
    close(10)

    print*, "Entrer votre choix d'approximation du flux : "
    print*, " 1) Approximation de Rusanov d'ordre 1"
    print*, " 2) Approximation de Rusanov + Mood d'ordre 2"
    read*, app%choix_approx_flux

    ! -- Ouverture du fichier d'écriture des résultats d'approximation
    if (app%choix_approx_flux == 1) then
        nameOUT = "out/RUSANOV_"//trim(date)//"_exp_"//trim(num)//".dat"
        open(unit = 10, file = nameOUT, action = "write")
    else if (app%choix_approx_flux == 2) then
        nameOUT = "out/MOOD_"//trim(date)//"_exp_"//trim(num)//".dat"
        open(unit = 10, file = nameOUT, action = "write")
    else
        print*,"The choice of approximation you made is not correct, please chose 1 or 2"
        stop
    end if

    print*, "Initialization ..."
    
    allocate(x_i(0:imax+1))
    allocate(ExactSol(0:imax+1, 2))
    allocate(app%Unp1(0:imax+1, 2))
    
    call initialize(imax, tmax, dx, x_i, app%Unp1, cfl, params)

    print*, "imax = ", imax
    print*, "tmax = ", tmax
    print*, "dx = " , dx
    print*, "cfl = ", cfl
    
    print*, "Initialisation completed. ------------------"

    ! -- Ouverture des fichier d'écriture de la solution exacte
    open(unit = 11, file = "out/Exact.dat", action = "write")

    ! -- Boucle en temps
    print*, "Time loop ..."
    tn = 0.
    iter = 0

    do while (tn <= tmax )

        iter = iter + 1

        ! -- Solution approchée
        call SolApprox(app, dt, cfl, dx, tn, imax, iter)

        ! -- Solution exacte
        call ExactSolFct(app, tn, x_i, ExactSol(:,1), ExactSol(:,2), sigma, params)

        do i = 0, imax+1
            write(10,*) x_i(i), app%Unp1(i, 1)
        end do
        write(10,*)
        write(10,*)

        do i = 0, imax+1
            write(11,*) x_i(i), ExactSol(i, 1)
        end do

        write(11,*)
        write(11,*)


    end do

    print*, "Time loop completed. ----------------"
    print*, " "
    print*, "Vitesse de l'onde de choc :", sigma, "m/s"

    print*, "Error computation ..."
    errL1 = Error_fct(app%Unp1(:,1), ExactSol(:,1), dx, 1)
    errL2 = Error_fct(app%Unp1(:,1), ExactSol(:,1), dx, 2)
    errLi = Error_fct(app%Unp1(:,1), ExactSol(:,1), dx, 3)
    print*, "---------------------------------"

    close(10)
    close(11)

    ! -- Expérience -- !
    !!-- Cette partie va permettre de visualiser l'évolution -- !!
    !! -- de la hauteur d'eau sur un capteur -- !
    print*,"Calcul évolution des capteurs ..."
    print*, dt
    dt=tmax/iter
    print*, dt
    call evolution_capteur(params, nameOUT, imax,iter,dt, date, num)
    print*, " "
    print*, "Done. --------------"



end program