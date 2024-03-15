program main

    use mod_constantes
    use mod_initialize
    use mod_flux
    use mod_ExactSol
    use mod_TimeSchemes
    use mod_ApproxSol
    use mod_error

    implicit none

    real(pr)                              :: tn, tmax, dt, x0, y0, dx
    real(pr), dimension(:,:), allocatable :: ApproxSol, ExactSol
    real(pr), dimension(:), allocatable   :: bi
    real(pr)                              :: errdt, errdt2

    integer :: imax

    !imax = 1
    imax = 1

    ! -- Initialisation
    allocate(ApproxSol(0:imax+1,1),ExactSol(1:imax,1))
    allocate(bi(1:imax+1))
    x0 = -2
    y0 = -2
    ApproxSol(0,1) = x0
    ApproxSol(1,1) = y0
    bi = 0.
    tn = 0.



    ! -- Ouverture des fichier d'écriture des résultats            print*, "passe"
    open(unit = 11, file = "out/test_RK2.dat", action = "write")

    write(11,*) tn, ApproxSol(0,1), x0
    !write(11,*) tn, ApproxSol(1,1), ApproxSol(2,1), x0, y0

    ! -- Boucle en temps
    print*, "Time loop ..."
    dt = 0.1
    tn = 0.1
    dx = 0.1
    tmax = 10

    do while (tn <= tmax )

        ApproxSol(1:imax,:) = RK2SSP(tn, ApproxSol, dt, imax, bi, dx)
        ExactSol(:,1) = exa(x0, y0, 2, tn, imax)

        !write(11,*) tn, ApproxSol(0, 1), ExactSol(1,1)
        write(11,*) tn, ApproxSol(0, 1), ApproxSol(1,1), ExactSol(1,1), ExactSol(2,1)

        tn = tn + dt

    end do

    print*, "Time loop completed. ----------------"

    errdt = abs(ApproxSol(0,1) - ExactSol(1,1)) !+ abs(ApproxSol(2,1) - ExactSol(2,1))

    print*, "Seconde Time loop ..."

    dt = dt/2.
    tn = 0.05
    ApproxSol(0,1) = x0
    ApproxSol(1,1) = y0
    bi = 0.

    do while (tn <= tmax )

        ApproxSol = RK2SSP(tn, ApproxSol, dt, imax, bi, dx)
        ExactSol(:,1) = exa(x0, y0, 2, tn, imax)

        ! write(11,*) tn, ApproxSol(1, 1), ExactSol(1,1)
        ! write(11,*) tn, ApproxSol(1, 1), ApproxSol(2,1), ExactSol(1,1), ExactSol(2,1)

        tn = tn + dt

    end do

    errdt2 = abs(ApproxSol(0,1) - ExactSol(1,1)) !+ abs(ApproxSol(2,1) - ExactSol(2,1))

    close(11)

    print*, log(errdt/errdt2)/log(2.)

    print*, " "
    print*, "Done. --------------"

contains

    function exa(x,y,ch,t,dim)result(fct)

        integer, intent(in)  :: ch, dim
        real(pr), intent(in) :: x, y, t

        real(pr), dimension(dim) :: fct

        ! -- Local
        

        if (ch==1) then
            fct(1) = x*cos(t) - y*sin(t)
            fct(2) = y*cos(t) + x*sin(t)
        else if (ch==2) then
            fct(1) = 2*x/(2-t**2*x)
        end if

    end function



end program