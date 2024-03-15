module mod_ApproxSol

    use mod_constantes
    use mod_flux
    use mod_TimeSchemes

    implicit none

    type approx_type

    integer      :: choix_approx_flux
    real(pr), dimension(:,:), allocatable :: Unp1

    end type

contains

    subroutine SolApprox(app, dt, cfl, dx, tn, imax, iter)

        ! - Externe
        type(approx_type), intent(inout)     :: app
        integer, intent(inout)               :: iter
        real(pr), intent(in)                 :: dx, cfl
        real(pr), intent(inout)              :: dt, tn
        integer, intent(in)                  :: imax
        

        ! - Interne
        integer                          :: i, test
        real(pr), dimension(1:imax+1, 2) :: Ud ! - Taille imax+1 (1:imax+1)
        real(pr), dimension(0:imax, 2)   :: Ug ! - Taille imax+1  (0:imax)
        real(pr), dimension(1:imax+1, 2) :: F  ! - Taille imax
        real(pr), dimension(0:imax+1, 2) :: Un, Un_temp ! - Taille imax+2 (0:imax+1)
        real(pr)                         :: unsurdx, eps
        real(pr), dimension(1:imax+1)    :: b_i

        unsurdx = 1./dx

        ! -- Tableau des bi
        do i=0, imax
            b_i(i+1) = f_bi(app%Unp1(i,1), app%Unp1(i+1,1), app%Unp1(i,2)/app%Unp1(i,1), app%Unp1(i+1,2)/app%Unp1(i+1,1))
        end do

        ! -- Calcule du pas de temps respectant la cfl
        dt = cfl * dx/(2 * MAXVAL(b_i))
        tn = tn + dt

        Un(0:imax+1,:) = app%Unp1(0:imax+1,:)

        select case(app%choix_approx_flux)

        case(1)

            Ug(0:imax,:) = Un(0:imax,:)
            Ud(1:imax+1,:) = Un(1:imax+1,:)

            F = FluxNum(imax, Ug, Ud, b_i)

            ! -- Schéma volumes finies
            app%Unp1(1:imax,:) = Un(1:imax,:) - dt*unsurdx*(F(2:imax+1,:) - F(1:imax,:))

            ! -- Condition de bord - Neumann
            app%Unp1(0,:) = app%Unp1(1,:)
            app%Unp1(imax+1,:) = app%Unp1(imax,:)

        case(2)

            ! app%Unp1(1:imax, :) = RK2SSP(tn, Un(0:imax+1, :), dt, imax, b_i, unsurdx)
            
            Un_temp(1:imax,:) = RK2SSP(tn, Un(0:imax+1,:), dt, imax, b_i, dx)
            
            ! -- Condition de bord - Neumann
            Un_temp(0,:) = app%Unp1(1,:)
            Un_temp(imax+1,:) = app%Unp1(imax,:)

            test = 0
            do i=0,imax+1
                if (Un_temp(i,1) < 0.0) then
                    test = 1
                end if
            end do
            eps = 0.00132
            do i=0,imax
                if (min(Un(i,1),Un(i+1,1)) - eps> Un_temp(i,1) .OR. max(Un(i,1),Un(i+1,1)) + eps < Un_temp(i,1)) then
                    test = 1
                end if
            end do

            if (test == 1) then
                Ug(0:imax,:) = Un(0:imax,:)
                Ud(1:imax+1,:) = Un(1:imax+1,:)
    
                F = FluxNum(imax, Ug, Ud, b_i)
    
                ! -- Schéma volumes finies
                app%Unp1(1:imax,:) = Un(1:imax,:) - dt*unsurdx*(F(2:imax+1,:) - F(1:imax,:))
    
                ! -- Condition de bord - Neumann
                app%Unp1(0,:) = app%Unp1(1,:)
                app%Unp1(imax+1,:) = app%Unp1(imax,:)

                print*, "Ordre 1 :", iter

            else
                app%Unp1 = Un_temp
            end if

        end select


    end subroutine SolApprox

end module
    