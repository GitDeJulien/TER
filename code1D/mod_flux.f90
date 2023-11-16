module mod_flux

    use mod_constantes
    use mod_vague_rupture

    implicit none

    type flux_type
        !choix du type d'approximation de flux
        integer      :: choix_approx_flux
        !les vecteurs d'entrées h et u (hauteur et vitesse init)
        real(pr), dimension(:), allocatable :: hnp1, unp1
        real(pr), dimension(:), allocatable :: f_h, f_q

        integer      :: n
    end type

contains


    subroutine sol_approx_tn(flux, dt, cfl, dx)

        ! - Externe
        type(flux_type), intent(inout)       :: flux
        real(pr), intent(in)                 :: dx, cfl
        real(pr), intent(inout)              :: dt
        

        ! - Interne
        integer                  :: i, imax
        real(pr)                 :: h, q
        real(pr), dimension(2)   :: Ugauche, Udroite
        real(pr),dimension(2)    :: F
        real(pr)                 :: unsurdx, det_A
        real(pr), dimension(2,2) :: A
        real(pr), dimension(size(flux%hnp1)) :: b_i


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

        select case(flux%choix_approx_flux)

        case(1)

        do i = 0, imax

            Ugauche(1) = flux%hnp1(i)
            Ugauche(2) = flux%hnp1(i)*flux%unp1(i)
            Udroite(1) = flux%hnp1(i+1)
            Udroite(2) = flux%hnp1(i+1)*flux%unp1(i+1)

            F = flux_num(flux%choix_approx_flux, Ugauche, Udroite, b_i(i+1), unsurdx, dt, det_A)

            flux%f_h(i+1) = F(1)
            flux%f_q(i+1) = F(2)

        end do

        case(2)

            do i = 0, imax

                call A_roe(flux%hnp1(i-1), flux%hnp1(i), flux%hnp1(i+1), flux%unp1(i-1), flux%unp1(i), &
                & flux%unp1(i+1), A, det_A, unsurdx)

                Ugauche(1) = flux%hnp1(i)
                Ugauche(2) = flux%hnp1(i)*flux%unp1(i)
                Udroite(1) = flux%hnp1(i+1)
                Udroite(2) = flux%hnp1(i+1)*flux%unp1(i+1)
    
                F = flux_num(flux%choix_approx_flux, Ugauche, Udroite, b_i(i+1), unsurdx, dt, det_A)
    
                flux%f_h(i+1) = F(1)
                flux%f_q(i+1) = F(2)
    
            end do

        end select


        ! -- Mise à jour des h et q au nouveau pas de temps sans les bords
        do i = 1, imax

            h = flux%hnp1(i) - dt*unsurdx*(flux%f_h(i+1) - flux%f_h(i))
            q = flux%hnp1(i)*flux%unp1(i) - dt*unsurdx*(flux%f_q(i+1) - flux%f_q(i))

            flux%hnp1(i) = h
            flux%unp1(i) = q/h

        end do


        ! condition de Neumann sur les bords
        flux%hnp1(0) = flux%hnp1(1)
        flux%unp1(0) = flux%unp1(1)
        flux%hnp1(imax+1) = flux%hnp1(imax)
        flux%unp1(imax+1) = flux%unp1(imax)


    end subroutine sol_approx_tn

! ###########################################################################################

    ! -- Calcul du flux numérique
    function flux_num(choix_flux, Ug, Ud, bi, unsurdx, dt, det_A)result(F)

        real(pr), dimension(2), intent(in)  :: Ug, Ud
        integer, intent(in)                 :: choix_flux
        real(pr), intent(in)                :: bi
        real(pr), intent(in)                :: unsurdx, dt, det_A
        real(pr), dimension(2)              :: F 

        select case(choix_flux)

        case(1)

            F(1) = 0.5*(Ud(2) + Ug(2)) - 0.5*bi*(Ud(1) - Ug(1))
            F(2) = (Ug(2)*Ug(2)/Ug(1) + g*Ug(1)*Ug(1)/2. + Ud(2)*Ud(2)/Ud(1) + g*Ud(1)*Ud(1)/2.)*0.5 - 0.5*bi*(Ud(2) - Ug(2))

        case(2)

            F(1) = 0.5*(Ud(2) + Ug(2)) - 0.5*det_A*(Ud(1) - Ug(1))
            F(2) = (Ug(2)*Ug(2)/Ug(1) + g*Ug(1)*Ug(1)/2. + Ud(2)*Ud(2)/Ud(1) + g*Ud(1)*Ud(1)/2.)*0.5 - 0.5*det_A*(Ud(2) - Ug(2))

        end select


    end function flux_num

! #############################################################################################

    ! -- Clacul du coefficient bi
    function f_bi(hg, hd, ug, ud)result(bi)
        real(pr), intent(in)                :: hg, hd
        real(pr), intent(in)                :: ud, ug
        real(pr)                            :: bi

        bi = max(abs(ug + sqrt(g*hg)), abs(ud + sqrt(g*hd)), abs(ug - sqrt(g*hg)), abs(ud - sqrt(g*hd)))

    end function f_bi

    ! -- Calcul de la matrice de Roe
    subroutine A_roe(hg, hm, hd, ug, um, ud, A, det_A, unsurdx)
        real(pr), intent(in)                  :: hg, hm, hd, ug, um, ud, unsurdx
        real(pr), dimension(:,:), intent(out) :: A 
        real(pr), intent(out)                 :: det_A

        A(1,1) = 1/2.*(ud-ug)+unsurdx*um
        A(1,2) = 1/2.*(hd-hg)+unsurdx*hm
        A(2,1) = unsurdx*g
        A(2,2) = 1/2.*(ud-ug)+unsurdx*um

        det_A = (1/2.*(ud-ug)+unsurdx*um)**2 - unsurdx*g*1/2.*(hd-hg)+unsurdx*hm

    end subroutine

! #############################################################################################

    function Error_fct(approx, exa, dx, norme)result(err)

        ! -- Code pour les erreur 
        ! - norme = 1 (L1)
        ! - norme = 2 (L2)
        ! - norme = 3 (Linfini)

        real(pr), dimension(:), intent(in) :: approx, exa
        real(pr), intent(in)               :: dx
        integer, intent(in)                :: norme

        real(pr)                           :: err

        ! -- Locale
        integer :: i

        err = 0.
        select case(norme)
        case(1) ! Norme L1
            do i=1,size(approx)
                err = err + dx*abs(exa(i)-approx(i))
            end do

        case(2)
            do i=1,size(approx)
                err = err + dx*abs(exa(i)-approx(i))**2
            end do
            err = sqrt(err)

        case(3)
            err = maxval(abs(exa(:)-approx(:)))
        end select

    end function





end module