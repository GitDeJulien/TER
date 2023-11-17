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
        real(pr), dimension(2)   :: Ugauche, Udroite, Umid
        real(pr),dimension(2)    :: F
        real(pr)                 :: unsurdx, det_J, det_J2
        real(pr), dimension(2,2) :: J, J2
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

            F = flux_num(flux%choix_approx_flux, Ugauche, Udroite, Umid, b_i(i+1), unsurdx, dt, det_J)

            flux%f_h(i+1) = F(1)
            flux%f_q(i+1) = F(2)

        end do

        ! case(2)

        !     do i = 1, imax

        !         call A_roe(flux%hnp1(i-1), flux%hnp1(i), flux%hnp1(i+1), flux%unp1(i-1), flux%unp1(i), flux%unp1(i+1)&
        !         , J, det_J, unsurdx)

        !         Ugauche(1) = flux%hnp1(i)
        !         Ugauche(2) = flux%hnp1(i)*flux%unp1(i)
        !         Udroite(1) = flux%hnp1(i+1)
        !         Udroite(2) = flux%hnp1(i+1)*flux%unp1(i+1)
    
        !         F = flux_num(flux%choix_approx_flux, Ugauche, Udroite, Umid, b_i(i+1), unsurdx, dt, det_J)
    
        !         flux%f_h(i+1) = F(1)
        !         flux%f_q(i+1) = F(2)
    
        !     end do

        ! case(3)

        !     do i = 1, imax

        !         call Jacobienne(flux%hnp1(i+1), flux%unp1(i+1), J, det_J, J2, det_J2)
        !         Umid(1) = flux%hnp1(i)
        !         Umid(2) = flux%hnp1(i)*flux%unp1(i)
        !         Ugauche(1) = flux%hnp1(i-1)
        !         Ugauche(2) = flux%hnp1(i-1)*flux%unp1(i-1)
        !         Udroite(1) = flux%hnp1(i+1)
        !         Udroite(2) = flux%hnp1(i+1)*flux%unp1(i+1)
    
        !         F = flux_num(flux%choix_approx_flux, Ugauche, Udroite, Umid, b_i(i+1), unsurdx, dt, det_J2)
    
        !         flux%f_h(i+1) = F(1)
        !         flux%f_q(i+1) = F(2)
    
        !     end do

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

        ! flux%unp1(1) = flux%unp1(2)
        ! flux%hnp1(imax) = flux%hnp1(imax-1)
        ! flux%unp1(imax) = flux%unp1(imax-1)


    end subroutine sol_approx_tn

! ###########################################################################################

    ! -- Calcul du flux numérique
    function flux_num(choix_flux, Ug, Ud, Um, bi, unsurdx, dt, det_J)result(F)

        real(pr), dimension(2), intent(in)  :: Ug, Ud, Um
        integer, intent(in)                 :: choix_flux
        real(pr), intent(in)                :: bi
        real(pr), intent(in)                :: unsurdx, dt, det_J
        real(pr), dimension(2)              :: F 

        select case(choix_flux)

        case(1)

            F(1) = 0.5*(Ud(2) + Ug(2)) - 0.5*bi*(Ud(1) - Ug(1))
            F(2) = (Ug(2)*Ug(2)/Ug(1) + g*Ug(1)*Ug(1)/2. + Ud(2)*Ud(2)/Ud(1) + g*Ud(1)*Ud(1)/2.)*0.5 - 0.5*bi*(Ud(2) - Ug(2))

        case(2)

            F(1) = 0.5*(Ud(2) + Ug(2)) - 0.5*det_J*(Ud(1) - Ug(1))
            F(2) = (Ug(2)*Ug(2)/Ug(1) + g*Ug(1)*Ug(1)/2. + Ud(2)*Ud(2)/Ud(1) + g*Ud(1)*Ud(1)/2.)*0.5 - 0.5*det_J*(Ud(2) - Ug(2))

        case(3)

            F(1) = 0.5*(Um(2) + Ud(2))- dt*unsurdx*0.5*(Ud(2)/Ud(1)*(Ud(1)-Um(1))+Ud(1)*(Ud(2)-Um(2)) - &
            dt*unsurdx*0.5*((Ud(2)*Ud(2)/(Ud(1)*Ud(1)) + g*Ud(1))*(Ud(1)-2*Um(1)+Ug(1)) + 2*Ud(2)*(Ud(2)-2*Um(2)+Ug(2))))

            F(2) = 0.5*(Um(2)*Um(2)/Um(1)+g*Um(1)**2/2. + Ud(2)*Ud(2)/Ud(1)+g*Ud(1)**2/2.)- &
            dt*unsurdx*0.5*(g*(Ud(1)-Um(1)+Ud(2)/Ud(1)*(Ud(2)-Um(2))) - dt*unsurdx*0.5*(2*g*Ud(2)/Ud(1)*(Ud(1)-2*Um(1)+Ug(1))&
            +(Ud(1)*g + Ud(2)*Ud(2)/(Ud(1)*Ud(1))) * (Ud(2)-2*Um(2)+Ug(2))))

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


    ! -- Calcul de la matrice de Roe
    subroutine Jacobienne(h, u, J, det_J, J2, det_J2 )
        real(pr), intent(in)                  :: h, u
        real(pr), dimension(:,:), intent(out) :: J, J2 
        real(pr), intent(out)                 :: det_J, det_J2

        J(1,1) = u
        J(1,2) = h
        J(2,1) = g
        J(2,2) = u

        J2(1,1) = u**2 + g*h
        J2(1,2) = 2*u*h
        J2(2,1) = 2*u*g
        J2(2,2) = h*g + u**2

        det_J = u**2 - g*h
        det_J2 = (u**2 + g*h) * (h*g + u**2) - 4*u*u*g*h

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