module mod_source
    use mod_constantes
    implicit none
    
contains

    function fond(x,forme) result(Z)             !Forme du fond
        ! Variables d'entrées  
        integer, intent(in) :: forme
        real(pr),intent(in) :: x

        ! Variable de sortie
        real(pr) :: Z

        ! Variables locale
        real(pr) :: alpha
        
        select case(forme)
        !------------Légere pente---------------!
        case(1)                 
            if(x<0._pr) then
                Z = 0._pr
            else
                alpha = 6._pr
                Z = -alpha*x       !En réalite Z = -tan(alpha)*x + hinit mais pas important pour la suite
            end if

        !------------Obstacle---------------!
            case(2)                 !Obstacle en x = ???
            if (x>12._pr .AND. x<25._PR) then
                Z = 25._pr
            else 
                Z = 0._pr
            end if

        !------------Fond plat (permet de vérifier si ça marche)---------------!
        case(3)                 !
            Z = 0._pr
        end select

    end function fond
    
    subroutine Build_source(h_i, x_i, F, forme, model)
        ! Variables d'entrée
        real(pr), dimension(0:), intent(inout)           :: h_i, x_i
        integer, intent(in)                             :: forme, model 

        ! Variable de sortie
        real(pr),dimension(:),allocatable,intent(out)   :: F

        ! Variables locales
        real(pr)                              :: hmoy, Zmoyp1, Zmoym1
        real(pr)                              :: xip1,xi,xim1,dx
        real(pr)                              :: hi,hip1
        integer                               :: i, imax

        imax = size(x_i) - 2 
        allocate(F(0:imax+1))

        select case(model)
        case(1)
            hmoy = h_i(0)
            xi = x_i(0)
            xip1 = x_i(1)
            Zmoym1 = fond(xi,forme)
            Zmoyp1 = (fond(xip1,forme) + fond(xi,forme))/2
            dx = x_i(1) - xip1
            F(0) = -g*hmoy*(Zmoyp1 - Zmoym1)/(2*dx)

            do i = 1, imax

                xi = x_i(i)
                xim1 = x_i(i-1)
                xip1 = x_i(i+1)

                hi = h_i(i)
                hip1 = h_i(i+1)

                hmoy = (hi + hip1)/2

                Zmoyp1 = (fond(xip1,forme) + fond(xi,forme))/2
                Zmoym1 = (fond(xi,forme) + fond(xim1,forme))/2

                dx = xip1 - xim1                                    ! 2*dx

                F(i) = -g*hmoy*(Zmoyp1 - Zmoym1)/dx
            end do
            hmoy = h_i(imax+1)
            xi = x_i(imax+1)
            xim1 = x_i(imax)
            Zmoym1 = (fond(xi,forme) + fond(xim1,forme))/2
            Zmoyp1 = fond(xi,forme)
            F(imax + 1) = -g*hmoy*(Zmoyp1 - Zmoym1)/dx
        end select
    end subroutine Build_source
end module mod_source