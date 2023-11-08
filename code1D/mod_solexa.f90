module mod_solexa
    use mod_constantes
    use mod_flux
    use mod_vague_rupture

    implicit none


    
contains
!!...........!!
! Calcul des deux valeurs propres !
! Cela va nous permettre de déterminer certaines zones pour notre sol exacte !
!!...........!!


    
subroutine val_propres(hexa_i,uexa_i,lambda1,lambda2,imax)
    !Variables d'entrée
    real(pr),dimension(:),intent(in) :: hexa_i, uexa_i
    integer                          :: imax
    !Variables de sortie
    real(pr), intent(out)            ::  lambda1,lambda2
    !Variables locales
    real(pr)                         ::  B0,Bn  

    B0 = hexa_i(0)*uexa_i(0)
    Bn = hexa_i(imax+1)*uexa_i(imax+1)
    !Calcul des deux valeurs propres
    lambda1 = B0/hexa_i(0) - sqrt(g*hexa_i(0))
    lambda2 = Bn/hexa_i(imax+1) - sqrt(g*hexa_i(imax+1))

end subroutine val_propres
!!...........!!
! Solution à l'instant tn !
! et pour xi avec i de 0 à imax + 1!
!!...........!!
subroutine sol_exa_tn(hexa_i,uexa_i,t,x_i,imax)          !!sol_approx_tn(flux, dt, cfl, dx)
    !Variables d'entrée
    integer,intent(in)                   :: imax
    real(pr),intent(in)                  :: t       
    real(pr),dimension(:), intent(in)    ::  x_i
    !Variables de sortie
    real(pr),dimension(:), intent(inout) :: hexa_i,uexa_i
    
    !Variables interne
    real(pr)                             :: lambda1,lambda2,L,A
    integer                              :: i

    call val_propres(hexa_i,uexa_i,lambda1,lambda2,imax)                             

    L = x_i(imax+1) - x_i(0)

    do i = 0, imax + 1
        !!...........!!
        ! xi dans la détente !
        !!...........!!
        if (x_i(i) >= lambda1*t .AND. x_i(i) <= lambda2*t) then
            A = sqrt(g*hexa_i(0))
            hexa_i(i) = (1._pr/(9*g))*(2*A - x_i(i))
            uexa_i(i) = 2*(A - sqrt(hexa_i(i)*g))
        !!...........!!
        ! xi à gauche de la détente ou à droite de la détente !
        !!...........!!
        else
            hexa_i(i) = hexa_i(i)
            uexa_i(i) = uexa_i(i)
        end if
        
    end do
 
    
end subroutine sol_exa_tn

   !!...........!!
    ! Evolution de la détente dans le temps !
    ! On calcule la valeur de l'intervalle en fonction du temps !
    !!...........!!

    !calcul à l'instant t!
    !!
!subroutine detente(arg1,  arg2)
  !  type1, intent(in) :: arg1
 !   type2, intent(out) ::  arg2

    
!end subroutine detente

!écriture dans un fichier!
!subroutine ecriture_detente(name_file, tmax)
 !   type1, intent(in) :: arg1
  !  type2, intent(out) ::  arg2

 
!end subroutine ecriture_detente
    
end module mod_solexa