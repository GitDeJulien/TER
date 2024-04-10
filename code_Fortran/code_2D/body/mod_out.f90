!-------------------------------------------------------------------
! module contenant la subroutine :
!  sortie
!
!  VERSION
!  	15/10/2022 par Luc Mieussens
!-------------------------------------------------------------------


module mod_out

    implicit none

contains

    !-------------------------------------------------------------------
    ! sortie dans un fichier au format ascii vtk, dataset unstructured grid,
    ! des valeurs contenues dans le tableau T (une valeur par maille)
    !-------------------------------------------------------------------
    subroutine out(iter, U, coord_noeud, noeud_maille, choix_sortie)

        use mod_constantes
        
        !--- entrees
        integer, intent(in) :: iter, choix_sortie
        real(pr), dimension(:,:), intent(in) :: U
        real(pr), dimension(:,:), intent(in) :: coord_noeud
        integer, dimension(:,:), intent(in) :: noeud_maille

        !--- locales
        character(len=30) :: it_ch
        integer :: i, nb_noeuds, nb_mailles
        

        write(it_ch,*) iter
        select case(choix_sortie)
        case(1)
            open(unit=20,file='out/out_'//trim(adjustl(it_ch))//'.vtk')

        end select
        nb_noeuds = size(coord_noeud,1)
        nb_mailles = size(U,1)

        write(20,'(1A26)') '# vtk DataFile Version 2.0'
        write(20,*) '2D Unstrucured Grid'
        write(20,*) 'ASCII'
        write(20,*) 'DATASET UNSTRUCTURED_GRID'

        write(20,*) 'POINTS',nb_noeuds,' double'
        do i=1,nb_noeuds
            write(20,*) coord_noeud(i,1),coord_noeud(i,2),0.0_pr
        end do

        write(20,*) 'CELLS ',nb_mailles,nb_mailles*4
        do i=1,nb_mailles
            write(20,*) 3,noeud_maille(i,1)-1, noeud_maille(i,2)-1, noeud_maille(i,3)-1
        end do

        write(20,*) 'CELL_TYPES ',nb_mailles
        do i=1,nb_mailles
            write(20,*) 5
        end do

        write(20,*) 'CELL_DATA',nb_mailles
        write(20,*) 'SCALARS Hauteur double'
        write(20,*) 'LOOKUP_TABLE default'
        do i=1,nb_mailles
            write(20,*) U(i,1)
        end do
        write(20,*)
        write(20,*) 'SCALARS xVelocity double'
        write(20,*) 'LOOKUP_TABLE default'
        do i=1,nb_mailles
            write(20,*) U(i,2)/U(i,1)
        end do
        write(20,*)
        write(20,*) 'SCALARS yVelocity double'
        write(20,*) 'LOOKUP_TABLE default'
        do i=1,nb_mailles
            write(20,*) U(i,3)/U(i,1)
        end do


        close(20)

    end subroutine out

end module mod_out