module mod_flux

    use mod_constantes

    implicit none


    contains


    subroutine FluxNum(EdgeNeighbor, ClEdge, U, ne, be, Fik, e, EdgeLength)

        ! -- In
        integer, intent(in) :: e
        real(pr), dimension(:), allocatable, intent(in)   :: EdgeLength, be
        real(pr), dimension(:,:), allocatable, intent(in) :: U
        real(pr), dimension(:,:), allocatable, intent(in) :: ne
        integer, dimension(:,:), allocatable, intent(in)  :: EdgeNeighbor
        integer, dimension(:), allocatable, intent(in)   :: ClEdge

        ! -- Out
        real(pr), dimension(:), allocatable, intent(out) :: Fik

        ! -- Local
        integer :: i, k
        real(pr) :: l

        allocate(Fik(1:3))

        Fik = 0

        i = EdgeNeighbor(e, 1)
        k = EdgeNeighbor(e, 2)
        l = EdgeLength(e)

        if (k/=0) then

            Fik(1) = l*(0.5*(U(i,2)*ne(e,1) + U(i,3)*ne(e,2) + &
            U(k,2)*(ne(e,1)) + U(k,3)*(ne(e,2))) &
            - 0.5*be(e)*(U(k,1) - U(i,1)))

            Fik(2) =  l*(0.5*((U(i,2)*U(i,2)/U(i,1) + &
            gravity*0.5*U(i,1)*U(i,1))*ne(e,1) + &
            U(i,2)*U(i,3)/U(i,1)*ne(e,2) + &
            (U(k,2)*U(k,2)/U(k,1) + &
            gravity*0.5*U(k,1)*U(k,1))*(ne(e,1)) + &
            U(k,2)*U(k,3)/U(k,1)*(ne(e,2))) - &
            0.5*be(e)*(U(k,2) - U(i,2)))

            Fik(3) = l*(0.5*((U(i,3)*U(i,3)/U(i,1) + &
            gravity*0.5*U(i,1)*U(i,1))*ne(e,2) + &
            U(i,2)*U(i,3)/U(i,1)*ne(e,1) + &
            (U(k,3)*U(k,3)/U(k,1) + &
            gravity*0.5*U(k,1)*U(k,1))*(ne(e,2)) + &
            U(k,2)*U(k,3)/U(k,1)*(ne(e,1))) - &
            0.5*be(e)*(U(k,3) - U(i,3)))



        else if (ClEdge(e) == 1) then !Bord Wall (Mure avec rebond)

            Fik(1) = l*(0.5*(U(i,2)*ne(e,1) + U(i,3)*ne(e,2)))

            Fik(2) =  l*(0.5*((U(i,2)*U(i,2)/U(i,1) + &
            gravity*0.5*U(i,1)*U(i,1))*ne(e,1) + &
            U(i,2)*U(i,3)/U(i,1)*ne(e,2) + &
            (gravity*0.5*U(i,1)*U(i,1))*(ne(e,1))) - &
            0.5*be(e)*(- U(i,2)))

            Fik(3) = l*(0.5*((U(i,3)*U(i,3)/U(i,1) + &
            gravity*0.5*U(i,1)*U(i,1))*ne(e,2) + &
            U(i,2)*U(i,3)/U(i,1)*ne(e,1) + &
            (gravity*0.5*U(i,1)*U(i,1))*(ne(e,2))) - &
            0.5*be(e)*(- U(i,3)))


        else if (ClEdge(e) == 2) then !Bord sortie (Neumann homogene)

            Fik(1) = l*(0.5*(U(i,2)*ne(e,1) + U(i,3)*ne(e,2) + &
            U(i,2)*(ne(e,1)) + U(i,3)*(ne(e,2))) &
            - 0.5*be(e)*(U(i,1) - U(i,1)))

            Fik(2) =  l*(0.5*((U(i,2)*U(i,2)/U(i,1) + &
            gravity*0.5*U(i,1)*U(i,1))*ne(e,1) + &
            U(i,2)*U(i,3)/U(i,1)*ne(e,2) + &
            (U(i,2)*U(i,2)/U(i,1) + &
            gravity*0.5*U(i,1)*U(i,1))*(ne(e,1)) + &
            U(i,2)*U(i,3)/U(i,1)*(ne(e,2))) - &
            0.5*be(e)*(U(i,2) - U(i,2)))

            Fik(3) = l*(0.5*((U(i,3)*U(i,3)/U(i,1) + &
            gravity*0.5*U(i,1)*U(i,1))*ne(e,2) + &
            U(i,2)*U(i,3)/U(i,1)*ne(e,1) + &
            (U(i,3)*U(i,3)/U(i,1) + &
            gravity*0.5*U(i,1)*U(i,1))*(ne(e,2)) + &
            U(i,2)*U(i,3)/U(i,1)*(ne(e,1))) - &
            0.5*be(e)*(U(i,3) - U(i,3)))


        else if (ClEdge(e) == 3) then !Bord entrée (Mur réflexif)

            Fik(1) = l*(0.5*(U(i,2)*ne(e,1) + U(i,3)*ne(e,2)))

            Fik(2) =  l*(0.5*((U(i,2)*U(i,2)/U(i,1) + &
            gravity*0.5*U(i,1)*U(i,1))*ne(e,1) + &
            U(i,2)*U(i,3)/U(i,1)*ne(e,2) + &
            (gravity*0.5*U(i,1)*U(i,1))*(ne(e,1))) - &
            0.5*be(e)*(- U(i,2)))

            Fik(3) = l*(0.5*((U(i,3)*U(i,3)/U(i,1) + &
            gravity*0.5*U(i,1)*U(i,1))*ne(e,2) + &
            U(i,2)*U(i,3)/U(i,1)*ne(e,1) + &
            (gravity*0.5*U(i,1)*U(i,1))*(ne(e,2))) - &
            0.5*be(e)*(- U(i,3)))

        end if

    end subroutine

    ! #######################################################################"

    function InitU(NumberOfCells, CellCenterCoord)result(U)

        ! -- In
        integer, intent(in) :: NumberOfCells
        real(pr), dimension(:,:), allocatable, intent(in) :: CellCenterCoord

        ! -- Out
        real(pr), dimension(:,:), allocatable :: U

        ! -- Local
        integer :: i

        allocate(U(1:NumberOfCells,3))

        print*, "Vector U initialisation ..."

        do i=1,NumberOfCells

            if (CellCenterCoord(i,1) < 0.) then
                U(i,1) = 0.08
                U(i,2) = 0.
                U(i,3) = 0.
            else
                U(i,1) = 0.02
                U(i,2) = 0.
                U(i,3) = 0.
            end if

        end do

        print*, "Height (meter) before the dam:", U(1,1)
        print*, "Height (meter) after the dam:", U(NumberOfCells,1)
        print*, "... done"
        print*, "----------------------------------------------------"

    end function

    ! ################################################################

    function dtComputng(cfl, be, CellArea, CellPerimeter)result(dt)

        ! -- In
        real(pr), intent(in)                            :: cfl
        real(pr), dimension(:), allocatable, intent(in) :: be
        real(pr), dimension(:), allocatable, intent(in) :: CellArea, CellPerimeter

        ! -- Out
        real(pr) :: dt

        ! -- Local
        real(pr) :: be_max, lambda_min

        be_max = maxval(be)
        lambda_min = minval(CellArea/CellPerimeter)

        dt = cfl*lambda_min/be_max

    end function

    ! ##################################################################

    function beComputing(NumberOfEdges, EdgeNeighbor, U, ne)result(b)

        ! -- In
        integer, intent(in) :: NumberOfEdges
        real(pr), dimension(:,:), allocatable, intent(in) :: U
        real(pr), dimension(:,:), allocatable, intent(in) :: ne
        integer, dimension(:,:), allocatable, intent(in)  :: EdgeNeighbor

        ! -- Out
        real(pr), dimension(:), allocatable :: b

        ! -- Local
        integer :: e, i, k

        allocate(b(1:NumberOfEdges))

        do e=1,NumberOfEdges

            i = EdgeNeighbor(e, 1)
            k = EdgeNeighbor(e, 2)

            if (k/=0) then
                b(e) = max(abs(U(i,2)/U(i,1)*ne(e,1) + U(i,3)/U(i,1)*ne(e,2) + sqrt(gravity*U(i,1))), &
                abs(U(i,2)/U(i,1)*ne(e,1) + U(i,3)/U(i,1)*ne(e,2) - sqrt(gravity*U(i,1))), &
                abs(U(k,2)/U(k,1)*(-ne(e,1)) + U(k,3)/U(k,1)*(-ne(e,2)) + sqrt(gravity*U(k,1))), &
                abs(U(k,2)/U(k,1)*(-ne(e,1)) + U(k,3)/U(k,1)*(-ne(e,2)) - sqrt(gravity*U(k,1))))
            else
                b(e) = max(abs(U(i,2)/U(i,1)*ne(e,1) + U(i,3)/U(i,1)*ne(e,2) + sqrt(gravity*U(i,1))), &
                abs(U(i,2)/U(i,1)*ne(e,1) + U(i,3)/U(i,1)*ne(e,2) - sqrt(gravity*U(i,1))))
            end if
        end do

    end function

    ! #####################################################################

    function topology(NumberOfCells, CellCenterCoord, shape)result(t)

        ! -- In
        integer, intent(in) :: NumberOfCells
        integer, intent(in) :: shape
        real(pr), dimension(:,:), intent(in) :: CellCenterCoord

        ! -- Out
        real(pr), dimension(:), allocatable :: t

        ! -- Local
        integer :: i

        allocate(t(1:NumberOfCells))
        
        if (shape == 1) then !GrandPave
            do i=1,NumberOfCells
                if (CellCenterCoord(i,1) > 0.94 .AND. CellCenterCoord(i,1) < 1.06) then
                    t(i) = 0.02
                else
                    t(i) = 0.0
                end if
            end do
        end if

    end function

end module
