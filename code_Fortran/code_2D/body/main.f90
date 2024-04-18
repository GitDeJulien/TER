program main


    use mod_constantes
    use mod_maillage
    use mod_flux
    use mod_out

    implicit none

    ! -- Mesh informations
    character(len=30)                      :: MeshFile
    !character(len=30)                      :: chiter
    character(len=40)                      :: ct
    integer                                :: NumberOfCells, NumberOfEdges
    integer, dimension(:), allocatable     :: ClEdge
    integer, dimension(:,:), allocatable   :: CellVertices, EdgeNeighbor, EdgeOfCell
    real(pr), dimension(:), allocatable    :: CellArea, CellPerimeter, EdgeLength
    real(pr), dimension(:,:), allocatable  :: NodeCoord, CellCenterCoord, NormalVectCoord

    ! -- Flux tools
    integer                               :: e, i, k, iter, nplot, nmax, cpt
    real(pr), dimension(:), allocatable   :: be, Fik
    real(pr), dimension(:,:), allocatable :: Un, Unp1

    ! -- TimeLoop tools
    real(pr) :: t, tmax, dt, cfl

    MeshFile = "maillage/maillage_0.04.mesh"

    call MeshGenerator(MeshFile, NumberOfCells, NumberOfEdges, NodeCoord, CellVertices, &
    ClEdge, EdgeLength, CellCenterCoord, CellArea, CellPerimeter, NormalVectCoord, &
    EdgeNeighbor, EdgeOfCell)

    allocate(Unp1(1:NumberOfCells,3))

    ! -- Time Loop
    tmax = 10. !seconde
    dt = 0.
    t = 0.
    cfl = 0.8
    iter = 0
    cpt = 0

    do while (t < tmax)

        if (t < 10e-5) then
            ! -- U initialization
            Un = InitU(NumberOfCells, CellCenterCoord)
        end if

        Unp1 = Un
        print*, "t=",t

        ! -- be computing
        be = beComputing(NumberOfEdges, EdgeNeighbor, Un, NormalVectCoord)
        !print*, be

        ! -- dt computing
        dt = dtComputng(cfl, be, CellArea, CellPerimeter)
        nmax = int(tmax/dt)
        nplot = int(nmax/200.)
        !print*, "nplot=",nplot

        do e = 1,NumberOfEdges

            ! -- Flux computing
            call FluxNum(EdgeNeighbor, ClEdge, Un, NormalVectCoord, be, Fik, e, EdgeLength)

            i = EdgeNeighbor(e, 1)
            k = EdgeNeighbor(e, 2)

            if (k/=0) then
                Unp1(k,:) = Unp1(k,:) + dt/CellArea(k) * Fik(:)
            end if

            Unp1(i,:) = Unp1(i,:) -  dt/CellArea(i) * Fik(:)

        end do !end edge loop

        !! Ecriture des résultat format vtk
        ! if (mod(iter,nplot) == 0) then
        !     call out(iter, Un, NodeCoord, CellVertices, 1)
        ! end if
        ! ----------------------------------------
        !! Ecriture des résultats dans un .dat
        if (mod(iter,nplot) == 0) then
            write(ct,'(I4)') cpt
            open(unit=20, file='out/sol.'//trim(adjustl(ct))//'.dat')
                do i=1,size(Un,1)
                    write(20,*) CellCenterCoord(i,1), CellCenterCoord(i,2), Un(i,1), Un(i,2)/Un(i,1), Un(i,3)/Un(i,1)
                end do
            close(20)
            cpt = cpt + 1
        end if
        !----------------------------------------

        Un = Unp1
        t = t + dt
        iter = iter + 1

    end do !end time loop
    print*,"dt=", dt

    deallocate(NodeCoord, CellVertices, &
    ClEdge, EdgeLength, CellCenterCoord, CellArea, CellPerimeter, NormalVectCoord, &
    EdgeNeighbor, EdgeOfCell)
    deallocate(be, Un, Fik, Unp1)


end program
