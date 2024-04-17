module mod_maillage

    use mod_constantes

    implicit none

contains

    subroutine MeshGenerator(MeshFile, NumberOfCells, NumberOfEdges, NodeCoord, CellVertices, &
        ClEdge, EdgeLength, CellCenterCoord, CellArea, CellPerimeter, NormalVectCoord, &
        EdgeNeighbor, EdgeOfCell)

        character(len=30), intent(in) :: MeshFile
        
        !Geometry
        integer, intent(out) :: NumberOfCells, NumberOfEdges
        real(pr), dimension(:), allocatable, intent(out) :: CellArea, CellPerimeter, EdgeLength
        real(pr), dimension(:,:), allocatable, intent(out) :: NodeCoord, CellCenterCoord, NormalVectCoord
        integer, dimension(:,:), allocatable, intent(out) :: CellVertices

        !Connectivity
        integer, dimension(:,:), allocatable, intent(out) :: EdgeNeighbor, EdgeOfCell 

        !CL
        integer, dimension(:), allocatable, intent(out) :: ClEdge

        !Local
        integer :: i, j
        integer :: NumberOfNodes, NumberOfSideEdges
        integer, dimension(:), allocatable :: ClSideEdge
        integer, dimension(:,:), allocatable :: SideEdgeVertices, EdgeVertices


        print*, "Mesh reading and creating connectivities..."
        call MeshReader(MeshFile, NumberOfNodes, NumberOfCells, NodeCoord, SideEdgeVertices, ClSideEdge, CellVertices)
        NumberOfSideEdges = size(ClSideEdge)

        print*,"-----------------------------------------"
        print*, "Number of nodes :", NumberOfNodes
        print*, "Number of cells :", NumberOfCells


        call Connectivity(NodeCoord, CellVertices, EdgeVertices, EdgeOfCell, EdgeNeighbor)
        NumberOfEdges = size(EdgeVertices, 1)
        print*, "Number of edges :", NumberOfEdges
        print*,"-----------------------------------------"
        print*, "... done"

        print*, "Centers, lengths and areas computing ..."

        !-- Edge conditions
        allocate(ClEdge(1:NumberOfEdges))
        do i=1,NumberOfEdges
            if (EdgeNeighbor(i,1)==0 .OR. EdgeNeighbor(i,2)==0) then

                do j=1,NumberOfSideEdges

                    if ( (SideEdgeVertices(j,1)==EdgeVertices(i,1) .AND. SideEdgeVertices(j,2)==EdgeVertices(i,2)) .OR. &
                    (SideEdgeVertices(j,1)==EdgeVertices(i,2) .AND. SideEdgeVertices(j,2)==EdgeVertices(i,1))) then

                        ClEdge(i) = ClSideEdge(j) 

                    end if
                end do
            end if
        end do

        ! -- Cell areas
        call AreasPerimeterCenterComputing(NumberOfCells, CellVertices, NodeCoord, CellCenterCoord, CellPerimeter, CellArea)

        call EdgeLengthComputing(NumberOfEdges, EdgeVertices, NodeCoord, EdgeLength)

        print*, "... done"
        print*, "Normal computing ..."

        call NormalComputing(NumberOfEdges, EdgeNeighbor, EdgeVertices, CellCenterCoord, NodeCoord, NormalVectCoord)

        print*, "... done"
        print*, "-------------------------------------------"
        print*, "The mesh have been well readed and all the necessary vectors have been generated."
        print*, "-------------------------------------------"



    end subroutine

! ###############################################################################################################

    subroutine MeshReader(MeshFile,NumberOfNodes, NumberOfCells, NodeCoord, SideEdgeVertices, ClSideEdge, CellVertices)

        ! -- In
        character(len=30), intent(in) :: MeshFile

        ! -- Out
        integer, intent(out) :: NumberOfCells, NumberOfNodes
        integer, dimension(:), allocatable, intent(out) :: ClSideEdge
        real(pr), dimension(:,:), allocatable, intent(out) :: NodeCoord
        integer, dimension(:,:), allocatable, intent(out) :: SideEdgeVertices, CellVertices

        !Local
        integer :: i, NumberOfSideEdges
        

        open(unit=10, file=trim(adjustl(MeshFile)))

        read(10,*)
        read(10,*)
        read(10,*)

        ! -- Vertices
        read(10,*)
        read(10,*) NumberOfNodes
        allocate(NodeCoord(1:NumberOfNodes,1:2))

        do i=1,NumberOfNodes
            read(10,*) NodeCoord(i,1), NodeCoord(i,2)
        end do

        ! -- Edges (here it is the edges of the domain)
        read(10,*)
        read(10,*) NumberOfSideEdges
        allocate(SideEdgeVertices(NumberOfSideEdges, 2), ClSideEdge(NumberOfSideEdges))

        do i=1,NumberOfSideEdges
            read(10,*) SideEdgeVertices(i,1), SideEdgeVertices(i,2), ClSideEdge(i)
        end do

        ! -- Cell
        read(10,*)
        read(10,*) NumberOfCells
        allocate(CellVertices(NumberOfCells, 3))

        do i=1,NumberOfCells
            read(10,*) CellVertices(i,1), CellVertices(i,2), CellVertices(i,3)
        end do

        
        close(10)

    end subroutine MeshReader

! ###########################################################################
    
    subroutine Connectivity(NodeCoord, CellVertices, EdgeVertices, EdgeOfCell, EdgeNeighbor)

        ! -- In
        real(pr), dimension(:,:), intent(in) :: NodeCoord
        integer, dimension(:,:), allocatable, intent(in) :: CellVertices

        ! -- Out
        integer, dimension(:,:), allocatable, intent(out) :: EdgeVertices, EdgeNeighbor, EdgeOfCell

        ! -- Local
        integer :: i, v1, v2, v3, cpt_edges, NbNodes, NbCells
        integer, dimension(:,:), allocatable :: flag

        NbNodes = size(NodeCoord,1)
        NbCells = size(CellVertices,1)

        cpt_edges = 0
        allocate(flag(1:NbNodes, 1:NbNodes))
        flag = 0
        do i=1,NbCells
            v1 = CellVertices(i,1)
            v2 = CellVertices(i,2)
            v3 = CellVertices(i,3)
            if (flag(v1,v2)==0 .OR. flag(v2,v1)==0) then
                cpt_edges = cpt_edges + 1
                flag(v1,v2)= cpt_edges
                flag(v2,v1)= cpt_edges
            end if
            if (flag(v2,v3)==0 .OR. flag(v3,v2)==0) then
                cpt_edges = cpt_edges + 1
                flag(v2,v3)= cpt_edges
                flag(v3,v2)= cpt_edges
            end if
            if (flag(v3,v1)==0 .OR. flag(v1,v3)==0) then
                cpt_edges = cpt_edges + 1
                flag(v3,v1)= cpt_edges
                flag(v1,v3)= cpt_edges
            end if
        end do

        allocate(EdgeVertices(1:cpt_edges,2))
        do i=1,NbCells
            v1 = CellVertices(i,1)
            v2 = CellVertices(i,2)
            v3 = CellVertices(i,3)
            EdgeVertices(flag(v1,v2),1)=v1
            EdgeVertices(flag(v1,v2),2)=v2
            EdgeVertices(flag(v2,v3),1)=v2
            EdgeVertices(flag(v2,v3),2)=v3
            EdgeVertices(flag(v3,v1),1)=v3
            EdgeVertices(flag(v3,v1),2)=v1
        end do

        allocate(EdgeNeighbor(1:cpt_edges,2))
        allocate(EdgeOfCell(1:NbCells,3))
        EdgeNeighbor = 0
        do i=1,NbCells
            v1 = CellVertices(i,1)
            v2 = CellVertices(i,2)
            v3 = CellVertices(i,3)
            EdgeOfCell(i,1) = flag(v1,v2)
            if (EdgeNeighbor(EdgeOfCell(i,1),1)/=0) then
                EdgeNeighbor(EdgeOfCell(i,1),2)=i
            else
                EdgeNeighbor(EdgeOfCell(i,1),1)=i
            end if

            EdgeOfCell(i,2) = flag(v2,v3)
            if (EdgeNeighbor(EdgeOfCell(i,2),1)/=0) then
                EdgeNeighbor(EdgeOfCell(i,2),2)=i
            else
                EdgeNeighbor(EdgeOfCell(i,2),1)=i
            end if

            EdgeOfCell(i,3) = flag(v3,v1)
            if (EdgeNeighbor(EdgeOfCell(i,3),1)/=0) then
                EdgeNeighbor(EdgeOfCell(i,3),2)=i
            else
                EdgeNeighbor(EdgeOfCell(i,3),1)=i
            end if
        end do

        deallocate(flag)

    end subroutine Connectivity


! ##########################################################################

    subroutine EdgeLengthComputing (NumberOfEdges, EdgeVertices, NodeCoord, EdgeLength)

        ! -- In
        integer, intent(in) :: NumberOfEdges
        integer, dimension(:,:), intent(in) :: EdgeVertices 
        real(pr), dimension(:,:), intent(in) :: NodeCoord

        ! -- Out
        real(pr), dimension(:), allocatable, intent(out) :: EdgeLength
        
        ! -- Local
        integer :: i!, v1, v2
        real(pr),dimension(2) :: A, B

        allocate(EdgeLength(1:NumberOfEdges))

        do i=1,NumberOfEdges
            ! v1 = EdgeVertices(i,1)
            ! v2 = EdgeVertices(i,2)

            A = NodeCoord(EdgeVertices(i,1), :)
            B = NodeCoord(EdgeVertices(i,2), :)

            !EdgeLength(i) = sqrt((NodeCoord(v1,1)-NodeCoord(v2,1))**2 + (NodeCoord(v1,2)-NodeCoord(v2,2))**2)

            EdgeLength(i) = sqrt(sum((B-A)**2))
        end do

    end subroutine EdgeLengthComputing 

! ##########################################################################

    subroutine AreasPerimeterCenterComputing (NumberOfCells, CellVertices, NodeCoord, CellCenterCoord,&
        CellPerimeter, CellArea)

        ! -- In
        integer, intent(in) :: NumberOfCells
        integer, dimension(:,:), intent(in):: CellVertices
        real(pr), dimension(:,:), intent(in) :: NodeCoord

        ! -- Out
        real(pr), dimension(:), allocatable, intent(out) :: CellArea, CellPerimeter
        real(pr), dimension(:,:), allocatable, intent(out) :: CellCenterCoord

        ! -- Local
        integer :: i, v1, v2, v3
        real(pr), dimension(2) :: A, B, C
        real(pr) :: diff, RA, RB, RC, dx, dy

        allocate(CellArea(1:NumberOfCells), CellPerimeter(1:NumberOfCells))
        allocate(CellCenterCoord(1:NumberOfCells, 2))

        do i=1,NumberOfCells
            v1 = CellVertices(i,1)
            v2 = CellVertices(i,2)
            v3 = CellVertices(i,3)

            A = NodeCoord(v1,:)
            B = NodeCoord(v2,:)
            C = NodeCoord(v3,:)

            diff = 2*( (B(1)-A(1))*(C(2)-A(2)) - (B(2)-A(2))*(C(1)-A(1)) ) 

            RA = sum(A**2); RB = sum(B**2); RC = sum(C**2)
            dx = (RB - RA)*(C(2) - A(2)) - (RC - RA)*(B(2) - A(2))
            dy = (RB - RA)*(C(1) - A(1)) - (RC - RA)*(B(1) - A(1))

            CellCenterCoord(i,1) = dx/diff
            CellCenterCoord(i,2) = - dy/diff

            CellPerimeter(i) = sqrt(sum(B-A)**2) + sqrt(sum(C-B)**2) + sqrt(sum(A-C)**2)
            CellArea(i) = abs(((B(1) - A(1))*(C(2)- A(2))) - ((B(2) - A(2))*(C(1) - A(1)))) / 2

        end do

    end subroutine AreasPerimeterCenterComputing


! ######################################################################

    subroutine NormalComputing (NumberOfEdges, EdgeNeighbor, EdgeVertices, CellCenterCoord,&
        NodeCoord, NormalVectCoord)

        ! -- In
        integer, intent(in) :: NumberOfEdges
        integer, dimension(:,:), intent(in) :: EdgeNeighbor, EdgeVertices
        real(pr), dimension(:,:), intent(in) :: CellCenterCoord, NodeCoord

        ! -- Out
        real(pr), dimension(:,:), allocatable, intent(out) :: NormalVectCoord

        ! -- Local
        integer :: i, n1, n2, A, B
        real(pr) :: x1, x2, y1, y2, norm
        real(pr), dimension(2) :: XiXk, XiXm
        real(pr), dimension(2) :: Xi, Xk, Xm

        allocate(NormalVectCoord(1:NumberOfEdges, 2))


        do i=1,NumberOfEdges
            n1 = EdgeVertices(i,1)
            n2 = EdgeVertices(i,2)

            x1 = NodeCoord(n1,1)
            y1 = NodeCoord(n1,2)
            x2 = NodeCoord(n2,1)
            y2 = NodeCoord(n2,2)

            A = EdgeNeighbor(i,1)
            B = EdgeNeighbor(i,2)

            Xi = CellCenterCoord(A,:)

            NormalVectCoord(i,1) = y1 - y2
            NormalVectCoord(i,2) = x2 - x1

            if (B/=0) then
                Xk = CellCenterCoord(B,:)
                XiXk = Xk - Xi

                if (NormalVectCoord(i,1)*XiXk(1) + NormalVectCoord(i,2)*XiXk(2) < 0) then
                    NormalVectCoord(i,:) = - NormalVectCoord(i,:)
                end if
            else
                Xm(1) = (x1 + x2) / 2
                Xm(2) = (y1 + y2) / 2

                XiXm = Xm - Xi

                if (NormalVectCoord(i,1)*XiXm(1) + NormalVectCoord(i,2)*XiXm(2) < 0) then
                    NormalVectCoord(i,:) = - NormalVectCoord(i,:)
                end if

            end if
            norm = sqrt((NormalVectCoord(i,1)**2+NormalVectCoord(i,2)**2))
            ! print*, "norm", norm

            NormalVectCoord(i,:) =  NormalVectCoord(i,:)/norm
        end do

    end subroutine



end module

