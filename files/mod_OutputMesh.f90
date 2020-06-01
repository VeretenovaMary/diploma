module module_OutputMesh
    use module_Mesh

    integer :: iPtGmsh = 0
    integer :: iLnGmsh = 0

contains


! ---------------------------------------------------------------
    subroutine output_cell_gmsh(funit,icell,mesh,title)
        implicit none
        ! input
        integer                 :: funit
        integer                 :: icell
        type(TMeshData),pointer :: mesh
        ! local
        integer :: i,iface
        character(*), optional :: title

        if(present(title)) then
            continue
        else
            title = ""        
        endif

        write(funit,"(a,a)") '//+ begin cell: ',trim(title)
        do i = 1,mesh%Cell(icell)%nfaces
            iface = mesh%Cell(icell)%iface(i)
            call output_face_gmsh(funit,iface,mesh)
        end do
        write(funit,"(a,a)") '//+ done cell: ',trim(title)

        return
    end subroutine output_cell_gmsh

! ---------------------------------------------------------------
    subroutine output_face_gmsh(funit,iface,mesh,title)
        implicit none
        ! input
        integer                 :: funit
        type(TMeshData),pointer :: mesh
        ! local
        integer       :: iface,inode
        integer       :: k,iflag,ipt0
        character(64) :: sfmt1,sfmt2
        character(*), optional :: title

        if(present(title)) then
            continue
        else
            title = ""        
        endif

        sfmt1 = "('Point(',i7,') = {',3(1x,1pE15.8e2,','),'1.0};')"
        sfmt2 = "('Line(',i7,') = {'i7,',',i7,'};')"

        iflag = 0
        ipt0  = iPtGmsh + 1
        do k = 1,mesh%Face(iface)%nnodes
            inode = mesh%Face(iface)%inode(k)

            iPtGmsh = iPtGmsh + 1
            write(funit,"(a,i7,a)") '//+ face point',k,trim(title)
            write(funit,sfmt1) iPtGmsh,mesh%Node(inode)%crd
            !write(funit,"(a)") '//+'
            if(iflag.eq.0) then
                iflag = iflag + 1
            else
                iLnGmsh = iLnGmsh + 1
                write(funit,"(a,i7,a)") '//+ face edge',k-1,trim(title)
                write(funit,sfmt2) iLnGmsh,iPtGmsh-1,iPtGmsh
            endif
        end do

        iLnGmsh = iLnGmsh + 1
        write(funit,"(a,i7,a,a)") '//+ face edge',mesh%Face(iface)%nnodes,' - last',trim(title)
        write(funit,sfmt2) iLnGmsh,iPtGmsh,ipt0

        return
    end subroutine output_face_gmsh

! ---------------------------------------------------------------
    subroutine output_line_gmsh(funit,pt1,pt2,title)
        implicit none
        ! input
        integer                 :: funit
        real(8), dimension(3)   :: pt1,pt2
        ! local
        character(64) :: sfmt1,sfmt2
        character(*), optional :: title

        if(present(title)) then
            continue
        else
            title = ""        
        endif

        sfmt1 = "('Point(',i7,') = {',3(1x,1pE15.8e2,','),'1.0};')"
        sfmt2 = "('Line(',i7,') = {'i7,',',i7,'};')"

        write(funit,"(a,a)") '//+ line point 1: ',trim(title)
        iPtGmsh = iPtGmsh + 1
        write(funit,sfmt1) iPtGmsh,pt1
        write(funit,"(a,a)") '//+ line point 2: ',trim(title)
        iPtGmsh = iPtGmsh + 1
        write(funit,sfmt1) iPtGmsh,pt2

        write(funit,"(a,a)") '//+ line from points: ',trim(title)
        iLnGmsh = iLnGmsh + 1
        write(funit,sfmt2) iLnGmsh,iPtGmsh-1,iPtGmsh

        return
    end subroutine output_line_gmsh

! ---------------------------------------------------------------
    subroutine output_point_gmsh(funit,pt,title)
        implicit none
        ! input
        integer                 :: funit
        real(8), dimension(3)   :: pt
        ! local
        character(64) :: sfmt1
        character(*), optional :: title

        if(present(title)) then
            continue
        else
            title = ""        
        endif

        sfmt1 = "('Point(',i7,') = {',3(1x,1pE15.8e2,','),'1.0};')"

        write(funit,"(a,a)") '//+ single point: ',trim(title)
        iPtGmsh = iPtGmsh + 1
        write(funit,sfmt1) iPtGmsh,pt

        return
    end subroutine output_point_gmsh

! -------------------------------------------------------------
    subroutine output_2DMesh_gmsh(mesh,FName)
        USE module_Mesh             ! KV (including additional grid data)

        IMPLICIT NONE
        ! input:
        type(TMeshData) :: mesh
        character(*)    :: FName
        ! local:
        integer :: i,j
        integer :: n1,n2,idx
!        integer  :: valpos_str
!        external :: valpos_str


        open(201,file=FName)
        write(201,"(a)") "$MeshFormat"
        write(201,"(a)") "2.2 0 8"
        write(201,"(a)") "$EndMeshFormat"
        write(201,"(a)") "$PhysicalNames"
        write(201,"(i2)") mesh%LocNBnds+mesh%LocNVols
        j = 0
        do i = 1,mesh%LocNBnds
            j = j + 1
            write(201,"(i2,i2,1x,a)") 2,mesh%LocBndIdx(i),mesh%LocBndTit(i)
!            write(201,"(i2,i2,1x,a)") 2,j,mesh%LocBndTit(i)
        enddo
        do i = 1,mesh%LocNVols
            j = j + 1
            write(201,"(i2,i2,1x,a)") 2,mesh%LocVolIdx(i),mesh%LocVolTit(i)
!            write(201,"(i2,i2,1x,a)") 2,j,mesh%LocVolTit(i)
        enddo
        write(201,"(a)") "$EndPhysicalNames"
        ! write nodes
        write(201,"(a)") "$Nodes"
        write(201,"(i5)") mesh%NNodes
        do i = 1,mesh%NNodes
            write(201,"(i5,3(1x,1pE12.5e2))") i,mesh%Node(i)%crd(1),mesh%Node(i)%crd(2),0.e0
        enddo
        write(201,"(a)") "$EndNodes"
        ! write elements
        write(201,"(a)") "$Elements"
        write(201,"(i5)") mesh%NBFaces+mesh%NCells
        j = 0
        do i = 1,mesh%NBFaces
            j = j + 1
            n1 = mesh%Face(i)%inode(1)
            n2 = mesh%Face(i)%inode(2)
!            idx = 1
!            do while(mesh%LocBndIdx(idx).ne.mesh%Face(i)%locIdx)
!                idx = idx + 1
!            enddo
            idx = mesh%Face(i)%locIdx
            write(201,"(i5,4(1x,i2),2(1x,i5))") j,1,2,idx,10+idx,n1,n2
        enddo

        do i = 1,mesh%NCells
            j = j + 1
!            idx = 5
            idx = mesh%Cell(i)%locIdx
            write(201,"(i5,4(1x,i2),4(1x,i5))") j,3,2,idx,10+idx,mesh%Cell(i)%inode
        enddo
        write(201,"(a)") "$EndElements"


        close(201)

        return
    end subroutine output_2DMesh_gmsh


end module module_OutputMesh
