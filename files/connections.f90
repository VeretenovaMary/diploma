! --------------------------------------------------------------------------
!subroutine ConnectFaceToFace(pmesh1,pmesh2,connector)
subroutine ConnectFaceToFace(pmesh1,pmesh2,connector,mapper)
!           meshInterpNodeByCell
!           meshInterpCellByNode
!           meshInterpNodeByNode
    use module_settings
    use module_geom
    use module_Mesh
    implicit none
    ! input:
    type(TMeshData) :: pmesh1    ! pull mesh (acceptor)
    type(TMeshData) :: pmesh2    ! push mesh (source)
    ! output:
    integer, dimension(pmesh1%NFaces) :: connector
    ! optional:
    optional :: mapper
    ! local:
    integer :: iface1,iface2,cface
    real(r8p), dimension(3) :: crd1,crd2
    real(r8p), dimension(3) :: pt1,pt2
    real(r8p) :: eps,dist

    do iface1 = 1,pmesh1%NFaces
        dist  = 1.0e10_r8p        ! distance between face centers
        cface = 0                 ! connected face index (searching)
        crd1 = pmesh1%Face(iface1)%center
        pt1  = crd1
        do iface2 = 1,pmesh2%NFaces
            crd2   = pmesh2%Face(iface2)%center
            pt2(1) = crd2(1)
            pt2(2) = 0.0e0_r8p
            pt2(3) = 0.0e0_r8p
            if(present(mapper)) then
                call mapper(crd2,pt2)
            endif
            eps = point_point_distance(pt1,pt2)
            if(eps.lt.dist) then
                dist  = eps
                cface = iface2
            end if
        end do

        if(cface.gt.0) then
            connector(iface1) = cface
        else
            write(*,*) " face connection has not been found"
            stop
        end if
    end do

    return
end subroutine ConnectFaceToFace

