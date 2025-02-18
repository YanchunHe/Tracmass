SUBROUTINE setup_grid
    ! =============================================================
    ! Set up the grid
    ! =============================================================
    ! Subroutine for defining the grid of the GCM. Run once
    ! before the loop starts.
    ! -------------------------------------------------------------
    ! The following arrays have to be populated:
    !
    !  dxdy - Area of horizontal cell-walls.
    !  dzt  - Height of k-cells in 4 dim
    !  kmt  - Number of k-cells from surface to seafloor.
    !
    ! The following might be needed to calculate
    ! dxdy, uflux, and vflux
    !
    !  dzu - Height of each u-gridcell.
    !  dzv - Height of each v-gridcell.
    !  dxv -
    !  dyu -
    ! -------------------------------------------------------------

    USE mod_precdef
    USE mod_param
    USE mod_getfile
    USE mod_grid
    USE mod_seedvars

    IMPLICIT none

    !INTEGER               :: jj, ii, ip, jp, kk

    !INTEGER, ALLOCATABLE, DIMENSION(:,:)       :: tmp2d
    !REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)    :: dzt
    !REAL(DP), ALLOCATABLE, DIMENSION(:,:)      :: dyu,dxv,dxdy

    ! dx and dy in u and v points
    dyu = get2DfieldNC(TRIM(topoDataDir)//TRIM(hgridFile), dyu_name,[imindom,jmindom,1,1],[imt,jmt,1,1],'st')
    dxv = get2DfieldNC(TRIM(topoDataDir)//TRIM(hgridFile), dxv_name,[imindom,jmindom,1,1],[imt,jmt,1,1],'st')

    ! Grid area
    dxdy = get2DfieldNC(TRIM(topoDataDir)//TRIM(hgridFile), 'parea',[imindom,jmindom,1,1],[imt,jmt,1,1],'st')

!   write(*,*) 'min(dxdy):',minval(dxdy)
!   write(*,*) 'max(dxdy):',maxval(dxdy)

    ! Sensitivity to dz value
    !ALLOCATE(dzt(imt,jmt,km,1))
    !tmp3d(:,:,:) = get3DfieldNC(TRIM(topoDataDir)//TRIM(zgridFile), dzt_name,[imindom,jmindom,1,1],[imt,jmt,km,1],'st')
    !dzt(:,:,:,1) = get3DfieldNC(TRIM(physDataDir)//TRIM(zgridFile), dzt_name,[imindom,jmindom,1,1],[imt,jmt,km,1],'st')

END SUBROUTINE setup_grid
