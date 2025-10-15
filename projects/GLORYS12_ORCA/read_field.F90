SUBROUTINE read_field

  !==========================================================================
  !
  ! Purpose
  ! -------
  !
  ! Read test model output to advect trajectories.
  ! Will be called by the loop each time step.
  !
  ! Method
  ! ------
  !
  ! Read velocities and optionally some tracers from netCDF files and
  ! update velocity fields for TRACMASS.
  !
  ! Updates the variables:
  !   uflux and vflux
  ! ==========================================================================


   USE mod_precdef
   USE mod_param
   USE mod_vel
   USE mod_time
   USE mod_grid
   USE mod_getfile
   USE mod_tracervars
   USE mod_tracers
   USE mod_calendar
   USE mod_swap

   USE netcdf

   IMPLICIT none

   INTEGER        :: kk, itrac

   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)  :: tmp3d
   CHARACTER (len=200)                      :: fieldFile, dateprefix

   ! Reassign the time index of uflux and vflux, dzt, ...
   CALL swap_time()

   ! Data files
   dateprefix = ' '

   ! Reading 2-time step variables
   ! In this case: velocities and tracers
   ! ===========================================================================
   ALLOCATE(tmp3d(imt,jmt,km))

   nctstep = currMon
   IF (l_onestep) nctstep = 1

   dateprefix = filledFileName(dateFormat, currYear, currMon, currDay)

   fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(uGridName)//TRIM(dateprefix)//TRIM(fileSuffix)
   uvel(1:imt,1:jmt,km:1:-1) = get3DfieldNC(fieldFile, ueul_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')

   IF (usgs_name/='') THEN
     tmp3d(1:imt,1:jmt,km:1:-1)  = get3DfieldNC(fieldFile, usgs_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')
     uvel(1:imt,1:jmt,1:km) = uvel(1:imt,1:jmt,1:km) + tmp3d(1:imt,1:jmt,1:km)
   END IF

   fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(vGridName)//TRIM(dateprefix)//TRIM(fileSuffix)
   vvel(1:imt,1:jmt,km:1:-1) = get3DfieldNC(fieldFile, veul_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')

   IF (vsgs_name/='') THEN
     tmp3d(1:imt,1:jmt,km:1:-1) = get3DfieldNC(fieldFile, vsgs_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')
     vvel(1:imt,1:jmt,1:km) = vvel(1:imt,1:jmt,1:km) + tmp3d(1:imt,1:jmt,1:km)
   END IF

   ! vertical flux
#if defined w_explicit
     
     ! Manually allocate wvel
     IF (.not. ALLOCATED(wvel)) THEN
       ALLOCATE ( wvel(imt,jmt,km))
       wvel(:,:,:) = 0.
     ENDIF
     
     !PRINT *,dateprefix
     fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(wGridName)//TRIM(dateprefix)//TRIM(fileSuffix)
     wvel(1:imt,1:jmt,km:1:-1) = get3DfieldNC(fieldFile, w_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')
#endif
   
   !! Tracers
   IF (l_tracers) THEN

     DO itrac = 1, numtracers

        ! Make sure the data array is empty
        tmp3d(:,:,:) = 0.d0

        IF (tracers(itrac)%action == 'read') THEN

            ! Read the tracer from a netcdf file
            IF (tracers(itrac)%varname == 'somxl010') THEN
                fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//'grid2D_'//TRIM(dateprefix)//TRIM(fileSuffix)
            ELSE
                fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(TGridName)//TRIM(dateprefix)//TRIM(fileSuffix)
            ENDIF
            
            IF (tracers(itrac)%dimension == '3D') THEN
                tmp3d(1:imt,1:jmt,km:1:-1) = get3DfieldNC(fieldFile, tracers(itrac)%varname,[imindom,jmindom,1,nctstep] &
                          ,[imt,jmt,km,1],'st')
            ELSE IF (tracers(itrac)%dimension == '2D') THEN
                tmp3d(1:imt,1:jmt,1) = get2DfieldNC(fieldFile, tracers(itrac)%varname,[imindom,jmindom,nctstep,1] &
                                        ,[imt,jmt,1,1],'st')
            END IF

        ELSE IF (tracers(itrac)%action == 'compute') THEN

            ! Compute the tracer from a function defined in mod_tracer.F90
            CALL compute_tracer(tracers(itrac)%name, tmp3d(1:imt,1:jmt,1:km))

        ELSE
            PRINT '(A34,I4)', 'No action defined for this tracer:', itrac
            STOP 10
        END IF

        ! Store the information
        IF (tracers(itrac)%dimension == '3D') THEN
          tracers(itrac)%data(:,1:jmt,:,2) = tracers(itrac)%scale*tmp3d(:,:,:) + tracers(itrac)%shift
        ELSE IF (tracers(itrac)%dimension == '2D') THEN
          tracers(itrac)%data(:,1:jmt,1,2) = tracers(itrac)%scale*tmp3d(:,:,1) + tracers(itrac)%shift
        END IF

     END DO
   END IF

   ! ===========================================================================

   ! uflux and vflux computation
   FORALL (kk = 1:km) uflux(:,:,kk,2)     = uvel(:,:,kk)*dyu(:,:)*dzu(:,:,kk,2)
   FORALL (kk = 1:km) vflux(:,1:jmt,kk,2) = vvel(:,1:jmt,kk)*dxv(:,1:jmt)*dzv(:,1:jmt,kk,2)
   
   WHERE (ISNAN(vflux)) vflux = 0.
   WHERE (ISNAN(uflux)) uflux = 0.
#if defined w_explicit
     FORALL (kk = 1:km) wflux(:,:,kk,2)     = wvel(:,:,kk)*dxdy(:,:)
     WHERE (ISNAN(wflux)) wflux = 0.
#endif

   !! Zero meridional flux at j=0 and j=jmt
   vflux(:,0  ,:,:) = 0.d0
   vflux(:,jmt,:,:) = 0.d0

   ! Reverse the sign of fluxes if trajectories are run backward in time.
   CALL swap_sign()

END SUBROUTINE read_field
