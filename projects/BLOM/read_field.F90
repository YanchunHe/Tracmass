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
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)  :: tmpuflux, tmpvflux, tmptracer
   !REAL(DP), ALLOCATABLE, DIMENSION(:,:)    :: tmp2d
   CHARACTER (len=200)                      :: fieldFile, dateprefix

   ! Reassign the time index of uflux and vflux, dzt, dzdt, hs, ...
   CALL swap_time()

   ! Data files
   dateprefix = ' '

   ! Reading 3-time step variables
   ! In this case: hs and zstot
   ! ===========================================================================
   IF (ints == 0) THEN

     ! 1 - Past
     IF (loopYears) THEN

       nctstep = prevMon
       IF (l_onestep) nctstep = 1

       dateprefix = filledFileName(dateFormat, prevYear, prevMon, prevDay)

       !fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(tGridName)//TRIM(dateprefix)//TRIM(fileSuffix)
       fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(dateprefix)//TRIM(fileSuffix)
       !WRITE(*,*) 'pass file:',fieldFile
       hs(1:imt,1:jmt,-1) = get2DfieldNC(fieldFile, hs_name,[imindom,jmindom,nctstep,1],[imt,jmt,1,1],'st')
       hs(imt+1,:,-1)     = hs(1,:,-1)

     END IF

     ! 2 - Present
     nctstep = currMon
     IF (l_onestep) nctstep = 1

     dateprefix = filledFileName(dateFormat, currYear, currMon, currDay)

     !fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(tGridName)//TRIM(dateprefix)//TRIM(fileSuffix)
     fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(dateprefix)//TRIM(fileSuffix)
       !WRITE(*,*) 'present file:',fieldFile
     hs(1:imt,1:jmt,0) = get2DfieldNC(fieldFile, hs_name,[imindom,jmindom,nctstep,1],[imt,jmt,1,1],'st')
     hs(imt+1,:,0)     = hs(1,:,0)

     WHERE (SUM(dzt(:,:,:,2),3) /= 0)
          zstot(1:imt,1:jmt,-1) = hs(1:imt,1:jmt,-1)/SUM(dzt(:,:,:,2),3) + 1
          zstot(1:imt,1:jmt, 0) = hs(1:imt,1:jmt, 0)/SUM(dzt(:,:,:,2),3) + 1
     ELSEWHERE
          zstot(:,:,-1) = 0.d0
          zstot(:,:, 0) = 0.d0
     END WHERE

   END IF

   ! 3 - Future
   IF (ints<intrun-1 .OR. loopYears) THEN

     nctstep = nextMon
     IF (l_onestep) nctstep = 1

     dateprefix = filledFileName(dateFormat, nextYear, nextMon, nextDay)

     !fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(tGridName)//TRIM(dateprefix)//TRIM(fileSuffix)
     fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(dateprefix)//TRIM(fileSuffix)
       !WRITE(*,*) 'future file:',fieldFile
     hs(1:imt,1:jmt,1) = get2DfieldNC(fieldFile, hs_name,[imindom,jmindom,nctstep,1],[imt,jmt,1,1],'st')
     hs(imt+1,:,1)     = hs(1,:,1)

   END IF

   ! Calculate SSH/depth
   WHERE (SUM(dzt(:,:,:,2),3) /= 0)
        zstot(1:imt,1:jmt,1)  = hs(1:imt,1:jmt,1)/SUM(dzt(:,:,:,2),3) + 1
   ELSEWHERE
        zstot(:,:,1) = 0.d0
   END WHERE

   WHERE (SUM(dzu(:,:,:,2),3) /= 0)
        zstou(1:imt,1:jmt) = 0.5*(hs(1:imt,1:jmt,0)+hs(2:imt+1,1:jmt,0))/SUM(dzu(:,:,:,2),3) + 1
   ELSEWHERE
        zstou = 0.d0
   END WHERE


   WHERE (SUM(dzv(:,1:jmt-1,:,2),3) /= 0)
        zstov(1:imt,1:jmt-1) = 0.5*(hs(1:imt,1:jmt-1,0)+hs(1:imt,2:jmt,0))/SUM(dzv(:,1:jmt-1,:,2),3) + 1
   ELSEWHERE
        zstov = 0.d0
   END WHERE

   ALLOCATE(tmptracer(imt,jmt,km),tmpuflux(imt,0:jmt,km), tmpvflux(imt,0:jmt,km))
   !ALLOCATE(tmp3d(imt,0:jmt,km), tmp2d(imt,0:jmt))

   ! Reading 2-time step variables
   ! In this case: velocities and tracers
   ! ===========================================================================
   ALLOCATE(tmp3d(imt,jmt,km))

   nctstep = currMon
   IF (l_onestep) nctstep = 1

   dateprefix = filledFileName(dateFormat, currYear, currMon, currDay)
   !write(*,*) 'currMon:',currMon
   !write(*,*) 'dateprefix:',dateprefix

   fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(dateprefix)//TRIM(fileSuffix)
   !write(*,*) 'reading fieldFile:', fieldFile

   !write(*,*) 'imindom:',imindom
   !write(*,*) 'jmindom:',jmindom
   !uvel(1:imt,1:jmt,km:1:-1) = get3DfieldNC(fieldFile, ueul_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')
   !write(*,*) 'line 71, read_field'
   uvel(1:imt,1:jmt,1:km) = get3DfieldNC(fieldFile, ueul_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')
   !write(*,*) 'line 73, read_field'

   IF (usgs_name/='') THEN
     !tmp3d(1:imt,1:jmt,km:1:-1)  = get3DfieldNC(fieldFile, usgs_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')
     tmp3d(1:imt,1:jmt,1:km)  = get3DfieldNC(fieldFile, usgs_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')
     uvel(1:imt,1:jmt,1:km) = uvel(1:imt,1:jmt,1:km) + tmp3d(1:imt,1:jmt,1:km)
   END IF
   uvel(1:imt,1:jmt,1:km) = uvel(1:imt,1:jmt,1:km)/1000.    !kg/s -> m3/s

   !vvel(1:imt,1:jmt,km:1:-1) = get3DfieldNC(fieldFile, veul_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')
   vvel(1:imt,1:jmt,1:km) = get3DfieldNC(fieldFile, veul_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')

   IF (vsgs_name/='') THEN
     !tmp3d(1:imt,1:jmt,km:1:-1) = get3DfieldNC(fieldFile, vsgs_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')
     tmp3d(1:imt,1:jmt,1:km) = get3DfieldNC(fieldFile, vsgs_name,[imindom,jmindom,1,nctstep],[imt,jmt,km,1],'st')
     vvel(1:imt,1:jmt,1:km) = vvel(1:imt,1:jmt,1:km) + tmp3d(1:imt,1:jmt,1:km)
   END IF
   vvel(1:imt,1:jmt,1:km) = vvel(1:imt,1:jmt,1:km)/1000.    !kg/s -> m3/s

    ! dz
    dzt(:,:,:,1) = get3DfieldNC(fieldFile, dzt_name,[imindom,jmindom,1,1],[imt,jmt,km,1],'st')
    !where (dzt < 100 ) dzt = 100.
    !where (dzt > 9.0E+036) dzt = 0.
    dzt(:,:,:,1) = dzt(:,:,:,1)/1.e4


   !! Tracers
   IF (l_tracers) THEN

     !write(*,*) 'numtracers:',numtracers
     DO itrac = 1, numtracers

        ! Make sure the data array is empty
        tmp3d(:,:,:) = 0.d0

        IF (tracers(itrac)%action == 'read') THEN

            ! Read the tracer from a netcdf file
            !fieldFile = TRIM(physDataDir)//TRIM(physPrefixForm)//TRIM(dateprefix)//TRIM(fileSuffix)
            !write(*,*) 'tracer name: ',trim(tracers(itrac)%varname)
            IF (tracers(itrac)%dimension == '3D') THEN
                tmp3d(1:imt,1:jmt,km:1:-1) = get3DfieldNC(fieldFile, tracers(itrac)%varname,[imindom,jmindom,1,nctstep] &
                          ,[imt,jmt,km,1],'st')
            ELSE IF (tracers(itrac)%dimension == '2D') THEN
                tmp3d(1:imt,1:jmt,1) = get2DfieldNC(fieldFile, tracers(itrac)%varname,[imindom,jmindom,nctstep,1] &
                                        ,[imt,jmt,1,1],'st')
            END IF
            !write(*,*) 'minval(tracer):',minval(tmp3d)
            !write(*,*) 'maxval(tracer):',maxval(tmp3d)

        ELSE IF (tracers(itrac)%action == 'compute') THEN

            ! Compute the tracer from a function defined in mod_tracer.F90
            CALL compute_tracer(tracers(itrac)%name, tmp3d(1:imt,1:jmt,1:km))

        ELSE
            PRINT '(A34,I4)', 'No action defined for this tracer:', itrac
            STOP 10
        END IF

        ! Store the information
        IF (tracers(itrac)%dimension == '3D') THEN
          !write(*,*) 'tracers(itrac)%scale:',tracers(itrac)%scale
          !write(*,*) 'tracers(itrac)%shift:',tracers(itrac)%shift
          !write(*,*) 'shape(tmp3d):',shape(tmp3d)
          !write(*,*) 'shape(tracers(itrac)%data):',shape(tracers(itrac)%data)
          !write(*,*) 'jmt:',jmt
          !write(*,*) 'itrac:',itrac

          !if (.not. allocated(tracers(itrac)%data)) then
            ! allocate the array
            !write(*,*) 'tracers(itrac) is not accocated, allocating it now'
            !allocate(tracers(itrac)%data(360, 385, 53, 2))
          !end if

          tmptracer = tracers(itrac)%scale*tmp3d(:,:,:) + tracers(itrac)%shift
          tracers(itrac)%data(:,1:jmt,:,2) = tmptracer
          !tracers(itrac)%data(:,1:jmt,:,2) = tracers(itrac)%scale*tmp3d(:,:,:) + tracers(itrac)%shift
          !tmptracer(:,1:jmt,:) = tracers(itrac)%scale*tmp3d(:,:,:) + tracers(itrac)%shift
        ELSE IF (tracers(itrac)%dimension == '2D') THEN
          tracers(itrac)%data(:,1:jmt,1,2) = tracers(itrac)%scale*tmp3d(:,:,1) + tracers(itrac)%shift
        END IF

     END DO
   END IF
   deallocate(tmp3d)
   deallocate(tmptracer)

   ! ===========================================================================

   ! uflux and vflux computation
   FORALL (kk = 1:km) uflux(:,:,kk,2)     = uvel(:,:,kk)
   FORALL (kk = 1:km) vflux(:,1:jmt,kk,2) = vvel(:,1:jmt,kk)

   ! dzdt calculation
!  IF (ints == 0 .AND. ( loopYears .EQV..FALSE.)) THEN
!     FORALL (kk = 1:km) dzdt(:,:,kk,2) = dzt(:,:,kk,2)*(zstot(:,:,1) - zstot(:,:,0))/tseas
!  ELSE IF (ints == intrun-1 .AND. ( loopYears .EQV..FALSE.)) THEN
!     FORALL (kk = 1:km) dzdt(:,:,kk,2) = dzt(:,:,kk,2)*(zstot(:,:,0) - zstot(:,:,-1))/tseas
!  ELSE
!     FORALL (kk = 1:km) dzdt(:,:,kk,2) = 0.5*dzt(:,:,kk,2)*(zstot(:,:,1) - zstot(:,:,-1))/tseas
!  END IF

   !! Zero meridional flux at j=0 and j=jmt
   vflux(:,0  ,:,:) = 0.d0
   vflux(:,jmt,:,:) = 0.d0

   ! Reverse the sign of fluxes if trajectories are run backward in time.
   CALL swap_sign()

END SUBROUTINE read_field
