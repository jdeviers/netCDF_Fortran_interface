PROGRAM openNetCDF
  USE netcdf
  implicit none

  CHARACTER(LEN=*),PARAMETER   :: INFILE = '../new.nc'
  INTEGER                      :: F_NCID
  INTEGER                      :: i,j,k

  INTEGER                      :: &
          S1_dim1_ID,S1_dim2_ID,  &
          S2_dim1_ID,S2_dim2_ID,  &
          S_xyz_ID, RI_ID

  INTEGER ::                       &
          S1_dim1_len,S1_dim2_len, &
          S2_dim1_len,S2_dim2_len, &
          S_xyz_len,RI_len

  INTEGER ::         &
          Sxyz1T_ID, &
          Sxyz2T_ID

  DOUBLE PRECISION,ALLOCATABLE :: &
          Sxyz1T(:,:,:,:),        &
          Sxyz2T(:,:,:,:)

  INTEGER ::             &
          Sxyz1T_str(4), &
          Sxyz2T_str(4)

  DOUBLE PRECISION :: tmp_Re,tmp_Im
  COMPLEX(8),ALLOCATABLE :: S1T(:,:,:),S2T(:,:,:)


! .. Open netCDF file
  CALL handle_err(NF90_OPEN(path=INFILE,mode=NF90_NOWRITE,ncid=F_NCID))

! .. Store dimension informations
  CALL handle_err(NF90_INQ_DIMID(ncid=F_NCID,name='S1_axis_1'     ,dimid=S1_dim1_ID))
  CALL handle_err(NF90_INQ_DIMID(ncid=F_NCID,name='S1_axis_2'     ,dimid=S1_dim2_ID))
  CALL handle_err(NF90_INQ_DIMID(ncid=F_NCID,name='S2_axis_1'     ,dimid=S2_dim1_ID))
  CALL handle_err(NF90_INQ_DIMID(ncid=F_NCID,name='S2_axis_2'     ,dimid=S2_dim2_ID))
  CALL handle_err(NF90_INQ_DIMID(ncid=F_NCID,name='S_xyz'         ,dimid=S_xyz_ID))
  CALL handle_err(NF90_INQ_DIMID(ncid=F_NCID,name='Real_Imaginary',dimid=RI_ID))

!  WRITE(*,*) 'Dimension ID of ','S1_axis_1',': ',S1_dim1_ID
!  WRITE(*,*) 'Dimension ID of ','S1_axis_2',': ',S1_dim2_ID
!  WRITE(*,*) 'Dimension ID of ','S2_axis_1',': ',S2_dim1_ID
!  WRITE(*,*) 'Dimension ID of ','S2_axis_2',': ',S2_dim2_ID
!  WRITE(*,*) 'Dimension ID of ','S1_xyz'   ,': ',S_xyz_ID
!  WRITE(*,*) 'Dimension ID of ','Real_Imag',': ',RI_ID

  CALL handle_err(NF90_INQUIRE_DIMENSION(ncid=F_NCID,dimid=S1_dim1_ID,len=S1_dim1_len))
  CALL handle_err(NF90_INQUIRE_DIMENSION(ncid=F_NCID,dimid=S1_dim2_ID,len=S1_dim2_len))
  CALL handle_err(NF90_INQUIRE_DIMENSION(ncid=F_NCID,dimid=S2_dim1_ID,len=S2_dim1_len))
  CALL handle_err(NF90_INQUIRE_DIMENSION(ncid=F_NCID,dimid=S2_dim2_ID,len=S2_dim2_len))
  CALL handle_err(NF90_INQUIRE_DIMENSION(ncid=F_NCID,dimid=S_xyz_ID  ,len=S_xyz_len))
  CALL handle_err(NF90_INQUIRE_DIMENSION(ncid=F_NCID,dimid=RI_ID     ,len=RI_len))



  ! .. Read in the entire Sxyz1T tensor, where the complex components are split into 2 reals
  CALL handle_err(NF90_INQ_VARID(ncid=F_NCID,name='Sxyz1T',varid=Sxyz1T_ID))
!  WRITE(*,*) 'Variable ID of ','Sxyz1T',': ',Sxyz1T_ID

  ALLOCATE( Sxyz1T(RI_len,S_xyz_len,S1_dim2_len,S1_dim1_len) ) ! Order of the dims must be REVERSED
!  PRINT*, SHAPE(Sxyz1T)
  Sxyz1T_str = (/1,1,1,1/)
  CALL handle_err(NF90_GET_VAR(ncid=F_NCID,varid=Sxyz1T_ID,values=Sxyz1T,start=Sxyz1T_str))

  WRITE(*,'(A,2E12.2)') 'First item in Sxyz1T(real and complex part): ',Sxyz1T(1:2,1,1,1)
  DEALLOCATE(Sxyz1T)

  ! .. Parsing the file and stepwise filling complex-valued S1T, This method is better.
  CALL handle_err(NF90_INQ_VARID(ncid=F_NCID,name='Sxyz1T',varid=Sxyz1T_ID))
  ALLOCATE( S1T(S1_dim1_len,S1_dim2_len,S_xyz_len) )
  DO i=1,S1_dim1_len
    DO j=1,S1_dim2_len
      DO k=1,S_xyz_len

        Sxyz1T_str = (/ 1,k,j,i /)
        CALL handle_err(NF90_GET_VAR(ncid=F_NCID,varid=Sxyz1T_ID,values=tmp_Re,start=Sxyz1T_str))

        Sxyz1T_str = (/ 2,k,j,i /)
        CALL handle_err(NF90_GET_VAR(ncid=F_NCID,varid=Sxyz1T_ID,values=tmp_Im,start=Sxyz1T_str))
        
        S1T(i,j,k) = CMPLX(tmp_Re,tmp_Im,KIND=8)
      END DO
    END DO
  END DO
  WRITE(*,'(A,(E12.2,SP,E12.2,"i"))') 'First item in S1T: ',S1T(1,1,1)
  DEALLOCATE(S1T)



  ! .. Read in the entire Sxyz2T tensor, where the complex components are split into 2 reals
  CALL handle_err(NF90_INQ_VARID(ncid=F_NCID,name='Sxyz2T',varid=Sxyz2T_ID))
!  WRITE(*,*) 'Variable ID of ','Sxyz2T',': ',Sxyz2T_ID

  ALLOCATE( Sxyz2T(RI_len,S_xyz_len,S2_dim2_len,S2_dim1_len) ) ! Order of the dims must be REVERSED, or use clever MAP
!  PRINT*, SHAPE(Sxyz2T)
  Sxyz2T_str = (/1,1,1,1/)
  CALL handle_err(NF90_GET_VAR(ncid=F_NCID,varid=Sxyz2T_ID,values=Sxyz2T,start=Sxyz2T_str))

  WRITE(*,'(A,2E12.2)') 'First item in Sxyz2T (real and complex part): ',Sxyz2T(1:2,1,1,1)
  DEALLOCATE(Sxyz2T)

  ! .. Parsing the file and stepwise filling complex-valued S2T. This method is better.
  CALL handle_err(NF90_INQ_VARID(ncid=F_NCID,name='Sxyz2T',varid=Sxyz2T_ID))
  ALLOCATE( S2T(S2_dim1_len,S2_dim2_len,S_xyz_len) )
  DO i=1,S2_dim1_len
    DO j=1,S2_dim2_len
      DO k=1,S_xyz_len

        Sxyz2T_str = (/ 1,k,j,i /)
        CALL handle_err(NF90_GET_VAR(ncid=F_NCID,varid=Sxyz2T_ID,values=tmp_Re,start=Sxyz2T_str))

        Sxyz2T_str = (/ 2,k,j,i /)
        CALL handle_err(NF90_GET_VAR(ncid=F_NCID,varid=Sxyz2T_ID,values=tmp_Im,start=Sxyz2T_str))
        
        S2T(i,j,k) = CMPLX(tmp_Re,tmp_Im,KIND=8)
      END DO
    END DO
  END DO
  WRITE(*,'(A,(E12.2,SP,E12.2,"i"))') 'First item in S2T: ',S2T(1,1,1)
  DEALLOCATE(S2T)

  CALL handle_err(NF90_CLOSE(ncid=F_NCID))


  
  contains
    SUBROUTINE handle_err(stat)
      INTEGER,INTENT(IN) :: stat

      IF (stat.NE.NF90_NOERR) THEN
        PRINT*, TRIM(NF90_STRERROR(stat))
        STOP "*** STOPPED ***"
      END IF

    END SUBROUTINE handle_err

END PROGRAM openNetCDF
