MODULE mod_readnetcdf
  USE mod_procedures
  USE mod_interfaces
  USE netcdf
  implicit none

  contains
!
  SUBROUTINE handle_err(stat)
    INTEGER,INTENT(IN) :: stat

    IF (stat.NE.NF90_NOERR) THEN
      PRINT*, TRIM(NF90_STRERROR(stat))
      STOP "*** STOPPED ***"
    END IF

  END SUBROUTINE handle_err
!
  SUBROUTINE OPEN_NC(infile,F_NCID)
    implicit none

    ! .. Parameters ..
    CHARACTER(LEN=*),INTENT(IN)  :: infile
    INTEGER,         INTENT(OUT) :: F_NCID

    CALL handle_err(NF90_OPEN(path=infile,mode=NF90_NOWRITE,ncid=F_NCID))

  END SUBROUTINE OPEN_NC
!
  SUBROUTINE CLOSE_NC(F_NCID)
    implicit none

    ! .. Parameters ..
    INTEGER,INTENT(IN) :: F_NCID

    CALL handle_err(NF90_CLOSE(ncid=F_NCID))

  END SUBROUTINE CLOSE_NC
!
  SUBROUTINE READNBFRAMES(F_NCID,NBFRAMES)
    implicit none

    ! .. Parameters ..
    INTEGER,INTENT(IN)  :: F_NCID
    INTEGER,INTENT(OUT) :: NBFRAMES
    ! .. NetCDF vars ..
    INTEGER :: frames_ID

    ! .. Find the ID of the frames variable 
    CALL handle_err(NF90_INQ_DIMID(ncid=F_NCID,name='frame',dimid=frames_ID))

    ! .. Find the dimension (size) of 'dimension' frame (through its ID)
    CALL handle_err(NF90_INQUIRE_DIMENSION(ncid=F_NCID,dimid=frames_ID,len=NBFRAMES))

  END SUBROUTINE READNBFRAMES
!
  SUBROUTINE check_NATOM(F_NCID,NATOM_top)
    implicit none

    ! .. Parameters ..
    INTEGER,INTENT(IN)  :: F_NCID,NATOM_top
    ! .. NetCDF vars ..
    INTEGER :: NATOM_ID,NATOM_nc

    CALL handle_err(NF90_INQ_DIMID(ncid=F_NCID,name='atom',dimid=NATOM_ID))
    CALL handle_err(NF90_INQUIRE_DIMENSION(ncid=F_NCID,dimid=NATOM_ID,len=NATOM_nc))

    IF (NATOM_top.NE.NATOM_nc) ERROR STOP "*** TOP AND NC HAVE DIFF NATOM: FATAL ERROR ***"

  END SUBROUTINE check_NATOM
!
  SUBROUTINE READCOORDSATFRAME(F_NCID,frame,COORDS) ! THIS FAILS: x,y,z are read vertically instead of horizontally
    implicit none

    ! .. Parameters ..
    INTEGER,INTENT(IN)             :: F_NCID,frame
    REAL,ALLOCATABLE,INTENT(INOUT) :: COORDS(:,:)
    ! .. NetCDF vars ..
    INTEGER :: coords_ID


    CALL handle_err(NF90_INQ_VARID(ncid=F_NCID,name='coordinates',varid=coords_ID))
    CALL handle_err(NF90_GET_VAR(ncid=F_NCID,varid=coords_ID,values=COORDS,start=(/1,1,frame/),count=(/3,8214,1/),map=(/3,1/) ))

  END SUBROUTINE READCOORDSATFRAME
!

  SUBROUTINE READROWS(F_NCID,MASK,COORDS,FRAME)
    implicit none

    INTEGER,            INTENT(IN)    :: F_NCID,FRAME
    INTEGER,ALLOCATABLE,INTENT(IN)    :: MASK(:)
    REAL,   ALLOCATABLE,INTENT(INOUT) :: COORDS(:,:)

    INTEGER :: coords_ID,i,ATINDEX

    CALL handle_err(NF90_INQ_VARID(ncid=F_NCID,name='coordinates',varid=coords_ID))

    DO i=1,UBOUND(MASK,1)
      ATINDEX=MASK(i)
      CALL handle_err(NF90_GET_VAR(ncid=F_NCID,varid=coords_ID,values=COORDS(i,:),start=(/1,ATINDEX,FRAME/),count=(/3,1,1/)))
    END DO
    !WRITE(*,*) COORDS ! This works; COORDS gets lost during export it seems.

  END SUBROUTINE READROWS

END MODULE mod_readnetcdf