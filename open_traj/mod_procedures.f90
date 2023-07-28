MODULE mod_procedures
  implicit none


  contains

  SUBROUTINE READVEC_R(VEC,VECSIZE,VECROW,CFMT)
    implicit none

    ! .. Params ..
    INTEGER,         INTENT(IN)    :: VECSIZE,VECROW
    REAL,            INTENT(INOUT) :: VEC(VECSIZE)
    CHARACTER(LEN=*),INTENT(IN)    :: CFMT
    ! .. Boundary variables ..
    INTEGER                        :: BCURSOR,ECURSOR

    BCURSOR=1 ; ECURSOR=BCURSOR+VECROW-1
    DO WHILE(.TRUE.)
      READ(10,FMT=TRIM(CFMT)) VEC(BCURSOR:ECURSOR)
      !WRITE(*,FMT=TRIM(CFMT)) VEC(BCURSOR:ECURSOR)
      BCURSOR=BCURSOR+VECROW
      ECURSOR=ECURSOR+VECROW
      IF (ECURSOR.GT.VECSIZE) ECURSOR=VECSIZE
      IF (BCURSOR.GT.ECURSOR) EXIT       
    END DO

  END SUBROUTINE READVEC_R

  SUBROUTINE READVEC_I(VEC,VECSIZE,VECROW,CFMT)
    implicit none

    ! .. Params ..
    INTEGER,         INTENT(IN)    :: VECSIZE,VECROW
    INTEGER,         INTENT(INOUT) :: VEC(VECSIZE)
    CHARACTER(LEN=*),INTENT(IN)    :: CFMT
    ! .. Boundary variables ..
    INTEGER                        :: BCURSOR,ECURSOR

    BCURSOR=1 ; ECURSOR=BCURSOR+VECROW-1
    DO WHILE(.TRUE.)
      READ(10,FMT=TRIM(CFMT)) VEC(BCURSOR:ECURSOR)
      !WRITE(*,FMT=TRIM(CFMT)) VEC(BCURSOR:ECURSOR)
      BCURSOR=BCURSOR+VECROW
      ECURSOR=ECURSOR+VECROW
      IF (ECURSOR.GT.VECSIZE) ECURSOR=VECSIZE
      IF (BCURSOR.GT.ECURSOR) EXIT       
    END DO

  END SUBROUTINE READVEC_I  

  SUBROUTINE READVEC_C(VEC,VECSIZE,VECROW,CFMT)
    implicit none

    ! .. Params ..
    INTEGER,         INTENT(IN)    :: VECSIZE,VECROW
    CHARACTER(LEN=*),INTENT(INOUT) :: VEC(VECSIZE)
    CHARACTER(LEN=*),INTENT(IN)    :: CFMT
    ! .. Boundary variables ..
    INTEGER                        :: BCURSOR,ECURSOR

    BCURSOR=1 ; ECURSOR=BCURSOR+VECROW-1
    DO WHILE(.TRUE.)
      READ(10,FMT=TRIM(CFMT)) VEC(BCURSOR:ECURSOR)
      !WRITE(*,FMT=TRIM(CFMT)) VEC(BCURSOR:ECURSOR)
      BCURSOR=BCURSOR+VECROW
      ECURSOR=ECURSOR+VECROW
      IF (ECURSOR.GT.VECSIZE) ECURSOR=VECSIZE
      IF (BCURSOR.GT.ECURSOR) EXIT       
    END DO

  END SUBROUTINE READVEC_C

  SUBROUTINE PRINTVEC_R(VEC,MESSAGE)
    implicit none

    REAL,                     INTENT(IN) :: VEC(:)
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: MESSAGE

    INTEGER :: i

    IF (PRESENT(MESSAGE)) WRITE(*,'(/,A)') MESSAGE
    DO i=1,UBOUND(VEC,1)
      WRITE(*,'(F16.8)') VEC(i)
    END DO

  END SUBROUTINE PRINTVEC_R

  SUBROUTINE PRINTVEC_I(VEC,MESSAGE)
    implicit none

    INTEGER,                  INTENT(IN) :: VEC(:)
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: MESSAGE

    INTEGER :: i

    IF (PRESENT(MESSAGE)) WRITE(*,'(/,A)') MESSAGE
    DO i=1,UBOUND(VEC,1)
      WRITE(*,'(I8)') VEC(i)
    END DO

  END SUBROUTINE PRINTVEC_I

  SUBROUTINE PRINTVEC_C(VEC,MESSAGE)
    implicit none

    CHARACTER(LEN=4),         INTENT(IN) :: VEC(:)
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: MESSAGE

    INTEGER :: i

    IF (PRESENT(MESSAGE)) WRITE(*,'(/,A)') MESSAGE
    DO i=1,UBOUND(VEC,1)
      WRITE(*,'(A)') VEC(i)
    END DO

  END SUBROUTINE PRINTVEC_C

  SUBROUTINE PRINTMAT_R(MAT,MESSAGE)
    implicit none

    REAL,                     INTENT(IN) :: MAT(:,:)
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: MESSAGE

    INTEGER :: i,j

    IF (PRESENT(MESSAGE)) WRITE(*,'(/,A)') MESSAGE
    DO i=1,UBOUND(MAT,1)
      WRITE(*,'(*(F16.8))') (MAT(i,j), j=1,UBOUND(MAT,2))
    END DO

  END SUBROUTINE PRINTMAT_R

  SUBROUTINE PRINTMAT_I(MAT,MESSAGE)
    implicit none

    INTEGER,                  INTENT(IN) :: MAT(:,:)
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: MESSAGE

    INTEGER :: i,j

    IF (PRESENT(MESSAGE)) WRITE(*,'(/,A)') MESSAGE
    DO i=1,UBOUND(MAT,1)
      WRITE(*,'(*(I8))') (MAT(i,j), j=1,UBOUND(MAT,2))
    END DO

  END SUBROUTINE PRINTMAT_I

  SUBROUTINE PRINTMAT_C(MAT,MESSAGE)
    implicit none

    CHARACTER(LEN=4),         INTENT(IN) :: MAT(:,:)
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: MESSAGE

    INTEGER :: i,j

    IF (PRESENT(MESSAGE)) WRITE(*,'(/,A)') MESSAGE
    DO i=1,UBOUND(MAT,1)
      WRITE(*,'(*(A4))') (MAT(i,j), j=1,UBOUND(MAT,2))
    END DO

  END SUBROUTINE PRINTMAT_C
!
  SUBROUTINE READMASKFILE(infile,MASK1)
    implicit none

    CHARACTER(LEN=*),   INTENT(IN)    :: infile
    INTEGER,ALLOCATABLE,INTENT(INOUT) :: MASK1(:)

    INTEGER :: CNT,i,io

    CNT=0
    OPEN(10,FILE=infile,STATUS='OLD',ACTION='READ',IOSTAT=io)
    DO WHILE (io.EQ.0)
      READ(10,*,IOSTAT=io)
      CNT=CNT+1
    END DO
    CNT=CNT-1
    ALLOCATE(MASK1(CNT))
    REWIND(10)
    DO i=1,CNT
      READ(10,*,IOSTAT=io) MASK1(i)
    END DO
    CLOSE(10)

  END SUBROUTINE READMASKFILE
!
  REAL FUNCTION DISTANCE(XYZ1,XYZ2)
    implicit none

    REAL,INTENT(IN) :: XYZ1(3),XYZ2(3)

    DISTANCE=SQRT( (XYZ2(1)-XYZ1(1))**2. + (XYZ2(2)-XYZ1(2))**2. + (XYZ2(3)-XYZ1(3))**2. )

  END FUNCTION DISTANCE
!
  REAL FUNCTION COULOMB(Q1,Q2,R12)
    implicit none

    REAL,INTENT(IN) :: Q1,Q2,R12

    COULOMB = (Q1*Q2) / (R12**2.)

  END FUNCTION COULOMB
!
  FUNCTION CALC_COM(XYZ)
    implicit none

    REAL,INTENT(IN)   :: XYZ(:,:)
    REAL,DIMENSION(3) :: CALC_COM
    REAL,DIMENSION(3) :: tmpCOM
    INTEGER           :: i

    tmpCOM=0.d0
    DO i=1,UBOUND(XYZ,1)
      tmpCOM(:) = tmpCOM(:) + XYZ(i,:)
    END DO
    CALC_COM(:) = tmpCOM(:)/UBOUND(XYZ,1)

  END FUNCTION CALC_COM
!
  FUNCTION WEIGHTED_COM(XYZ,WEIGHTS)
    ! Can do mass-weighting, or find the COM of charges
    implicit none

    REAL,INTENT(IN)   :: XYZ(:,:)
    REAL,INTENT(IN)   :: WEIGHTS(:)
    REAL,DIMENSION(3) :: WEIGHTED_COM
    REAL,DIMENSION(3) :: tmpCOM
    INTEGER           :: i

    IF (UBOUND(XYZ,1).NE.UBOUND(WEIGHTS,1)) ERROR STOP "XYZ and WEIGHTS have different length"

    tmpCOM=0.d0
    DO i=1,UBOUND(XYZ,1)
      tmpCOM(:) = tmpCOM(:) + ( XYZ(i,:)*WEIGHTS(i) )
    END DO
    WEIGHTED_COM(:) = tmpCOM(:)/UBOUND(XYZ,1)

  END FUNCTION WEIGHTED_COM
!
  REAL FUNCTION F1F2COULOMB(XYZ1,CHG1,XYZ2,CHG2)
    implicit none

    REAL,INTENT(IN) :: XYZ1(:,:),XYZ2(:,:)
    REAL,INTENT(IN) :: CHG1(:),CHG2(:)
    REAL            :: tmpCoulomb
    INTEGER         :: i,j

    tmpCoulomb=0.d0
    DO i=1,UBOUND(XYZ1,1)
      DO j=1,UBOUND(XYZ2,1)
        tmpCoulomb = tmpCoulomb + COULOMB( CHG1(i),CHG2(j),DISTANCE( XYZ1(i,:),XYZ2(j,:) ) )
        !PRINT*,tmpCoulomb
      END DO
    END DO
    F1F2COULOMB = tmpCoulomb

  END FUNCTION F1F2COULOMB
!
  END MODULE mod_procedures