PROGRAM prmtop
!
! -- Source: https://ambermd.org/prmtop.pdf
!
  USE mod_procedures
  USE mod_interfaces
  USE mod_readprmtop
  USE mod_readnetcdf
  implicit none

  CHARACTER(LEN=*),PARAMETER :: &
    intop  = '/home/jean/Documents/SUP_binding/SUPinRandomLocations/Daniel_test/trajectories/noWAT.6pu0_Cl.prmtop'

  CHARACTER(LEN=80)          :: KEYWORD

  LOGICAL                    :: M1=.FALSE.,M2=.FALSE.
  INTEGER,ALLOCATABLE        :: MASK1(:),MASK2(:)
  INTEGER                    :: i,j,k,FRAME

  ! ----------------------------- PRMTOP INFO ---------------------------------
  INTEGER :: NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,NHPARM,NPARM
  INTEGER :: NNB,NRES,NBONA,NTHETA,NPHIA,NUMBND,NUMANG,NPTRA,NATYP,NPHB
  INTEGER :: IFPERT,NBPER,NGPER,NDPER,MBPER,MGPER,MDPER,IFBOX,NMXRS,IFCAP
  INTEGER :: IPTRES,NSPM,NSPSOL,NATCAP,IPOL
  REAL    :: OLDBETA,BOX(3),CUTCAP,XCAP,YCAP,ZCAP
  CHARACTER(LEN=80) :: RADIUS_SET


  CHARACTER(LEN=4),ALLOCATABLE :: &
    ATOM_NAME(:),RESIDUE_LABEL(:),RESIDUES(:),AMBER_ATOM_TYPE(:),             &
    TREE_CHAIN_CLASSIFICATION(:)

  REAL,            ALLOCATABLE :: &
    CHARGE(:),MASS(:),BOND_FORCE_CONSTANT(:),BOND_EQUIL_VALUE(:),             &
    ANGLE_FORCE_CONSTANT(:),ANGLE_EQUIL_VALUE(:),DIHEDRAL_FORCE_CONSTANT(:),  &
    DIHEDRAL_PERIODICITY(:),DIHEDRAL_PHASE(:),SCEE_SCALE_FACTOR(:),           &
    SNCB_SCALE_FACTOR(:),SOLTY(:),LENNARD_JONES_ACOEF(:),                     &
    LENNARD_JONES_BCOEF(:),HBOND_ACOEF(:),HBOND_BCOEF(:),HBCUT(:),RADII(:),   &
    POLARIZABILITY(:)

  INTEGER,         ALLOCATABLE :: &
    ATOMIC_NUMBER(:),ATOM_TYPE_INDEX(:),NUMBER_EXCLUDED_ATOMS(:),             &
    NONBONDED_PARM_INDEX(:),RESIDUE_POINTER(:),RESID(:),BONDS_INC_HYDROGEN(:),&
    BONDS_WITHOUT_HYDROGEN(:),ANGLES_INC_HYDROGEN(:),                         &
    ANGLES_WITHOUT_HYDROGEN(:),DIHEDRALS_INC_HYDROGEN(:),                     &
    DIHEDRALS_WITHOUT_HYDROGEN(:),EXCLUDED_ATOMS_LIST(:),JOIN_ARRAY(:),       &
    IROTAT(:),ATOMS_PER_MOLECULE(:)
  ! ---------------------------------------------------------------------------


  ! ----------------------------- NETCDF INFO ---------------------------------
  CHARACTER(LEN=*),PARAMETER :: &
    intraj = '/home/jean/Documents/SUP_binding/SUPinRandomLocations/Daniel_test/trajectories/6pu0_Cl_1.noWAT.nc'

  INTEGER          :: F_NCID
  INTEGER          :: NBFRAMES
  REAL,ALLOCATABLE :: COORDS1(:,:),COORDS2(:,:)
  REAL,ALLOCATABLE :: CHG1(:),CHG2(:)
  !REAL :: CLMB12

  ! ---------------------------------------------------------------------------


    ! READ THE AMBER PRMTOP
    CALL READTOP(intop,NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,NHPARM,NPARM, &
    NNB,NRES,NBONA,NTHETA,NPHIA,NUMBND,NUMANG,NPTRA,NATYP,NPHB,IFPERT,NBPER,            &
    NGPER,NDPER,MBPER,MGPER,MDPER,IFBOX,NMXRS,IFCAP,IPTRES,                             &
    NSPM,NSPSOL,NATCAP,IPOL,OLDBETA,BOX,CUTCAP,XCAP,YCAP,ZCAP,RADIUS_SET,               &
    ATOM_NAME,RESIDUE_LABEL,RESIDUES,RESID,AMBER_ATOM_TYPE,TREE_CHAIN_CLASSIFICATION,   &                     
    CHARGE,MASS,BOND_FORCE_CONSTANT,BOND_EQUIL_VALUE,ANGLE_FORCE_CONSTANT,              &
    ANGLE_EQUIL_VALUE,DIHEDRAL_FORCE_CONSTANT,DIHEDRAL_PERIODICITY,                     &
    DIHEDRAL_PHASE,SCEE_SCALE_FACTOR,SNCB_SCALE_FACTOR,SOLTY,LENNARD_JONES_ACOEF,       &
    LENNARD_JONES_BCOEF,HBOND_ACOEF,HBOND_BCOEF,HBCUT,RADII,POLARIZABILITY,             &                     
    ATOMIC_NUMBER,ATOM_TYPE_INDEX,NUMBER_EXCLUDED_ATOMS,                                &
    NONBONDED_PARM_INDEX,RESIDUE_POINTER,BONDS_INC_HYDROGEN,                            &
    BONDS_WITHOUT_HYDROGEN,ANGLES_INC_HYDROGEN,                                         &
    ANGLES_WITHOUT_HYDROGEN,DIHEDRALS_INC_HYDROGEN,                                     &
    DIHEDRALS_WITHOUT_HYDROGEN,EXCLUDED_ATOMS_LIST,JOIN_ARRAY,                          &
    IROTAT,ATOMS_PER_MOLECULE)


    ! OPEN NETCDF TRAJ
    CALL OPEN_NC(intraj,F_NCID)

    ! Check the prmtop and the netcdf have the same NATOM
    CALL check_NATOM(F_NCID,NATOM)

    CALL READNBFRAMES(F_NCID,NBFRAMES)
    WRITE(*,*) NBFRAMES

    ! Check whether files containing atom masks are provided.
    SELECT CASE (IARGC())
      CASE (1) 
        M1=.TRUE.
      CASE (2)
        M1=.TRUE.
        M2=.TRUE.
    END SELECT


    IF (M1) THEN

      CALL GETARG(1,KEYWORD) ! The first keyword is always mask1
      CALL READMASKFILE('mask1.dat',MASK1)
      ALLOCATE(COORDS1(UBOUND(MASK1,1),3))

      ! Add more vectors/tensors with properties of MASKX atoms
      ALLOCATE(CHG1(UBOUND(MASK1,1)))

      IF (M2) THEN

        CALL GETARG(2,KEYWORD) ! The second keyword is always mask2
        CALL READMASKFILE('mask2.dat',MASK2)
        ALLOCATE(COORDS2(UBOUND(MASK2,1),3))

        ! Add more vectors/tensors with properties of MASKX atoms
        ALLOCATE(CHG2(UBOUND(MASK2,1)))

      END IF

    ELSE

      ALLOCATE(MASK1(NATOM)) ! IF NO KEYWORD: the mask is all atoms
      DO i=1,NATOM
        MASK1(i) = i
      END DO
      ALLOCATE(COORDS1(UBOUND(MASK1,1),3))

    END IF



    WRITE(*,*) SHAPE(ATOM_NAME),SHAPE(RESIDUES),SHAPE(COORDS1)
    !CALL PRINTVEC(MASK1,MESSAGE='MASK 1:')
    !DO FRAME=1,NBFRAMES
    DO FRAME=1,5
      CALL READROWS(F_NCID,MASK1,COORDS1,FRAME)
      IF (M2) CALL READROWS(F_NCID,MASK2,COORDS2,FRAME)

      WRITE(*,'(/,A,I0,A)') 'FRAME ',FRAME,':'

      WRITE(*,'(/,A)') 'MASK 1:'
      DO j=1,UBOUND(MASK1,1)
        k=MASK1(j)
        WRITE(*,'(A,A,I0,4F12.6)') ATOM_NAME(k),RESIDUES(k),RESID(k),COORDS1(j,:),CHARGE(k)/18.2223
        CHG1(j) = CHARGE(k)/18.2223
      END DO

      CALL PRINTVEC(CHG1,'MASK1 ATOMS CHARGES:')
      WRITE(*,'(A,3F8.3,/)') 'COM of MASK1 atoms = ',CALC_COM(COORDS1)

      IF (M2) THEN
        WRITE(*,'(/,A)') 'MASK 2:'
        DO j=1,UBOUND(MASK2,1)
          k=MASK2(j)
          WRITE(*,'(A,A,I0,4F12.6)') ATOM_NAME(k),RESIDUES(k),RESID(k),COORDS1(j,:),CHARGE(k)/18.2223
          CHG2(j) = CHARGE(k)/18.2223
        END DO

        CALL PRINTVEC(CHG2,'MASK2 ATOMS CHARGES:')
        WRITE(*,'(A,3F8.3,/)') 'COM of MASK2 atoms = ',CALC_COM(COORDS2)
      END IF


      IF (M2) THEN
        ! Coulombic interaction between the 2 fragments:
        WRITE(*,'(A,E12.3,/)') 'Coulomb interaction b/w frags 1 and 2 = ',F1F2COULOMB(COORDS1,CHG1,COORDS2,CHG2)
      END IF

    END DO


    !CALL PRINTMAT(COORDS,'LAST COORDS MAT: ')
    DEALLOCATE(COORDS1)
    IF (ALLOCATED(COORDS2)) DEALLOCATE(COORDS2)


    ! CLOSE NETCDF TRAJ
    CALL CLOSE_NC(F_NCID)

END PROGRAM prmtop