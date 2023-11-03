SUBROUTINE SYST
  USE mar_data
  LOGICAL :: FOR_ALL,LIN
  INTEGER :: I
  LOGICAL :: ENTER_NAME
  CHARACTER*200 :: NAM

  IF (NDAT<1) THEN
     WRITE(6,*) "<E> No data set in memory!"
     RETURN
  END IF
  IF (NDAT>1) THEN
     FOR_ALL=LIN("Change system response for all data sets? ",.TRUE.)
  ELSE
     FOR_ALL=.FALSE.
  ENDIF
  SYS_RESP(ACT_DAT)=LIN("Use system response? ",SYS_RESP(ACT_DAT))
  IF (SYS_RESP(ACT_DAT)) THEN
     IF (ENTER_NAME(NAM,'sys')) THEN
!     WRITE(6,'(A,$)') ' Name of system response file (incl. extension): '
!     READ(5,*) SYS_NAME(ACT_DAT)
        SYS_NAME(ACT_DAT)=TRIM(ADJUSTL(NAM))
     ELSE
        RETURN
     END IF
  END IF
  IF (FOR_ALL) THEN
     DO I=1,NDAT
        SYS_RESP(I)=SYS_RESP(ACT_DAT)
        SYS_NAME(I)=SYS_NAME(ACT_DAT)
     ENDDO
  ENDIF
  CALL LOAD_SYS
END SUBROUTINE SYST

SUBROUTINE LOAD_SYS
  USE mar_data
  CHARACTER*200 :: NAM,DUM
  INTEGER UN,J1,I,J2
  CHARACTER :: DUM1

  CALL GET_UN(UN)
  WRITE(6,*) '<I> Reading system responses'
  IF (ALLOCATED(SYS_N)) DEALLOCATE(SYS_N,SYS_D,SYS_T,SYS_PER,SYS_TYPE)
  ALLOCATE(SYS_N(NDAT),SYS_PER(NDAT),SYS_TYPE(NDAT))
  J1=0
  DO I=1,NDAT
     IF (SYS_RESP(I)) THEN
        NAM=TRIM(ADJUSTL(SYS_NAME(I)))
        OPEN(UNIT=UN,FILE=NAM,STATUS='OLD',ERR=905)
        DO J2=1,4
           READ(UN,*) DUM
        END DO
        READ(UN,'(A1,I5)',ERR=905,END=905) DUM1,SYS_N(I)
        READ(UN,*) DUM
        READ(UN,'(A1,G14.6)',ERR=905,END=905) DUM1,SYS_PER(I)
        READ(UN,*) DUM
        READ(UN,'(A1,I5)',ERR=905,END=905) DUM1,SYS_TYPE(I)
        CLOSE(UN)
        J1=MAX(SYS_N(I),J1)
     ENDIF
  ENDDO
  ALLOCATE(SYS_D(NDAT,J1),SYS_T(NDAT,J1))
  DO I=1,NDAT
     IF (SYS_RESP(I)) THEN
        NAM=TRIM(ADJUSTL(SYS_NAME(I)))
        OPEN(UNIT=UN,FILE=NAM,STATUS='OLD',ERR=905)
        DO J2=1,11
           READ(UN,*) DUM
        END DO
        DO J2=1,SYS_N(I)
           READ(UN,*,ERR=905,END=905) SYS_T(I,J2),SYS_D(I,J2)
        END DO
        CLOSE(UN)
        IF (SYS_T(I,1).NE.0.0) THEN
           WRITE(6,*) "<I> Warning: System response does not start at t=0.0 s"
           WRITE(6,*) "<I> This will produce a shift in the data"
        END IF
     ENDIF
  ENDDO
  GOTO 999
905 WRITE(6,*) '<E> Error reading system response for data set ',I
  IF (ALLOCATED(SYS_T)) DEALLOCATE(SYS_T,SYS_D)
  IF (ALLOCATED(SYS_N)) DEALLOCATE(SYS_N,SYS_PER,SYS_TYPE)
999 CALL FREE_UN(UN)
END SUBROUTINE LOAD_SYS

SUBROUTINE INTERPOLATE_SYST
  USE mar_para
  LOGICAL :: LIN
  INTEGER :: I,J,UN
  LOGICAL :: ENTER_NAME,CPOL,POL
  CHARACTER*200 :: NAM,DUM
  CHARACTER (LEN=15) :: ISYS_NAME
  REAL (OP), DIMENSION(:), ALLOCATABLE :: ISYS_T,ISYS_D
  REAL (OP) :: ISYS_PER,IDT,Y,T
  INTEGER :: ISYS_N,ISYS_TYPE,IPERIODS,PPP
  REAL (OP), DIMENSION(:), ALLOCATABLE :: YAB,PK,QK,AS,BS,CS,DS
  CHARACTER :: DUM1

  WRITE(6,*) 'Loading Martin system response file:'
  IF (ENTER_NAME(NAM,'sys')) THEN
     ISYS_NAME=TRIM(ADJUSTL(NAM))
  ELSE
     RETURN
  END IF
  CALL GET_UN(UN)
  WRITE(6,*) '<I> Reading system responses'
  NAM=TRIM(ADJUSTL(ISYS_NAME))
  OPEN(UNIT=UN,FILE=NAM,STATUS='OLD',ERR=905)
  DO I=1,4
     READ(UN,*) DUM
  END DO
  READ(UN,'(A1,I5)',ERR=905,END=905) DUM1,ISYS_N
  READ(UN,*) DUM
  READ(UN,'(A1,G14.6)',ERR=905,END=905) DUM1,ISYS_PER
  READ(UN,*) DUM
  READ(UN,'(A1,I5)',ERR=905,END=905) DUM1,ISYS_TYPE
  REWIND(UN)
  ALLOCATE(ISYS_D(ISYS_N),ISYS_T(ISYS_N))
  DO I=1,11
     READ(UN,*) DUM
  END DO
  DO I=1,ISYS_N
     READ(UN,*,ERR=905,END=905) ISYS_T(I),ISYS_D(I)
  END DO
  CLOSE(UN)
  IF (ISYS_N.LE.2) THEN
     WRITE(6,*) "<E> Not enough data points in SYS-file"
     GOTO 999
  END IF
  IDT=ISYS_T(2)-ISYS_T(1)
  CALL RIN('Enter new sampling interval (s):',IDT)
  IPERIODS=0
  CALL IIN('Repeat sequence how often?:',IPERIODS)
  CPOL=LIN('Change polarity every other sequence?',.TRUE.)
  POL=.NOT.CPOL
  WRITE(6,'(A,$)') ' Enter file name for interpolated system response file: '
  READ(5,*) NAM
  ISYS_PER=ISYS_PER/1000.0_OP
  PPP=INT(ISYS_PER/IDT+0.5)
  WRITE(6,*) "<I> ",PPP," data points per period"
  ALLOCATE(YAB(ISYS_N),PK(ISYS_N-1),QK(ISYS_N-1),AS(ISYS_N-1),BS(ISYS_N-1),CS(ISYS_N-1),DS(ISYS_N-1))
  YAB(1)=0
  YAB(ISYS_N)=(ISYS_D(ISYS_N)-ISYS_D(ISYS_N-1))/(ISYS_T(ISYS_N)-ISYS_T(ISYS_N-1))
  PK=20.0_OP
  QK=20.0_OP
  CALL RASPL(ISYS_N,ISYS_T,ISYS_D,PK,QK,YAB,AS,BS,CS,DS)
  OPEN(UNIT=UN,FILE=NAM,STATUS='UNKNOWN',ERR=910)
  
  DO I=1,IPERIODS
     IF (CPOL) POL=.NOT.POL
     DO J=1,PPP
        T=(J-1)*IDT
        CALL YSPLINE(ISYS_N,ISYS_T,T,AS,BS,CS,DS,Y,PK,QK)
        IF (.NOT.POL) Y=-Y 
        WRITE(UN,'(2G14.6)',ERR=910) T+(I-1)*ISYS_PER,Y
     END DO
  END DO
  CLOSE(UN)
  GOTO 999
905 WRITE(6,*) '<E> Error reading system response file'
  GOTO 999
910 WRITE(6,*) '<E> Error writing interpolated system response file'
  
999 IF (ALLOCATED(PK)) DEALLOCATE(YAB,PK,QK,AS,BS,CS,DS)
  CALL FREE_UN(UN)
END SUBROUTINE INTERPOLATE_SYST
