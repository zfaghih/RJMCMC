LOGICAL FUNCTION ENTER_NAME(FILE_NAME,EXT)
  CHARACTER(LEN=200) FILE_NAME,FILNAM
  CHARACTER(LEN=3) EXT
  INTEGER I
  LOGICAL LBLANK,EX

  WRITE(6,*)'Enter * to see all available '//EXT//'-files'
1 WRITE(6,'(A,$)') ' Enter file name (without extension): '
  READ(5,'(A)') FILNAM
! Check whether directory requested and if so, do it:
  I=INDEX(FILNAM,'*')
! directory also requested when blank file name
  LBLANK=FILNAM(1:4).EQ.'    '
  IF(I/=0.OR.LBLANK) THEN              ! directory requested
     CALL SYSTEM('ls -k *.'//EXT)      ! file directory
     GOTO1
  ENDIF
  I=INDEX(FILNAM,'.')
  IF (I==0) THEN
     FILE_NAME=TRIM(ADJUSTL(FILNAM))//'.'//EXT
  ELSE
     FILE_NAME=TRIM(ADJUSTL(FILNAM))
  ENDIF
  INQUIRE(FILE=FILE_NAME,EXIST=EX)
  ENTER_NAME=.TRUE.
  IF (.NOT.EX) THEN
     WRITE(6,*) '<E> File does not exist'
  ENDIF
END FUNCTION ENTER_NAME

SUBROUTINE FREE_UN(N)
  USE mar_data
  INTEGER N
  UN_NO(N)=.FALSE.
END SUBROUTINE FREE_UN

SUBROUTINE GET_UN(N) 
  USE mar_data
  INTEGER I,J,N
  J=0
  DO I=100,20,-1
     IF (.NOT.UN_NO(I).AND.J==0) THEN
        J=I
     END IF
  END DO
  N=J
  UN_NO(J)=.TRUE.
END SUBROUTINE GET_UN

SUBROUTINE SAVE(WHAT)
  USE mar_data
  CHARACTER*200 NAM
  CHARACTER*4 FTYPE
  INTEGER UN,I,J,WHAT,OUT,I1,I2,NOFF
  LOGICAL OUT_MOD,LIN
  REAL (OP) DEPTH

  UN=-1
  IF (NLAY==0.AND.NDAT==0) THEN
     WRITE(6,*) "<E> No data in memory"
     RETURN
  ENDIF
  IF (WHAT==2.AND.NLAY==0) WHAT=3
  IF (NDAT==0.AND.WHAT==2) WHAT=1
  IF (WHAT==1.AND.NLAY==0) THEN
     WRITE(6,*) "<E> No model in memory"
     RETURN
  ENDIF
  IF (WHAT==3.AND.NDAT==0) THEN
     WRITE(6,*) "<E> No data sets in memory"
     RETURN
  ENDIF
     
  WRITE(6,'(A,$)') ' Enter file name: '
  READ(5,*) NAM
  I=INDEX(NAM,'.')
  IF (I==0) THEN
     NAM=TRIM(ADJUSTL(NAM))//'.mar'
  ELSE
     NAM=TRIM(ADJUSTL(NAM))
  ENDIF
  CALL GET_UN(UN)
  OUT=0
  IF (WHAT==3.AND.NDAT>1) THEN
     OUT=ACT_DAT
     OUT_MOD=LIN('Write out model? ',.FALSE.)
     IF (OUT_MOD) WHAT=2
  ENDIF
  IF (WHAT==1) THEN
     FTYPE='MODL'
  ELSEIF (WHAT==3) THEN
     FTYPE='DATA'
  ELSE
     FTYPE='CMPL'
  ENDIF
  OPEN(UNIT=UN,FILE=NAM,STATUS='UNKNOWN',ERR=900)
  WRITE(UN,'(A)',ERR=900) '#######################'
  WRITE(UN,'(A)',ERR=900) '# Martin-file'
  WRITE(UN,'(A)',ERR=900) '#'
  WRITE(UN,'(A)',ERR=900) '# File version:'
  WRITE(UN,'(A1,I5)',ERR=900) '#',FVER
  WRITE(UN,'(A)',ERR=900) '# File Type:'
  WRITE(UN,'(A2,A4)',ERR=900) '# ',FTYPE
  IF (WHAT<3) THEN
     WRITE(UN,'(A)',ERR=900) '#######################'
     WRITE(UN,'(A)',ERR=900) '# Number of layers'
     WRITE(UN,'(A1,I5)',ERR=900) '#',NLAY
     WRITE(UN,'(A)',ERR=900) '# Depth and resistivity of water'
     WRITE(UN,'(A1,2G14.6)',ERR=900) '#',W_DEPTH,W_RES
     WRITE(UN,'(A)',ERR=900) '# No, resistivity, anisotropy, thickness and depth of layer'
     DEPTH=0
     DO I=1,NLAY-1
        WRITE(UN,'(A1,I5,4G14.6)',ERR=900) '#',I,RES(I),ANIS(I),TNESS(I),DEPTH
        DEPTH=DEPTH+TNESS(I)
     END DO
     WRITE(UN,'(A1,I5,2G14.6,A14,G14.6)',ERR=900) '#',I,RES(NLAY),ANIS(NLAY),' ',DEPTH
     IF (NDAT>0.AND.WHAT==1) THEN
        WRITE(UN,'(A)',ERR=900) '#######################'
        WRITE(UN,'(A)',ERR=900) '# Calibration factors for data sets:'
        WRITE(UN,'(A1,I5)',ERR=900) '#',NDAT
        WRITE(UN,'(A)',ERR=900) '# Offset/Calibration-factors:'
        DO I=1,NDAT
           DEPTH=SQRT(X(I)*X(I)+Y(I)*Y(I)+Z(I)*Z(I))
           WRITE(UN,'(A2,2G15.6)',ERR=900)'# ',DEPTH,CALF(I)
        END DO
     ENDIF
     WRITE(6,*) '<I> Saved model'
  ENDIF
  IF (WHAT>1.AND.NDAT>0) THEN
     WRITE(UN,'(A)',ERR=900) '#######################'
     WRITE(UN,'(A)',ERR=900) '# Number of data sets'
     IF (OUT/=0) THEN
        I1=OUT
        I2=OUT
     ELSE
        I1=1
        I2=NDAT
     END IF
     WRITE(UN,'(A1,I5)',ERR=900) '#',I2-I1+1
     WRITE(UN,'(A)',ERR=900) '# Total misfit'
     WRITE(UN,'(A1,G14.6)',ERR=900) '#',TCHI
     WRITE(UN,'(A)',ERR=900) '# Tx-designation'
     WRITE(UN,'(A2,A15)',ERR=900) '# ',TX_DESIG
     WRITE(UN,'(A)',ERR=900) '# Tx-coordinates'
     WRITE(UN,'(A1,3G14.6)',ERR=900) '#',TX_X,TX_Y,TX_Z
     WRITE(UN,'(A)',ERR=900) '# Tx dipole length / current'
     WRITE(UN,'(A1,2G14.6)',ERR=900) '#',DL,CUR
     WRITE(UN,'(A)',ERR=900) '# Horizontal Tx angle in degree'
     WRITE(UN,'(A1,G14.6)',ERR=900) '#',0.0
     WRITE(UN,'(A)',ERR=900) '# Frequency domain data?'
     WRITE(UN,'(A2,L1)',ERR=900) "# ",LFD
     DO I=I1,I2
        WRITE(UN,'(A)',ERR=900) '#######################'
        WRITE(UN,'(A5,A15,I3)',ERR=900) '# Rx ',TX_DESIG,I
        WRITE(UN,'(A)',ERR=900) '# Type of data set'
        WRITE(UN,'(A2,A4)',ERR=900) '# ',TYPE(I)
        IF (.NOT.LFD) THEN
           WRITE(UN,'(A)',ERR=900) '# Impulse, step-off response, time derivative'
           WRITE(UN,'(A2,L1,A1,L1,A1,L1)',ERR=900) '# ',IMP(I),' ', &
                STEP_OFF(I),' ',DERIVATE(I)
        ENDIF
        WRITE(UN,'(A)',ERR=900) '# Number of data points and misfit'
        IF (LFD) THEN
           WRITE(UN,'(A1,I5,G14.6)',ERR=900) '#',NT(I)/2,CHI(I)
        ELSE
           WRITE(UN,'(A1,I5,G14.6)',ERR=900) '#',NT(I),CHI(I)
        ENDIF
        WRITE(UN,'(A)',ERR=900) '# Rx-offset'
        WRITE(UN,'(A1,3G14.6)',ERR=900) '#',X(I),Y(I),Z(I)
        WRITE(UN,'(A)',ERR=900) '# Sensor angles'
        WRITE(UN,'(A1,2G14.6)',ERR=900) '#',ANG1(I),ANG2(I)
        IF (.NOT.LFD) THEN
           WRITE(UN,'(A)',ERR=900) '# System-response'
           WRITE(UN,'(A2,L1,A1,A15)',ERR=900)'# ',SYS_RESP(I),' ',SYS_NAME(I)
           WRITE(UN,'(A)',ERR=900) '# Delay'
           WRITE(UN,'(A1,G14.6)',ERR=900) '#',DELAY(I)
        ENDIF
        WRITE(UN,'(A)',ERR=900) '# Calibration-factor'
        WRITE(UN,'(A2,L1,A1,G15.6)',ERR=900)'# ',LFCAL(I),' ',CALF(I)
        IF (LFD) THEN
           WRITE(UN,'(A)',ERR=900) '# Frequency in Hz, real part, c. real part, error in %, imag. part, c. imag. part., error in %'
           NOFF=NT(I)/2
           DO J=1,NOFF
              WRITE(UN,'(7G14.6)',ERR=900) TIME(I,J),VOLT(I,J),CVOLT(I,J), &
                   ERROR(I,J),VOLT(I,J+NOFF),CVOLT(I,J+NOFF),ERROR(I,J+NOFF)
           END DO
        ELSE
           WRITE(UN,'(A)',ERR=900) '# Time in s,     voltage, calculated voltage, error in %'
           DO J=1,NT(I)
              WRITE(UN,'(4G14.6)',ERR=900) TIME(I,J),VOLT(I,J),CVOLT(I,J),ERROR(I,J)
           END DO
        ENDIF
        WRITE(6,*) "<I> Saved data set ",I
     ENDDO
  ENDIF
  WRITE(UN,'(A)',ERR=900) '#######################'
  CLOSE(UN,ERR=900)
  GOTO 999
900 WRITE(6,*) '<E> Error writing file'
999 IF (UN/=-1) CALL FREE_UN(UN)
  RETURN
END SUBROUTINE SAVE

SUBROUTINE LOAD(WHAT,NAM)
  USE mar_data
  CHARACTER*200 :: NAM,DUM
  CHARACTER :: DUM1,DUM11
  CHARACTER*2 :: DUM2
  CHARACTER*4 :: FTYPE
  LOGICAL :: ENTER_NAME,TEST,LDUM,READ_DEPTH,LIN
  INTEGER :: UN,I,RFVER,IDUM,WHAT,J1,J2,NOFF,WHAT2
  INTEGER, DIMENSION(:),ALLOCATABLE :: IBUF
  REAL (OP), DIMENSION(:),ALLOCATABLE :: RBUF
  LOGICAL, DIMENSION(:),ALLOCATABLE :: LBUF
  REAL (OP), DIMENSION(:,:),ALLOCATABLE :: RMBUF
  CHARACTER (LEN=4), DIMENSION(:),ALLOCATABLE :: S4BUF
  CHARACTER (LEN=15), DIMENSION(:),ALLOCATABLE :: S15BUF
  REAL (OP) :: RDUM,RDUM2

  PRINT*,NAM
  NOFF=0
  UN=-1
  WHAT2=WHAT
  IF (NDAT==0) THEN
     WHAT2=MIN(WHAT2,3)
  ELSE
     IF (WHAT2==2) THEN
        WRITE(6,*) '<I> Will overwrite data sets in memory'
        WRITE(6,*) '<I> Enter non-sense name to avoid this'
     ENDIF
  END IF
  !IF (ENTER_NAME(NAM,'mar')) THEN
     CALL GET_UN(UN)
     OPEN(UNIT=UN,FILE=NAM,STATUS='OLD',ERR=900)
     READ(UN,*,ERR=900,END=910) DUM
     READ(UN,*,ERR=900,END=910) DUM
     READ(UN,*,ERR=900,END=910) DUM
     READ(UN,*,ERR=900,END=910) DUM
     READ(UN,'(A1,I5)',ERR=900,END=910) DUM1,RFVER
     READ(UN,*,ERR=900,END=910) DUM
     READ(UN,'(A2,A4)',ERR=900,END=910) DUM2,FTYPE
     READ(UN,*,ERR=900,END=910) DUM
     IF (WHAT2<3.AND.FTYPE.NE.'DATA') THEN
        READ(UN,*,ERR=900,END=910) DUM
        READ(UN,'(A1,I5)',ERR=900,END=910) DUM1,NLAY
        IF (ALLOCATED(TNESS)) DEALLOCATE(TNESS,RES,ANIS)
        ALLOCATE(TNESS(NLAY-1),RES(NLAY),ANIS(NLAY))
        READ(UN,*,ERR=900,END=910) DUM
        READ_DEPTH=.TRUE.
        IF (WHAT2==1.AND.NLAY>0) THEN
           READ_DEPTH=LIN('Read water informations from file? ',.FALSE.)
        ENDIF
        IF (READ_DEPTH) THEN
           READ(UN,'(A1,2G14.6)',ERR=900,END=910) DUM1,W_DEPTH,W_RES
        ELSE
           READ(UN,'(A1,2G14.6)',ERR=900,END=910) DUM1,RDUM,RDUM2
        ENDIF
        READ(UN,*,ERR=900,END=910) DUM
        DO I=1,NLAY-1
!           READ(UN,'(A1,I5,3G14.6)',ERR=900,END=910) DUM1,IDUM,RES(I),ANIS(I),TNESS(I)
           READ(UN,*,ERR=900,END=910) DUM1,IDUM,RES(I),ANIS(I),TNESS(I)
        END DO
!        READ(UN,'(A1,I5,2G14.6)',ERR=900,END=910) DUM1,IDUM,RES(NLAY),ANIS(NLAY)
        READ(UN,*,ERR=900,END=910) DUM1,IDUM,RES(NLAY),ANIS(NLAY)
        IF (FTYPE=='MODL'.AND.NDAT>0) THEN
           READ(UN,*,ERR=920,END=920) DUM
           READ(UN,*,ERR=920,END=920) DUM
           READ(UN,'(A1,I5)',ERR=920,END=920) DUM1,IDUM
           READ(UN,*,ERR=920,END=920) DUM
           IF (IDUM==NDAT) THEN
              DO I=1,NDAT
                 READ(UN,'(A2,2G15.6)',ERR=920) DUM2,RDUM,CALF(I)
              ENDDO
           END IF
920        CONTINUE
        ELSEIF (WHAT2==1.AND.FTYPE/='MODL') THEN
           READ(UN,*,ERR=900,END=921) DUM
           READ(UN,*,ERR=900,END=921) DUM
           READ(UN,'(A1,I5)',ERR=900,END=921) DUM1,IDUM
           DO I=1,8
              READ(UN,*,ERR=900,END=921) DUM
           END DO
           LFD=.FALSE.
           IF (RFVER.GT.2) THEN
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,'(A2,L1)',ERR=900,END=910) DUM2,LFD       
           ENDIF
              
           IF (IDUM==NDAT) THEN
              DO J1=1,NDAT
                 DO I=1,5
                    READ(UN,*,ERR=900,END=910) DUM
                 END DO
                 IF (RFVER>1.AND..NOT.LFD) THEN
                    READ(UN,*,ERR=900,END=910) DUM
                    READ(UN,*,ERR=900,END=910) DUM
                 ENDIF   
                 READ(UN,'(A1,I5,G14.6)',ERR=900,END=910) DUM1,J2,RDUM
                 DO I=1,5
                    READ(UN,*,ERR=900,END=910) DUM
                 END DO
                 IF (.NOT.LFD) THEN
                    DO I=1,4
                       READ(UN,*,ERR=900,END=910) DUM
                    END DO
                 ENDIF
                 
                 READ(UN,'(A2,L1,A1,G15.6)',ERR=900,END=910) DUM2,LDUM,DUM1,CALF(J1)
                 DO I=1,J2+1
                    READ(UN,*,ERR=900,END=910) DUM
                 ENDDO
              ENDDO
              WRITE(6,*) '<I> Reading calibration factors'
           END IF
921        CONTINUE
        END IF
        IF (NLAY>0) WRITE(6,*) '<I> Loaded model with ',NLAY,' layers'
        READ(UN,*,ERR=922,END=922) DUM
922     CALL RESET_FIXES
     ELSEIF (FTYPE.NE.'DATA') THEN
        READ(UN,*,ERR=900,END=910) DUM
        READ(UN,'(A1,I5)',ERR=900,END=910) DUM1,IDUM
        DO I=1,4+IDUM
           READ(UN,*,ERR=900,END=910) DUM
        END DO
     ENDIF
     IF (WHAT2>1.AND.FTYPE.NE.'MODL') THEN
        IF (WHAT2.EQ.2) TNT=0
        READ(UN,*,ERR=900,END=910) DUM
        READ(UN,'(A1,I5)',ERR=900,END=910) DUM1,IDUM
        IF (IDUM==0) THEN
           WRITE(6,*) '<I> No data sets found'
        ELSE
           IF (WHAT2==3) THEN
              WRITE(6,*) '<I> Appending ',IDUM,' data sets'
              DO I=1,8
                 READ(UN,*,ERR=900,END=910) DUM
              END DO
              IF (RFVER.GT.2) THEN
                 DO I=1,4
                    READ(UN,*,ERR=900,END=910) DUM
                 END DO
              ENDIF
              ALLOCATE(IBUF(NDAT),RBUF(NDAT),LBUF(NDAT),S4BUF(NDAT),S15BUF(NDAT))
              IBUF=NT
              RBUF=X
              LBUF=SYS_RESP
              S4BUF=TYPE
              S15BUF=SYS_NAME
              IF (ALLOCATED(NT)) DEALLOCATE(NT,X,SYS_RESP,TYPE,SYS_NAME)
              ALLOCATE(NT(NDAT+IDUM),X(NDAT+IDUM),TYPE(NDAT+IDUM),SYS_NAME(NDAT+IDUM),SYS_RESP(NDAT+IDUM))
              NT(1:NDAT)=IBUF
              X(1:NDAT)=RBUF
              SYS_RESP(1:NDAT)=LBUF
              TYPE(1:NDAT)=S4BUF
              SYS_NAME(1:NDAT)=S15BUF
              RBUF=Y
              LBUF=LFCAL
              IF (ALLOCATED(Y)) DEALLOCATE(Y,LFCAL)
              ALLOCATE(Y(NDAT+IDUM),LFCAL(NDAT+IDUM))
              Y(1:NDAT)=RBUF
              LFCAL(1:NDAT)=LBUF
              RBUF=Z
              LBUF=IMP
              IF (ALLOCATED(Z)) DEALLOCATE(Z,IMP)
              ALLOCATE(Z(NDAT+IDUM),IMP(NDAT+IDUM))
              Z(1:NDAT)=RBUF
              IMP(1:NDAT)=LBUF
              RBUF=CHI
              LBUF=STEP_OFF
              IF (ALLOCATED(CHI)) DEALLOCATE(CHI,STEP_OFF)
              ALLOCATE(CHI(NDAT+IDUM),STEP_OFF(NDAT+IDUM))
              CHI(1:NDAT)=RBUF
              STEP_OFF(1:NDAT)=LBUF
              RBUF=DELAY
              LBUF=DERIVATE
              IF (ALLOCATED(DELAY)) DEALLOCATE(DELAY,DERIVATE)
              ALLOCATE(DELAY(NDAT+IDUM),DERIVATE(NDAT+IDUM))
              DELAY(1:NDAT)=RBUF
              DERIVATE(1:NDAT)=LBUF
              RBUF=CALF
              IF (ALLOCATED(CALF)) DEALLOCATE(CALF)
              ALLOCATE(CALF(NDAT+IDUM))
              CALF(1:NDAT)=RBUF
              RBUF=ANG1
              IF (ALLOCATED(ANG1)) DEALLOCATE(ANG1)
              ALLOCATE(ANG1(NDAT+IDUM))
              ANG1(1:NDAT)=RBUF
              RBUF=ANG2
              IF (ALLOCATED(ANG2)) DEALLOCATE(ANG2)
              ALLOCATE(ANG2(NDAT+IDUM))
              ANG2(1:NDAT)=RBUF
              DEALLOCATE(RBUF,IBUF,LBUF,S4BUF,S15BUF)
           ELSE
              WRITE(6,*) '<I> Reading ',IDUM,' data sets'
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,'(A1,G14.6)',ERR=900,END=910) DUM1,TCHI
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,'(A2,A15)',ERR=900,END=910) DUM2,TX_DESIG
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,*,ERR=900,END=910) DUM1,TX_X,TX_Y,TX_Z
!              READ(UN,'(A1,3G14.6)',ERR=900,END=910) DUM1,TX_X,TX_Y,TX_Z
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,'(A1,2G14.6)',ERR=900,END=910) DUM1,DL,CUR
              LFD=.FALSE.
              IF (RFVER.GT.2) THEN
                 READ(UN,*,ERR=900,END=910) DUM
                 READ(UN,*,ERR=900,END=910) DUM
                 READ(UN,*,ERR=900,END=910) DUM
                 READ(UN,'(A2,L1)',ERR=900,END=910) DUM2,LFD       
              ENDIF
              IF (ALLOCATED(NT)) THEN
                 DEALLOCATE(TYPE,NT,CHI,X,Y,Z,ANG1,ANG2,SYS_RESP,SYS_NAME)
                 DEALLOCATE(DELAY,LFCAL,CALF,TIME,VOLT,CVOLT,ERROR)
                 DEALLOCATE(IMP,DERIVATE,STEP_OFF)                 
              ENDIF
              NDAT=0
              ALLOCATE(TYPE(IDUM),NT(IDUM),CHI(IDUM),X(IDUM),Y(IDUM),Z(IDUM))
              ALLOCATE(SYS_RESP(IDUM),SYS_NAME(IDUM),DELAY(IDUM),LFCAL(IDUM))
              ALLOCATE(ANG1(IDUM),ANG2(IDUM),CALF(IDUM),IMP(IDUM))
              ALLOCATE(DERIVATE(IDUM),STEP_OFF(IDUM))
              MAXNT=0
           END IF
           DO I=1,IDUM
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,'(A2,A4)',ERR=900,END=910) DUM2,TYPE(NDAT+I)
              IF (RFVER>1.AND..NOT.LFD) THEN
                 READ(UN,*,ERR=900,END=910) DUM
                 READ(UN,'(A2,L1,A1,L1,A1,L1))',ERR=900,END=910) DUM2,&
                      IMP(NDAT+I),DUM1,STEP_OFF(NDAT+I),DUM11,DERIVATE(NDAT+I)
              ELSE
                 IF (NDAT+I>1) THEN
                    IMP(NDAT+I)=IMP(1)
                    STEP_OFF(NDAT+I)=STEP_OFF(1)
                 ELSE
                    IMP(1)=.FALSE.
                    STEP_OFF(1)=.FALSE.
                 ENDIF
                 IF ((TYPE(NDAT+I)=='EZBP'.OR.TYPE(NDAT+I)=="EZBY".OR. &
                      TYPE(NDAT+I)=="EZBX").AND..NOT.LFD) THEN
                    DERIVATE(NDAT+I)=.TRUE.
                 ELSE
                    DERIVATE(NDAT+I)=.FALSE.
                 ENDIF
              ENDIF
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,'(A1,I5,G14.6)',ERR=900,END=910) DUM1,NT(NDAT+I),CHI(NDAT+I)
              IF (LFD) THEN
                 NOFF=NT(NDAT+I)
                 NT(NDAT+I)=2*NOFF
              END IF
                 
              READ(UN,*,ERR=900,END=910) DUM
!              READ(UN,'(A1,3G14.6)',ERR=900,END=910) DUM1,X(NDAT+I),Y(NDAT+I),Z(NDAT+I)
              READ(UN,*,ERR=900,END=910) DUM1,X(NDAT+I),Y(NDAT+I),Z(NDAT+I)
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,'(A1,2G14.6)',ERR=900,END=910) DUM1,ANG1(NDAT+I),ANG2(NDAT+I)
              IF (.NOT.LFD) THEN
                 READ(UN,*,ERR=900,END=910) DUM
                 READ(UN,'(A2,L1,A1,A15))',ERR=900,END=910) DUM2,SYS_RESP(NDAT+I),DUM1,SYS_NAME(NDAT+I)
                 READ(UN,*,ERR=900,END=910) DUM
                 READ(UN,'(A1,1G14.6)',ERR=900,END=910) DUM1,DELAY(NDAT+I)
              ELSE
                 SYS_RESP(NDAT+I)=.FALSE.
              END IF
              READ(UN,*,ERR=900,END=910) DUM
              READ(UN,'(A2,L1,A1,G15.6)',ERR=900,END=910) DUM2,LFCAL(NDAT+I),DUM1,CALF(NDAT+I)
              READ(UN,*,ERR=900,END=910) DUM
              IF (ALLOCATED(TIME)) THEN
                 ALLOCATE(RMBUF(NDAT+I-1,MAXNT))
                 RMBUF=TIME
                 DEALLOCATE(TIME)
                 ALLOCATE(TIME(NDAT+I,MAX(MAXNT,NT(NDAT+I))))
                 DO J1=1,NDAT+I-1
                    DO J2=1,NT(J1)
                       TIME(J1,J2)=RMBUF(J1,J2)
                    END DO
                 END DO
                 RMBUF=VOLT
                 DEALLOCATE(VOLT)
                 ALLOCATE(VOLT(NDAT+I,MAX(MAXNT,NT(NDAT+I))))
                 DO J1=1,NDAT+I-1
                    DO J2=1,NT(J1)
                       VOLT(J1,J2)=RMBUF(J1,J2)
                    END DO
                 END DO
                 RMBUF=CVOLT
                 DEALLOCATE(CVOLT)
                 ALLOCATE(CVOLT(NDAT+I,MAX(MAXNT,NT(NDAT+I))))
                 DO J1=1,NDAT+I-1
                    DO J2=1,NT(J1)
                       CVOLT(J1,J2)=RMBUF(J1,J2)
                    END DO
                 END DO
                 RMBUF=ERROR
                 DEALLOCATE(ERROR)
                 ALLOCATE(ERROR(NDAT+I,MAX(MAXNT,NT(NDAT+I))))
                 DO J1=1,NDAT+I-1
                    DO J2=1,NT(J1)
                       ERROR(J1,J2)=RMBUF(J1,J2)
                    END DO
                 END DO
                 DEALLOCATE(RMBUF)
              ELSE
                 ALLOCATE(TIME(1,NT(1)),VOLT(1,NT(1)),CVOLT(1,NT(1)),ERROR(1,NT(1)))
              ENDIF
              MAXNT=MAX(MAXNT,NT(NDAT+I))
              IF (LFD) THEN
                 DO J1=1,NOFF
                    READ(UN,*,ERR=900,END=910) TIME(NDAT+I,J1), &
                         VOLT(NDAT+I,J1),CVOLT(NDAT+I,J1),ERROR(NDAT+I,J1), &
                         VOLT(NDAT+I,J1+NOFF),CVOLT(NDAT+I,J1+NOFF), &
                         ERROR(NDAT+I,J1+NOFF)
                    TIME(NDAT+I,J1+NOFF)=TIME(NDAT+I,J1)
                 END DO
              ELSE
                 DO J1=1,NT(NDAT+I)
                    READ(UN,*,ERR=900,END=910) TIME(NDAT+I,J1),VOLT(NDAT+I,J1),CVOLT(NDAT+I,J1),ERROR(NDAT+I,J1)
                 END DO
              ENDIF
              WRITE(6,*) '<I> Read data set ',I
              TNT=TNT+NT(NDAT+I)
           END DO
           NDAT=NDAT+IDUM
           ACT_DAT=1
        ENDIF
     END IF
	 WRITE(6,*) 'system response name ',SYS_NAME,SYS_RESP
     CLOSE(UN,ERR=900)
     IF (FTYPE.NE.'MODL'.AND.WHAT2>1) THEN
        TEST=.FALSE.
        DO I=1,NDAT
           TEST=TEST.OR.SYS_RESP(I)
        ENDDO
        IF (TEST) CALL LOAD_SYS
     ENDIF
  !ENDIF
  GOTO 999
900   WRITE(6,*) 'what,ftype ',WHAT2,FTYPE
      WRITE(6,*) '<E> Error reading file'
  GOTO 999
910 WRITE(6,*) '<E> Error unexpected end of file'
999 IF (UN/=-1) CALL FREE_UN(UN)
  RETURN
! The next lines are just to avoid compiler warnings
  DUM=DUM
  LDUM=LDUM
  DUM11=DUM11
  DUM1=DUM1
  DUM2=DUM2
  RDUM=RDUM
  RDUM2=RDUM2
END SUBROUTINE LOAD

