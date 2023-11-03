SUBROUTINE ENTER_MODEL
  USE mar_data
  INTEGER :: I
  LOGICAL :: ENTER,LIN
  REAL (OP) :: RES_IN,ANI_IN,D1,DMAX,DF,GAMMA,T,D2

  IF (NLAY==0) THEN
     NLAY=20
     IF (NDAT==0) THEN
        W_DEPTH=1000
     ELSE
        W_DEPTH=TX_Z
     ENDIF
     W_RES=0.333333
  END IF
  IF (ALLOCATED(TNESS)) DEALLOCATE(TNESS,ANIS,RES)
  CALL IIN('Enter number of layers ',NLAY)
  IF (NLAY<1) THEN
     WRITE(6,*) "<E> Number of layers must not be less than 1!"
     RETURN
  ENDIF
  ALLOCATE(TNESS(NLAY-1),ANIS(NLAY),RES(NLAY))
  DO I=1,NLAY
    ANIS(I)=1.0
  ENDDO
  CALL RIN('Enter water depth in m ',W_DEPTH)
  CALL RIN('Enter water resistivity in Ohmm ',W_RES) 
  ENTER=(NLAY<=6)
  IF (NLAY>6) ENTER=LIN('Enter individual values for each layer ',ENTER)
  IF (ENTER) THEN
     CALL RINN('Enter resistivities in Ohmm: ',RES,NLAY)
     CALL RINN('Enter anisotropies:          ',ANIS,NLAY)
     CALL RINN('Enter thicknesses in m:      ',TNESS,NLAY-1)
  ELSE
     RES_IN=1.0
     ANI_IN=1.0
     CALL RIN('Enter resistivity in Ohmm: ',RES_IN)
     CALL RIN('Enter anisotropy         : ',ANI_IN)
     DO I=1,NLAY
        RES(I)=RES_IN
        ANIS(I)=ANI_IN
     ENDDO
     D1=0.5
     DMAX=300
     CALL RIN('Enter Thickness of first layer in m: ',D1) 
6    CALL RIN('Enter Depth of last layer in m: ',DMAX) 
     GAMMA=1.0/FLOAT(NLAY-2) 
     DF=D1 
     DMAX=DMAX-(NLAY-2)*DF 
     IF (DMAX<=0) THEN 
        WRITE(6,*)'<E> Please use a bottom depth below ',(NLAY-2)*DF
        DMAX=(NLAY-2)*DF
        GOTO 6
     ENDIF
     T=(ABS(DMAX)/D1)**GAMMA
     TNESS(1)=D1
     DO I = 1,NLAY-2                                     
        D2 = D1*T 
        TNESS(I+1) = D2 - D1 + DF                              
        D1 = D2
     ENDDO
  ENDIF
  CALL PRINT_MODEL(6)
  CALL RESET_FIXES
END SUBROUTINE ENTER_MODEL

SUBROUTINE PRINT_MODEL(UNIT)
  USE mar_data
  REAL (OP) :: DEPTH
  INTEGER I,UNIT
  DEPTH=0
  WRITE(UNIT,*) ' Water depth: ',W_DEPTH,' m, resistivity:',W_RES,' Ohmm'
  WRITE(UNIT,*) ' No.  Resistivity   Anisotropy    Thickness    Depth'
  DO I=1,NLAY-1
     WRITE(UNIT,'(I4,4G14.6)') I,RES(I),ANIS(I),TNESS(I),DEPTH
     DEPTH=DEPTH+TNESS(I)
  ENDDO
  WRITE(UNIT,'(I4,2G14.6,A14,2G14.6)') I,RES(I),ANIS(I),' ',DEPTH
  DO I=1,NDAT
     WRITE(UNIT,'(A,I4,A,G14.6)') '  Calibration factor, data set',I,': ',CALF(I)
  ENDDO
END SUBROUTINE PRINT_MODEL

SUBROUTINE CHANGE_PAR
  USE mar_data
  CHARACTER COM*1,CSTRING*10
  INTEGER I,IBSTR
  WRITE(6,*) 'Change Parameter:'
  WRITE(6,*) ' R - Resistivity    T - Thickness'
  WRITE(6,*) ' A - Anisotropy     C - Calibration'
  WRITE(6,*) ' D - Water depth    W - Water resistivity'
  WRITE(6,*) ' L - List model     E - EXIT'       
  WRITE(6,*) ' R, T, A and C have to be followed by the number of the layer'
20 CSTRING(1:10)='          '
  write(6,*)
  CALL CIN('Change parameter> ',CSTRING)
  CSTRING=TRIM(ADJUSTL(CSTRING))
  COM=CSTRING(1:1)
  CALL LCCONV(COM,1)
  IF(COM.EQ.'R'.OR. COM.EQ.'T'.OR.COM.EQ.'A'.OR.COM.EQ.'C')THEN
     IBSTR=INDEX(CSTRING,' ')
     IF (IBSTR.LE.2) THEN
        IF (((COM=='R'.OR.COM=='A').AND.NLAY==1).OR.(NLAY==2.AND.COM=='T') &
             & .OR.(COM=='C')) THEN
           IF (COM=='C') THEN
              IF (ACT_DAT>0.AND.ACT_DAT<=NDAT) THEN
                 I=ACT_DAT
              ELSE
                 WRITE(6,*) '<E> Wrong active model number'
                 GOTO 20
              END IF
           ELSE
              I=1
           END IF
        ELSE
           WRITE(6,*) '<E> More than one parameter required'
           GOTO 20
        END IF
     ELSEIF(IBSTR.GT.10) THEN
        WRITE(6,*) '<E> More than one parameter required'
        GOTO 20
     ENDIF
     READ(CSTRING(2:10),*,ERR=20)I
  ENDIF
  IF(COM.EQ.'L') THEN
     CALL PRINT_MODEL(6)
     GOTO 20
  ELSEIF(COM.EQ.'E') THEN
     GOTO 999
  ELSEIF(COM.EQ.'D') THEN
     CALL RIN('Enter new water depth ',W_DEPTH)
  ELSEIF(COM.EQ.'W') THEN
     CALL RIN('Enter new water resistivity ',W_RES)
  ELSEIF(COM.EQ.'C') THEN
     CALL RIN('Enter new calibration factor ',CALF(I))
  ELSEIF(COM.EQ.'R') THEN
     IF (I>NLAY.OR.I<=0) GOTO 900
     WRITE(6,*) "Enter new resisitvity for layer ",I
     CALL RIN('Resistivity ',RES(I))
     IF (RES(I)<=0) GOTO 910
  ELSEIF(COM.EQ.'T') THEN
    IF (I>NLAY-1.OR.I<=0) GOTO 900
     WRITE(6,*) "Enter new thickness for layer ",I
     CALL RIN('Thickness ',TNESS(I))
     IF (TNESS(I)<=0) GOTO 910
  elseif(COM.eq.'A') then
     IF (I>NLAY.OR.I<=0) GOTO 900
     WRITE(6,*) "Enter new anisotropy for layer ",I
     CALL RIN('Anisotropy ',ANIS(I))
     IF (ANIS(I)<=0) GOTO 910
  ELSE
     WRITE(6,*) '<E> Wrong input'
     GOTO 20
  ENDIF
  GOTO 20
900 WRITE(6,*) '<E> Layer number out of range'
  GOTO 20
910 WRITE(6,*) '<E> Value out of range'
  GOTO 20
999 RETURN
END SUBROUTINE CHANGE_PAR

SUBROUTINE PRINT_FIXES
  USE mar_data
  INTEGER I
  WRITE(6,*) "Fixes (T=fixed, F=free):"
  WRITE(6,*) "Resistivities: ",(LFIXR(I),I=1,NLAY)
  WRITE(6,*) "Anisotropies:  ",(LFIXA(I),I=1,NLAY)
  WRITE(6,*) "Thicknesses:   ",(LFIXT(I),I=1,NLAY-1)
  WRITE(6,*) "Water depth: ",LFIXWD,"  Water resistivity: ",LFIXWR
  WRITE(6,*) "Calibration factors: ",(LFCAL(I),I=1,NDAT)
END SUBROUTINE PRINT_FIXES

SUBROUTINE RESET_FIXES
  USE mar_data
  INTEGER I
  IF (ALLOCATED(LFIXR)) DEALLOCATE(LFIXR,LFIXT,LFIXA)
  ALLOCATE(LFIXR(NLAY),LFIXA(NLAY),LFIXT(NLAY-1))
  DO I=1,NLAY
     LFIXR(I)=.FALSE.
     LFIXA(I)=.TRUE.
     IF (I<NLAY) LFIXT(I)=.FALSE.
  END DO
  LFIXWD=.TRUE.
  LFIXWR=.TRUE.
!  WRITE(6,*) "<I> Resetting fixes"
END SUBROUTINE RESET_FIXES

SUBROUTINE FIX_PAR
  USE mar_data
  CHARACTER COM*1,CSTRING*10
  INTEGER I,IBSTR
  LOGICAL ALL
  
  WRITE(6,*) 'Fix/Free Parameter:'
  WRITE(6,*) ' R - Resistivity    T - Thickness'
  WRITE(6,*) ' A - Anisotropy     C - Calibration'
  WRITE(6,*) ' D - Water depth    W - Water resistivity'
  WRITE(6,*) ' L - List fixes     Z - Resetting fixes'
  WRITE(6,*) ' E - Exit'
  WRITE(6,*) ' R, T and A have to be followed by the number of the layer'
  WRITE(6,*) ' C has to be followed by the number of the data set'
  WRITE(6,*) ' If no data set/layer number is given, fixes for all layers/data sets are changed.'
20 CSTRING(1:10)='e         '
  write(6,*)
  CALL CIN('Fix> ',CSTRING)
  CSTRING=TRIM(ADJUSTL(CSTRING))
  COM=CSTRING(1:1)
  CALL LCCONV(COM,1)
  IF(COM.EQ.'R'.OR. COM.EQ.'T'.OR.COM.EQ.'A'.OR.COM.EQ.'C')THEN
     IBSTR=INDEX(CSTRING,' ')
     ALL=.FALSE.
     IF (IBSTR.LE.2.OR.IBSTR.GT.10) THEN
        ALL=.TRUE.
     ELSE
        READ(CSTRING(2:10),*,ERR=20)I
     ENDIF
  ENDIF
  IF(COM.EQ.'L') THEN
     CALL PRINT_FIXES
     GOTO 20
  ELSEIF(COM.EQ.'E') THEN
     GOTO 999
  ELSEIF(COM.EQ.'D') THEN
     LFIXWD=.NOT.LFIXWD
  ELSEIF(COM.EQ.'W') THEN
     LFIXWR=.NOT.LFIXWR
  ELSEIF(COM.EQ.'C') THEN
     IF (ALL) THEN
        DO I=1,NDAT
           LFCAL(I)=.NOT.LFCAL(I)
        END DO
     ELSE
        IF (I>0.AND.I<=NDAT) THEN
           LFCAL(I)=.NOT.LFCAL(I)
        ELSE
           WRITE(6,*) '<E> Wrong data set number'
        ENDIF
     END IF
  ELSEIF(COM.EQ.'R') THEN
     WRITE(6,*) "<I> Fixes for the resistivities will be ignored in Occam's inversion"
     IF (ALL) THEN
        DO I=1,NLAY
           LFIXR(I)=.NOT.LFIXR(I)
        END DO
     ELSE
        IF (I>NLAY.OR.I<=0) GOTO 900
        LFIXR(I)=.NOT.LFIXR(I)
     END IF
  ELSEIF(COM.EQ.'T') THEN
     WRITE(6,*) "<I> Fixes for the thicknesses will be ignored in Occam's inversion"
     IF (ALL) THEN
        DO I=1,NLAY-1
           LFIXT(I)=.NOT.LFIXT(I)
        END DO
     ELSE
        IF (I>NLAY-1.OR.I<=0) GOTO 900
        LFIXT(I)=.NOT.LFIXT(I)
     ENDIF
  ELSEIF(COM.EQ.'A') then
     IF (ALL) THEN
        DO I=1,NLAY
           LFIXA(I)=.NOT.LFIXA(I)
        END DO
     ELSE
        IF (I>NLAY.OR.I<=0) GOTO 900
        LFIXA(I)=.NOT.LFIXA(I)
     END IF
  ELSEIF(COM.EQ.'Z') THEN
     CALL RESET_FIXES
  ELSE
     WRITE(6,*) '<E> Wrong input'
     GOTO 20
  ENDIF
  GOTO 20
900 WRITE(6,*) '<E> Layer number out of range'
  GOTO 20
999 RETURN
END SUBROUTINE FIX_PAR
