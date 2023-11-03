SUBROUTINE WEIGHTS
  USE mar_data
  USE inv_para
  USE occam,only:DTRANS,DDTRANS,SCF
  INTEGER :: I,J
  REAL (OP) :: MAXV,VAL,DSCF

  IF (ALLOCATED(WTS)) DEALLOCATE(WTS)
  ALLOCATE(WTS(NDAT,MAXNT))

!  IF (DTRT>1) THEN
!     RS=-10  
!     RV=10000  
!     DO I=1,NDAT
!        DO J=1,NT(I)
!           IF (ABS(VOLT(I,J))>RS.AND.(VOLT(I,J)/=0)) RS=ABS(VOLT(I,J))  
!           IF (ABS(VOLT(I,J))<RV.AND.(VOLT(I,J)/=0)) RV=ABS(VOLT(I,J))  
!        ENDDO
!     ENDDO
!     SCF=ABS(RS/(10*(ALOG10(RS)-ALOG10(RV))))         
!  ENDIF

  ERR_MUL=1_OP
  IF (.NOT.L_FLOOR) THEN
     DO I=1,NDAT
        DO J=1,NT(I)
           IF (ERROR(I,J)/=0) ERR_MUL=MAX(ERR_MUL,ABS(SMALL_ERR/ERROR(I,J)))
        ENDDO
     ENDDO
  ENDIF

  TCHI=0
  DO I=1,NDAT
     DSCF=SCF(I)
     CHI(I)=0
     IF (L_EHIGH) THEN
        MAXV=0.0
        DO J=1,NT(I)
           MAXV=MAX(MAXV,ABS(VOLT(I,J)))
        ENDDO
     ENDIF
     DO J=1,NT(I)
        IF (L_EHIGH) THEN
           VAL=MAXV
        ELSE
           IF (VOLT(I,J)<0) THEN
              VAL=MIN(VOLT(I,J),-ETA**2)
           ELSE
              VAL=MAX(VOLT(I,J),ETA**2)
           ENDIF
        ENDIF
        IF (LUSE_ERR) THEN
           WTS(I,J)=100.0/MAX(SMALL_ERR,ERROR(I,J)*ERR_MUL)/VAL/DDTRANS(VAL,DSCF)
        ELSE
           WTS(I,J)=100/VAL/DDTRANS(VAL,DSCF)
        ENDIF
        WTS(I,J)=ABS(WTS(I,J))
        CHI(I)=CHI(I)+((DTRANS(VOLT(I,J),DSCF)-DTRANS(CVOLT(I,J),DSCF))*WTS(I,J))**2
!        WRITE(6,'(I2,5G14.6)') J,ERROR(I,J)*ERR_MUL,VOLT(I,J),WTS(I,J), &
!             ((DTRANS(VOLT(I,J),DSCF)-DTRANS(CVOLT(I,J),DSCF))*WTS(I,J))**2
     END DO
     TCHI=TCHI+CHI(I)
     CHI(I)=SQRT(CHI(I)/NT(I))
  END DO
  TCHI=SQRT(TCHI/TNT)
END SUBROUTINE WEIGHTS
