SUBROUTINE MDED_TX_FHT(SMT)
  USE mar_data
  USE calc_mod

! Routine to calculate the CED forward response for Er and BP
! Originally from Mogilatov; edited by A. Haroon 2013

  INTEGER IJK,I,NDIP,J,NF,IDUMP,L
  REAL (OP) :: ERR,RMIN,SMT,DDUM,REC_POS
  REAL (OP) :: RBUF, DBUF, COEFF1, COEFF2
  COMPLEX (OP),DIMENSION(:),ALLOCATABLE :: FRESP,FRESP_TMP
  REAL (OP), DIMENSION(:),ALLOCATABLE :: OMEGAIR
  REAL (OP), DIMENSION(:),ALLOCATABLE :: CTHICK2,CCON2,CANIS2
  REAL (OP) :: ANG
  LOGICAL :: NEED_STEPOFF
  CHARACTER (LEN=4) :: CTYPE
  DATA ERR/.001/

  REAL (OP) :: HC(250,2),U,SQROOT(2)
  INTEGER :: NC0(2),NC(2),NN
  COMPLEX (OP) :: DC
  COMPLEX (OP) :: F1_N,F2_N,F3_N
  
  EXTERNAL TMFNL,PMFNL,EZML

  COEFF1=0_OP
  COEFF2=0_OP
  RBUF=0_OP
  DBUF=0_OP
  ANG = cos(45*PI/180.0)

!-----------------------
! IMPORTANT PARAMETERS
! Water thickness: W_DEPTH
! Water resistivity: W_RES
! Layer thickness: TNESS - for n layers, n-1 parameters
! Layer resistivity: RES - for n layers, n parameters
! Layer anisotropy: ANIS - for n layers, n parameters
! Current: CUR
! NT: Number of Timepoints
!-----------------------


! Find which Data is calculated - Er, Bp of jz
  COMP=0
  DO IJK=1,NDAT
     CTYPE=TYPE(IJK)
     IF (CTYPE=="MDER") COMP=1
     IF (CTYPE=="MDEZ") COMP=3
  END DO
  IF (COMP==0) RETURN
  NEED_STEPOFF=.FALSE.
   
!-------------------------------------------------------------------  
!                      CREATE NEW MODEL  
!-------------------------------------------------------------------
! Layers have to be included if the DED and/or receiver are not placed 
! at the surface or the seafloor
  CALL FILTER16(4,NC(1),NC0(1),SQROOT(1),HC(1,1)) 
  CALL FILTER16(5,NC(2),NC0(2),SQROOT(2),HC(1,2)) 

  NF = NC_FRT + NT_F - 1  ! NC_FRT: the number of filter coefficients
  
  IF (LFD) NF = NT_F
  
  ALLOCATE(FRESP(NF),FRESP_TMP(NF),OMEGAIR(NF))
  
  DO IJK=1,NDAT

!   ALLOCATE(FRESPNE_1(NF),FRESPNE_2(NF),FRESPNW_1(NF),FRESPNW_2(NF),FRESPW_test(NF))
    DO I = 1,NF
      IF (LFD) THEN
        OMEGAIR(I) = FTIME(I)*2.0*PI
      ELSE
        IDUMP = -NC_FRT+NC0_FRT+I
 !       WRITE(6,*) IDUMP ,SMT,I
        OMEGAIR(I) = (SQRT10**(1-IDUMP))*QSHIFT/SMT
      END IF
    END DO
    
    CDSEA = W_DEPTH
    DC = CMPLX(0.0_OP,0.0_OP)
  
    IF (-TX_Z>W_DEPTH) THEN
      WRITE(6,*) "<E> MCED Transmitter is located in the air!"
      RETURN
    END IF
    
    REC_POS = TX_Z+Z(IJK)
    IF (-REC_POS>W_DEPTH) THEN
      WRITE(6,*) "<E> MCED Receiver is located in the air!"
      STOP
      RETURN
    ENDIF
! Include Layer Boundaries for TX    
    TX_ZC = TX_Z
    RX_ZC = Z(IJK)
    
    IF (TYPE(IJK)=='MDEZ') THEN
        TX_ZC = Z(IJK)+TX_Z
        RX_ZC = -Z(IJK)
    ENDIF
    
     ALLOCATE(CTHICK2(NLAY+2),CCON2(NLAY+3),CANIS2(NLAY+3))
     IF (TX_ZC.LE.0) THEN
        CTHICK2(1)=W_DEPTH+TX_ZC
        CTHICK2(2)=-TX_ZC
        DO I=1,2
           CANIS2(I)=1_OP
           CCON2(I)=1_OP/W_RES
        END DO
        DO I=1,NLAY
           CCON2(I+2)=1.0/RES(I)
           CANIS2(I+2)=ANIS(I)
           IF (I/=NLAY) CTHICK2(I+2)=TNESS(I)
        END DO
     ELSE
        DDUM=0_OP
        I=0
        DO WHILE (TX_ZC>DDUM.AND.I<NLAY-1)
           I=I+1
           DDUM=DDUM+TNESS(I)
        END DO
        CTHICK2(1)=W_DEPTH
        CANIS2(1)=1_OP
        CCON2(1)=1_OP/W_RES
        IF (I==NLAY-1.AND.DDUM<TX_ZC) THEN
           NLU=1+NLAY
           NLL=1
           DDUM=0_OP
           DO I=1,NLAY
              CCON2(I+1)=1_OP/RES(I)
              CANIS2(I+1)=ANIS(I)
              IF (I/=NLAY) THEN
                 CTHICK2(I+1)=TNESS(I)
                 DDUM=DDUM+TNESS(I)
              ELSE
                 CTHICK2(I+1)=TX_ZC-DDUM
              ENDIF
           END DO
           CCON2(NLAY+2)=1./RES(NLAY)
           CANIS2(NLAY+2)=ANIS(NLAY)
        ELSE
           NLU=I+1
           NLL=NLAY-I+1
           CCON2(NLU)=1.0/RES(I)
           CANIS2(NLU)=ANIS(I)
           CTHICK2(NLU)=TNESS(I)-DDUM+TX_ZC
           CCON2(NLU+1)=1.0/RES(I)
           CANIS2(NLU+1)=ANIS(I)
           CTHICK2(NLU+1)=DDUM-TX_ZC
           DO L=1,I-1
              CCON2(L+1)=1.0/RES(L)
              CANIS2(L+1)=ANIS(L)
              CTHICK2(L+1)=TNESS(L)
           END DO
           DO L=I+1,NLAY
              CCON2(L+2)=1.0/RES(L)
              CANIS2(L+2)=ANIS(L)
              IF (L/=NLAY) CTHICK2(L+2)=TNESS(L)
           END DO
        ENDIF
     ENDIF

     ! Insert layer boundary for RX
     CZ=RX_ZC+TX_ZC

     IDUMP=0
     XATM=-W_DEPTH
     ALLOCATE(CTHICK(NLAY+2),CCON(NLAY+3),CANIS(NLAY+3))
     DO L=1,NLAY+1
        IF (IDUMP.EQ.0.AND.XATM.LT.CZ.AND.XATM+CTHICK2(L).GE.CZ) THEN
           IDUMP=1
           CCON(L)=CCON2(L)
           CCON(L+1)=CCON2(L)
           CANIS(L)=CANIS2(L)
           CANIS(L+1)=CANIS2(L)
           CTHICK(L)=CZ-XATM
           CTHICK(L+1)=CTHICK2(L)-CZ+XATM
        ELSE
           CTHICK(L+IDUMP)=CTHICK2(L)
           CCON(L+IDUMP)=CCON2(L)
           CANIS(L+IDUMP)=CANIS2(L)
        END IF
        XATM=XATM+CTHICK2(L)
     END DO
     IF (IDUMP.EQ.0) THEN
        CTHICK(NLAY+2)=CZ-XATM
        CCON(NLAY+2)=CCON2(NLAY+2)
        CANIS(NLAY+2)=CANIS2(NLAY+2)
     END IF
     CCON(NLAY+3)=CCON2(NLAY+2)
     CANIS(NLAY+3)=CANIS2(NLAY+2)
     
     CZ=RX_ZC

     CCON2=CCON
     CTHICK2=CTHICK
     CANIS2=CANIS
     DEALLOCATE(CCON,CTHICK,CANIS)

     ! Check for and kick out zero thickness layers
     IDUMP=0
     DO I=1,NLAY+2
        IF (CTHICK2(I).LT.ETA) IDUMP=IDUMP+1
     ENDDO
     CNL=NLAY+3-IDUMP
     ALLOCATE(CTHICK(CNL-1),CCON(CNL),CANIS(CNL))
     
     IDUMP=0
     DO L=1,NLAY+3
        IF (L.LT.NLAY+3) THEN
           IF (CTHICK2(L).LT.ETA) THEN
              IDUMP=IDUMP+1
           ELSE
              CCON(L-IDUMP)=CCON2(L)
              CANIS(L-IDUMP)=CANIS2(L)
              CTHICK(L-IDUMP)=CTHICK2(L)
           END IF
        ELSE
           CCON(L-IDUMP)=CCON2(L)
           CANIS(L-IDUMP)=CANIS2(L)          
        END IF
     END DO
     DEALLOCATE(CTHICK2,CCON2,CANIS2)
     XATM=-W_DEPTH-TX_ZC
     NLU=0
     NLRX=1
     DIFF=ABS(XATM-CZ)
     DIFF_TX=ABS(XATM)
     !WRITE(6,*) ""
     DO L=1,CNL-1
        !WRITE(6,*) "Vorher: ",DIFF,XATM,NLRX,NLU
        !IF (XATM+ETA<0) NLU=NLU+1
        XATM=XATM+CTHICK(L)
        IF (ABS(XATM)<DIFF_TX) THEN
           NLU=NLU+1
           DIFF_TX=ABS(XATM)
        ENDIF
        IF (ABS(XATM-CZ)<DIFF) THEN
           NLRX=NLRX+1
           DIFF=ABS(XATM-CZ)
        END IF        
        !WRITE(6,*) "Nachher: ",DIFF,XATM,NLRX,NLU
     END DO
     NLL=CNL-NLU
     
     !WRITE(6,*) (1_OP/CCON(I),I=1,CNL)
     !WRITE(6,*) (CANIS(I),I=1,CNL)
     !WRITE(6,*) (CTHICK(I),I=1,CNL-1)
     !WRITE(6,*) IJK,CNL,NLU,NLL,NLRX,CZ
    
    
    
!******************************************************************
!------------------------------------------------------------------

!     Start doing Forward Calculation for each Transmitter Arm
!     Arms are labeled after compass direction (East and  West)

!**************************************************************************
! Left Transmitter Arm
!**************************************************************************

    FRESP=(0.0_OP,0.0_OP)
    DO I=0,NF
        IF (I==0) THEN
           IF (NEED_STEPOFF) THEN
              CYCLE
           ELSE
              OMEGA=0.0_OP
              DC=(0.0_OP,0.0_OP)
           ENDIF
        ELSE
           IF (LFD) THEN
              OMEGA=FTIME(I)*2.0_OP*PI
           ELSE
              IDUMP=-NC_FRT+NC0_FRT+I
              OMEGA=(SQRT10**(1-IDUMP))*QSHIFT/SMT
           ENDIF
        ENDIF
        
        X_TMP = (X(IJK) + DL/2)
        Y_TMP = Y(IJK)

!        WRITE(6,*) X_TMP, Y_TMP
        RMIN=MIN(10000.0,SQRT(X_TMP*X_TMP+Y_TMP*Y_TMP+CZ*CZ))
        NDIP=INT(REAL(DIPOLE_PARTS)*DL/RMIN)
        NDIP=MAX(NDIP,1)



        DO J = 1,NDIP
          XATM = X_TMP + DL/2.0 + REAL(DL)/(2.0*NDIP) - REAL(DL)/REAL(NDIP)*J
          CRHO = SQRT(XATM*XATM + Y_TMP*Y_TMP)
          
          COEFF1 = XATM*XATM/CRHO/CRHO
          COEFF2 = -Y_TMP*Y_TMP/CRHO/CRHO
          

          F1_N = CMPLX(0.0_OP,0.0_OP)
   !       F2_N = CMPLX(0.0_OP,0.0_OP)
          IF (TYPE(IJK) == 'MDER') THEN          
            DO NN = 1,NC(1)
              U = SQROOT(1)**(NN-NC0(1))/CRHO
              CALL PMFNL(U,F2_N)
              CALL TMFNL(U,F3_N)
!                 F3=0.0
              F1_N = F1_N - (COEFF1*F2_N - COEFF2*F3_N) * HC(NN,1)/CRHO
            END DO
          
            DO NN = 1,NC(2)
              U = SQROOT(2)**(NN - NC0(2))/CRHO
              CALL PMFNL(U,F2_N)
              CALL TMFNL(U,F3_N)
!                 F3=0.0
              F1_N = F1_N - U*(COEFF1*F3_N - COEFF2*F2_N) * HC(NN,2)
            END DO
          ELSEIF(TYPE(IJK) == 'MDEZ') THEN
 !           COMP = 2
            DO NN=1,NC(1)
              U = SQROOT(1)**(NN-NC0(1))/CRHO
              CALL EZML(U,F2_N)
              F1_N = F1_N + F2_N*HC(NN,1)
            END DO
          ENDIF
          
          IF (I==0) THEN
            DC = DC + F1_N/NDIP/CRHO
          ELSE
            FRESP(I) = FRESP(I) + F1_N/NDIP/CRHO
          ENDIF
  
         END DO
        
         IF (I>0.AND..NOT.LFD) THEN
           IF (STEP_OFF(IJK)) THEN
              FRESP(I)=DC-FRESP(I)
           ENDIF
           IF (.NOT.SYS_RESP(IJK).OR.SYS_TYPE(IJK)==0) &
                FRESP(I)=FRESP(I)/CMPLX(0.0,OMEGA)
           IF (IMP(IJK)) &
                FRESP(I)=FRESP(I)*CMPLX(0.0,OMEGA)
           IF (DERIVATE(IJK)) &
                FRESP(I)=FRESP(I)*CMPLX(0.0,OMEGA)
         ENDIF
   
       END DO
       
       ! Save response of left transmitter arm ind FRESP_TMP
       DO I = 1,NF
        FRESP_TMP(I) = -FRESP(I)
       ENDDO

!**************************************************************************
! Right Transmitter Arm
!**************************************************************************

    FRESP=(0.0_OP,0.0_OP)
    DO I=0,NF
        IF (I==0) THEN
           IF (NEED_STEPOFF) THEN
              CYCLE
           ELSE
              OMEGA=0.0_OP
              DC=(0.0_OP,0.0_OP)
           ENDIF
        ELSE
           IF (LFD) THEN
              OMEGA=FTIME(I)*2.0_OP*PI
           ELSE
              IDUMP=-NC_FRT+NC0_FRT+I
              OMEGA=(SQRT10**(1-IDUMP))*QSHIFT/SMT
           ENDIF
        ENDIF
        
        X_TMP = X(IJK) - DL/2
        Y_TMP = Y(IJK)
        
        RMIN=MIN(10000.0,SQRT(X_TMP*X_TMP+Y_TMP*Y_TMP+CZ*CZ))
        NDIP=INT(REAL(DIPOLE_PARTS)*DL/RMIN)
        NDIP=MAX(NDIP,1)



        DO J = 1,NDIP

          
          XATM = X_TMP + DL/2.0 + REAL(DL)/(2.0*NDIP) - REAL(DL)/REAL(NDIP)*J
          CRHO = SQRT(XATM*XATM + Y_TMP*Y_TMP)
          
          COEFF1 = XATM*XATM/CRHO/CRHO
          COEFF2 = -Y_TMP*Y_TMP/CRHO/CRHO
          

          F1_N = CMPLX(0.0_OP,0.0_OP)
          IF (TYPE(IJK) == 'MDER') THEN          
            DO NN = 1,NC(1)
              U = SQROOT(1)**(NN-NC0(1))/CRHO
              CALL PMFNL(U,F2_N)
              CALL TMFNL(U,F3_N)
!                 F3=0.0
              F1_N = F1_N - (COEFF1*F2_N - COEFF2*F3_N) * HC(NN,1)/CRHO
            END DO
          
            DO NN = 1,NC(2)
              U = SQROOT(2)**(NN - NC0(2))/CRHO
              CALL PMFNL(U,F2_N)
              CALL TMFNL(U,F3_N)
!                 F3=0.0
              F1_N = F1_N - U*(COEFF1*F3_N - COEFF2*F2_N) * HC(NN,2)
            END DO
          ELSEIF(TYPE(IJK) == 'MDEZ') THEN
            DO NN=1,NC(1)
              U=SQROOT(1)**(NN-NC0(1))/CRHO
              CALL EZML(U,F2_N)
              F1_N = F1_N + F2_N *HC(NN,1)
            END DO
          ENDIF
          
          IF (I==0) THEN
            DC=DC+F1_N/NDIP/CRHO
          ELSE
            FRESP(I)=FRESP(I)+F1_N/NDIP/CRHO
          ENDIF
  
         END DO
        
         IF (I>0.AND..NOT.LFD) THEN
           IF (STEP_OFF(IJK)) THEN
              FRESP(I)=DC-FRESP(I)
           ENDIF
           IF (.NOT.SYS_RESP(IJK).OR.SYS_TYPE(IJK)==0) &
                FRESP(I)=FRESP(I)/CMPLX(0.0,OMEGA)
           IF (IMP(IJK)) &
                FRESP(I)=FRESP(I)*CMPLX(0.0,OMEGA)
           IF (DERIVATE(IJK)) &
                FRESP(I)=FRESP(I)*CMPLX(0.0,OMEGA)
         ENDIF
     
        END DO
        ! Close loop of right Arm
        
        ! Add right arm to left arm
        DO I = 1,NF
          FRESP_TMP(I) = FRESP_TMP(I) + FRESP(I)
        ENDDO
       
       
        FRESP = (0.0_OP,0.0_OP)
        DO I = 1,NF
          FRESP(I) = FRESP_TMP(I)
!          1.03 fits somehow better to the Control transient...I don't know why
          IF (LFD) THEN
            EOUT(IJK,I)=REAL(FRESP(I))
            EOUT(IJK,I+NT_F)=AIMAG(FRESP(I))
          ELSE
            CALL FRTSIN10(FRESP,NF,IJK)
          ENDIF 
        END DO
       
    DEALLOCATE(CTHICK,CCON,CANIS)
  END DO

  DEALLOCATE(FRESP,FRESP_TMP)
        
 END SUBROUTINE MDED_TX_FHT