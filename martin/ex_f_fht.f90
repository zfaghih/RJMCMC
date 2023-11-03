SUBROUTINE EX_FD_FHT(SMT)
  USE mar_data
  USE calc_mod

  INTEGER :: IJK,I,NDIP,J,NF,IDUMP,L
  REAL (OP) :: RMIN,SMT,DBUF,RBUF,DDUM
  COMPLEX (OP),DIMENSION(:),ALLOCATABLE :: FRESP

  REAL (OP) :: HC(300,2),U,SQROOT(2),XATM,TX_ZC,RX_ZC,DIFF,DIFF_TX
  INTEGER :: NC0(2),NC(2),NN
  LOGICAL :: FOUND_ONE,NEED_STEPOFF,LBUF
  
  COMPLEX (OP) :: F1,F2,F3,DC
  REAL (OP), PARAMETER :: SPD = 20.0_OP
  REAL (OP) :: COEFF1,COEFF2
  REAL (OP), DIMENSION(:),ALLOCATABLE :: CTHICK2,CCON2,CANIS2

  EXTERNAL TMFNL,PMFNL

  COEFF1=0_OP
  COEFF2=0_OP
  RBUF=0_OP
  DBUF=0_OP

! Set filter coefficients for needed Hankel transforms
! Set 1 is J1 (needed for poloroidal fields)
! Set 2 is J1 prime (needed for toroidal fields)

  CALL FILTER16(4,NC(1),NC0(1),SQROOT(1),HC(1,1))
  CALL FILTER16(5,NC(2),NC0(2),SQROOT(2),HC(1,2))
 
  NF=NC_FRT+NT_F-1 !! NC_FRT the number of filter coefficients
  
  FOUND_ONE=.FALSE.
  NEED_STEPOFF=.FALSE.
  DO IJK=1,NDAT
     IF (TYPE(IJK)/='EXEX'.AND.TYPE(IJK)/='EXEY' &
          .AND.TYPE(IJK)/='EXBZ'.AND.TYPE(IJK)/='EXEZ') CYCLE
     FOUND_ONE=.TRUE.
     NEED_STEPOFF=NEED_STEPOFF.OR.STEP_OFF(IJK)
  END DO

  IF (.NOT.FOUND_ONE) RETURN

  LBUF=.FALSE.
  IF (NEED_STEPOFF.AND.W_DEPTH==0.0_OP) THEN
     LBUF=.TRUE.
     DBUF=W_DEPTH
     RBUF=W_RES
     W_DEPTH=1E-6_OP
     W_RES=1e08_OP
  ENDIF

  ALLOCATE(FRESP(NF))
  CDSEA=W_DEPTH

  IF (LFD) NF=NT_F

  DO IJK=1,NDAT
     FRESP=(0.0_OP,0.0_OP)
     DC=(0.0_OP,0.0_OP)
     IF (TYPE(IJK)/='EXEX'.AND.TYPE(IJK)/='EXEY' &
          .AND.TYPE(IJK)/='EXBZ'.AND.TYPE(IJK)/='EXEZ') CYCLE

! Constructing new model above and below TX
     TX_ZC=TX_Z
     RX_ZC=Z(IJK)
     IF (TYPE(IJK)=='EXBZ'.OR.TYPE(IJK)=='EXEZ') THEN
        TX_ZC=Z(IJK)+TX_Z
        RX_ZC=-Z(IJK)
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
     
     DO I=0,NF
        IF (I==0) THEN
           IF (.NOT.NEED_STEPOFF) THEN
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
        RMIN=MIN(10000.0,SQRT(X(IJK)*X(IJK)+Y(IJK)*Y(IJK)+CZ*CZ))
        NDIP=INT(REAL(DIPOLE_PARTS)*DL/RMIN)
        NDIP=MAX(NDIP,1)
        
        DO J=1,NDIP
           XATM=X(IJK)+DL/2.0+REAL(DL)/(2.0*NDIP)-REAL(DL)/REAL(NDIP)*J
           CRHO=SQRT(XATM*XATM+Y(IJK)*Y(IJK))
           IF (TYPE(IJK)=='EXEX') THEN
              COEFF1=XATM*XATM/CRHO/CRHO
              COEFF2=-Y(IJK)*Y(IJK)/CRHO/CRHO
           ELSEIF (TYPE(IJK)=='EXEY') THEN
              COEFF1=XATM*Y(IJK)/CRHO/CRHO
              COEFF2=COEFF1
           ELSEIF (TYPE(IJK)=='EXBZ') THEN
              COEFF1=Y(IJK)*Y(IJK)/CRHO/CRHO
           ENDIF
          F1=CMPLX(0.0_OP,0.0_OP)
           IF (TYPE(IJK)=='EXEX'.OR.TYPE(IJK)=='EXEY') THEN
              COMP=1
              DO NN=1,NC(1)
                 U=SQROOT(1)**(NN-NC0(1))/CRHO
                 CALL PMFNL(U,F2)
                 CALL TMFNL(U,F3)
!                 F3=0.0
                 F1=F1-(COEFF1*F2-COEFF2*F3)*HC(NN,1)/CRHO
              END DO
              DO NN=1,NC(2)
                 U=SQROOT(2)**(NN-NC0(2))/CRHO
                 CALL PMFNL(U,F2)
                 CALL TMFNL(U,F3)
!                 F3=0.0
                 F1=F1-U*(COEFF1*F3-COEFF2*F2)*HC(NN,2)
              END DO
           ELSEIF (TYPE(IJK).EQ.'EXEZ') THEN
              COMP=2
              DO NN=1,NC(1)
                 U=SQROOT(1)**(NN-NC0(1))/CRHO
                 CALL EZML(U,F2)
                 F1=F1+F2*HC(NN,1)
              END DO
           ELSEIF (TYPE(IJK).EQ.'EXBZ') THEN
              COMP=1
              DO NN=1,NC(1)
                 U=SQROOT(1)**(NN-NC0(1))/CRHO
                 CALL PMFNL(U,F2)
                 F1=F1+F2*HC(NN,1)*U*U
              END DO
              F1=F1/CMPLX(0.0,OMEGA)*COEFF1
           ENDIF          
           IF (I==0) THEN
              DC=DC+F1/NDIP/CRHO
           ELSE
              FRESP(I)=FRESP(I)+F1/NDIP/CRHO
           ENDIF
        END DO
        IF (I>0.AND..NOT.LFD) THEN
           IF (STEP_OFF(IJK)) THEN
              FRESP(I)=DC-FRESP(I)
           ENDIF
           IF (.NOT.SYS_RESP(IJK)) &
                FRESP(I)=FRESP(I)/CMPLX(0.0,OMEGA)
           IF (IMP(IJK)) &
                FRESP(I)=FRESP(I)*CMPLX(0.0,OMEGA)
           IF (DERIVATE(IJK)) &
                FRESP(I)=FRESP(I)*CMPLX(0.0,OMEGA)
        ENDIF
        IF (LFD) THEN
           EOUT(IJK,I)=REAL(FRESP(I))
           EOUT(IJK,I+NT_F)=AIMAG(FRESP(I))
        ELSE
           CALL FRTSIN10(FRESP,NF,IJK)
        ENDIF 
     END DO
     DEALLOCATE(CTHICK,CCON,CANIS)
  END DO

  IF (LBUF) THEN
     W_DEPTH=DBUF
     W_RES=RBUF
  ENDIF
  DEALLOCATE(FRESP)

END SUBROUTINE EX_FD_FHT

SUBROUTINE TMFNL(ALDA,ANSTM)
  USE calc_mod
  USE mat_funcs,only:CTANH,CSECH
!
  INTEGER :: J,K
  REAL (OP) :: AL2,FN,ALDA,MU
  COMPLEX (OP) :: YOLD,THETA,TNFNEW,YMINUS,YPLUS,ANSTM,IW
  COMPLEX (OP), DIMENSION(:), ALLOCATABLE :: QL,QU

  ALLOCATE(QL(NLL),QU(NLU+1))
  MU=PI*4.D-07
  AL2=ALDA*ALDA
  IW=CMPLX(0.0,OMEGA)

! First lower part:
  FN=CANIS(CNL)
  THETA=CSQRT(AL2*FN*FN+IW*CCON(CNL)*MU)
  QL(1)=THETA/CCON(CNL)
!
  K=1
  DO J=CNL-1,NLU+1,-1
     K=K+1
     YOLD=QL(K-1)
     FN=CANIS(J)
     THETA=CSQRT(AL2*FN*FN+IW*CCON(J)*MU)
     TNFNEW=CTANH(THETA*CTHICK(J))
     QL(K)=THETA/CCON(J)*(YOLD*CCON(J)+THETA*TNFNEW)/(THETA+YOLD*TNFNEW*CCON(J))
  END DO
  YMINUS=QL(K)

  THETA=CSQRT(AL2*FN*FN+IW*1D-8*MU)
  QU(1)=THETA*1D+8

  DO J=1,NLU
     YOLD=QU(J)
     FN=CANIS(J)
     THETA=CSQRT(AL2*FN*FN+IW*CCON(J)*MU)
     TNFNEW=CTANH(THETA*CTHICK(J))
     QU(J+1)=THETA/CCON(J)*(YOLD*CCON(J)+THETA*TNFNEW)/(THETA+YOLD*TNFNEW*CCON(J))
  ENDDO
  YPLUS=QU(NLU+1)
  IF (COMP==1) THEN
     ANSTM = YPLUS*YMINUS/(2_OP*PI*(YPLUS+YMINUS))
  ELSEIF (COMP==2) THEN
     ANSTM = MU*IW/(2_OP*PI*(YPLUS+YMINUS))
  ENDIF
  
  IF (CZ<0) THEN
! Rxs higher than Tx        
     DO J=NLU,NLRX,-1
        FN=CANIS(J)
        THETA=CSQRT(AL2*FN*FN+IW*CCON(J)*MU)
        TNFNEW=THETA*CTHICK(J)
        ANSTM=ANSTM*CSECH(TNFNEW)*QU(J)*CCON(J)/(THETA*CTANH(TNFNEW)+QU(J)*CCON(J))
     ENDDO
  ELSE
! Rxs deeper than Tx        
     DO J=NLU+1,NLRX-1
        FN=CANIS(J)
        THETA=CSQRT(AL2*FN*FN+IW*CCON(J)*MU)
        TNFNEW=THETA*CTHICK(J)
        ANSTM=ANSTM*CSECH(TNFNEW)*QL(CNL-J)*CCON(J)/(THETA*CTANH(TNFNEW)+QL(CNL-J)*CCON(J))
     ENDDO
  ENDIF
  
  DEALLOCATE(QL,QU)

  RETURN
END SUBROUTINE TMFNL

SUBROUTINE PMFNL(ALDA,ANSTM)
  USE calc_mod
  USE mat_funcs,only:CTANH,CSECH
!
  INTEGER :: J,K
  COMPLEX (OP) :: QOLD,THETA,TNFNEW,QMINUS,QPLUS,ANSTM,IW
  REAL (OP) :: ALDA,AL2,MU
  COMPLEX (OP), DIMENSION(:), ALLOCATABLE :: QL,QU

  ALLOCATE(QL(NLL),QU(NLU+1))
!
  MU=PI*4.D-07
  AL2=ALDA*ALDA
  IW=CMPLX(0.0,OMEGA)

! First lower part:
  THETA=CSQRT(AL2+IW*CCON(CNL)*MU)
  QL(1)=MU/THETA
  K=1
  DO J=CNL-1,NLU+1,-1
     K=K+1
     QOLD=QL(K-1)
     THETA=CSQRT(AL2+IW*CCON(J)*MU)
     TNFNEW=CTANH(THETA*CTHICK(J))
     QL(K)=MU/THETA*(THETA*QOLD+MU*TNFNEW)/(MU+THETA*QOLD*TNFNEW)
  END DO
  QMINUS=QL(K)

! Now upper part:
  THETA=CSQRT(AL2+IW*1D-08*MU)
  QU(1)=MU/THETA
  DO J=1,NLU
     THETA=CSQRT(AL2+IW*MU*CCON(J))
     TNFNEW=CTANH(THETA*CTHICK(J))
     QU(J+1)=MU/THETA*(THETA*QU(J)+MU*TNFNEW)/(MU+QU(J)*THETA*TNFNEW)
  END DO
  QPLUS=QU(NLU+1)

  IF (COMP==1) THEN
     ANSTM = IW*QPLUS*QMINUS/(2_OP*PI*(QPLUS+QMINUS))
  ELSEIF (COMP==2) THEN
     ANSTM = MU/(2_OP*PI*(QPLUS+QMINUS))
  ENDIF

  IF (CZ<0_OP) THEN
! Rxs higher than Tx       
     DO J=NLU,NLRX,-1
        THETA=CSQRT(AL2+IW*MU*CCON(J))
        TNFNEW=THETA*CTHICK(J)
        ANSTM=ANSTM*THETA*CSECH(TNFNEW)*QU(J)/(THETA*QU(J)+MU*CTANH(TNFNEW))
     ENDDO
  ELSE
! Rxs deeper than Tx        
     DO J=NLU+1,NLRX-1
        THETA=CSQRT(AL2+IW*MU*CCON(J))
        TNFNEW=THETA*CTHICK(J)
        ANSTM=ANSTM*THETA*CSECH(TNFNEW)*QL(CNL-J)/(THETA*QL(CNL-J)+MU*CTANH(TNFNEW))
     ENDDO
  ENDIF
  RETURN
END SUBROUTINE PMFNL

