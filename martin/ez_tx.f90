SUBROUTINE EZ_TX_FHT(SMT)
  USE mar_data
  USE calc_mod

! routine for calculating fields of an EZ-TX in a stratified medium
! written Feb, 2006 by C. Scholl

  INTEGER IJK,I,NDIP,J,NF,K,IDUMP,L,IDUMP2,NLDIP,KK,N1
  REAL (OP) :: ERR,RS,RMIN,SMT,DDUM,D1,D2
  COMPLEX (OP),DIMENSION(:),ALLOCATABLE :: FRESP
  REAL (OP), DIMENSION(:),ALLOCATABLE :: OMEGAIR,DPIECE1,DPIECE2,DPL,DPZ
  REAL (OP), DIMENSION(:),ALLOCATABLE :: CTHICK2,CCON2,CANIS2
  REAL (OP) :: CTZ,SINT,COST
  LOGICAL :: FOUND
  CHARACTER (LEN=4) :: CTYPE
  DATA ERR/.001/

  REAL (OP) :: HC(250,2),U,SQROOT(2)
  INTEGER :: NC0(2),NC(2),NN,FCN
  COMPLEX (OP) :: F1,F2,DC,FKERN
  
  EXTERNAL EZML
     
  COMP=0
  DO IJK=1,NDAT
     CTYPE=TYPE(IJK)
     IF (CTYPE=="EZEZ") COMP=3
     IF (CTYPE=="EZEY") COMP=2
     IF (CTYPE=="EZEX") COMP=2
     IF (CTYPE=="EZER") COMP=2
     IF (CTYPE=="EZBP") COMP=1
     IF (CTYPE=="EZBY") COMP=1
     IF (CTYPE=="EZBX") COMP=1
  ENDDO
  IF (COMP==0) RETURN

  CALL FILTER16(3,NC(1),NC0(1),SQROOT(1),HC(1,1))
  CALL FILTER16(4,NC(2),NC0(2),SQROOT(2),HC(1,2))

! First the TX-Dipole is divided into one dipole for each layer crossed
  D1=TX_Z-DL/2.0
  D2=TX_Z+DL/2.0
  IF (-D1>W_DEPTH) THEN
     WRITE(6,*) "<E> Transmitting dipole in air!"
     RETURN
  ENDIF
  NLDIP=0
  IF (D1<0) NLDIP=NLDIP+1
  DDUM=0.0
  DO J=1,NLAY-1
     IF ((D1>=DDUM.AND.D1<=DDUM+TNESS(J)).OR.(NLDIP>0.AND.D2>=DDUM)) THEN
        NLDIP=NLDIP+1
     ENDIF
     DDUM=DDUM+TNESS(J)
  ENDDO
  IF (D2>DDUM) NLDIP=NLDIP+1
  IDUMP=1
  ALLOCATE(DPIECE1(NLDIP),DPIECE2(NLDIP),DPL(NLDIP),DPZ(NLDIP))
  IF (D1<0) THEN
     DPIECE1(1)=D1
     DPIECE2(1)=MIN(D2,0.0)
     IDUMP=2
  ENDIF
  DO J=IDUMP,NLDIP
     IF (J>1) THEN
        DPIECE1(J)=DPIECE2(J-1)
     ELSE
        DPIECE1(J)=D1
     ENDIF
     DDUM=0.0
     DO I=1,NLAY-1
        IF (DPIECE1(J)>=DDUM.AND.DPIECE1(J)<=DDUM+TNESS(I)) THEN
           IF (D2>=DDUM+TNESS(I)) THEN
              DPIECE2(J)=DDUM+TNESS(I)
           ELSE
              DPIECE2(J)=D2
           ENDIF
        ENDIF
        DDUM=DDUM+TNESS(I)
     END DO
     IF (D2>DDUM.AND.J.EQ.NLDIP) DPIECE2(J)=D2
  ENDDO
  DO J=1,NLDIP
     DPL(J)=DPIECE2(J)-DPIECE1(J)
     DPZ(J)=DPIECE1(J)+DPL(J)/2.0
  END DO

  DO IJK=1,NDAT
     CTYPE=TYPE(IJK)
     COMP=0
     IF (CTYPE=="EZEZ") COMP=3
     IF (CTYPE=="EZEY") COMP=2
     IF (CTYPE=="EZEX") COMP=2
     IF (CTYPE=="EZER") COMP=2
     IF (CTYPE=="EZBP") COMP=1
     IF (CTYPE=="EZBY") COMP=1
     IF (CTYPE=="EZBX") COMP=1
     IF (COMP==0) CYCLE

     NF=NC_FRT+NT_F-1 !! NC_FRT the number of filter coefficients
     IF (LFD) NF=NT_F
     ALLOCATE(FRESP(NF),OMEGAIR(NF))
     DO I=1,NF
        IF (LFD) THEN
           OMEGAIR(I)=FTIME(I)*2.0*PI
        ELSE
           IDUMP=-NC_FRT+NC0_FRT+I
           OMEGAIR(I)=(SQRT10**(1-IDUMP))*QSHIFT/SMT
        END IF
     END DO
     CDSEA=W_DEPTH
     DC=CMPLX(0.0_OP,0.0_OP)
     DO K=1,NF
        OMEGA=OMEGAIR(K)
        FRESP(K)=CMPLX(0.0_OP,0.0_OP)
        DO KK=1,NLDIP
           CZ=Z(IJK)+(TX_Z-DPZ(KK))
           RS=SQRT(X(IJK)*X(IJK)+CZ*CZ+Y(IJK)*Y(IJK))
           RMIN=MIN(10000.0,RS)
           NDIP=INT(REAL(DIPOLE_PARTS)*DPL(KK)/RMIN)
           NDIP=MAX(NDIP,1)
           DO J=1,NDIP
              CTZ=DPZ(KK)+DPL(KK)/2.0+REAL(DPL(KK))/(2.0*NDIP)-REAL(DPL(KK))/REAL(NDIP)*J
              CZ=-CTZ+Z(IJK)+TX_Z
              CRHO=SQRT(X(IJK)*X(IJK)+Y(IJK)*Y(IJK))
! Build new model 
              IF (-CTZ>W_DEPTH) THEN
                 WRITE(6,*) "<E> Transmitting dipole in air!"
                 RETURN
              ENDIF
              ALLOCATE(CTHICK2(NLAY+1),CCON2(NLAY+2),CANIS2(NLAY+2))
              IF (CTZ<=0) THEN
                 NLU=1
                 NLL=NLAY+1
                 CTHICK2(1)=W_DEPTH+CTZ
                 CTHICK2(2)=-CTZ
                 CANIS2(1)=1.0
                 CANIS2(2)=1.0
                 CCON2(1)=1.0/W_RES
                 CCON2(2)=1.0/W_RES
                 DO I=1,NLAY
                    CCON2(I+2)=1.0/RES(I)
                    CANIS2(I+2)=ANIS(I)
                    IF (I/=NLAY) CTHICK2(I+2)=TNESS(I)
                 END DO
              ELSE
                 DDUM=0
                 I=0
                 DO WHILE (CTZ.GE.DDUM.AND.I<NLAY-1)
                    I=I+1
                    DDUM=DDUM+TNESS(I)
                 END DO
                 CTHICK2(1)=W_DEPTH
                 CANIS2(1)=1.0
                 CCON2(1)=1.0/W_RES
                 IF (I==NLAY-1.AND.DDUM<CTZ) THEN
                    NLU=1+NLAY
                    NLL=1
                    DDUM=0.0
                    DO I=1,NLAY
                       CCON2(I+1)=1.0/RES(I)
                       CANIS2(I+1)=ANIS(I)
                       IF (I/=NLAY) THEN
                          CTHICK2(I+1)=TNESS(I)
                          DDUM=DDUM+TNESS(I)
                       ELSE
                          CTHICK2(I+1)=CTZ-DDUM
                       ENDIF
                    END DO
                    CCON2(NLAY+2)=1./RES(NLAY)
                    CANIS2(NLAY+2)=ANIS(NLAY)
                 ELSE
                    NLU=I+1
                    NLL=NLAY-I+1
                    CCON2(NLU)=1.0/RES(I)
                    CANIS2(NLU)=ANIS(I)
                    CTHICK2(NLU)=TNESS(I)-DDUM+CTZ
                    CCON2(NLU+1)=1.0/RES(I)
                    CANIS2(NLU+1)=ANIS(I)
                    CTHICK2(NLU+1)=DDUM-CTZ
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
              ! Check for and kick out zero thickness layer
              CNL=NLAY+2
              
              IDUMP=0
              DO I=1,NLAY+1
                 IF (CTHICK2(I)<ETA) IDUMP=I
              ENDDO
              
              IF (IDUMP==0) THEN
                 CNL=NLAY+2
                 ALLOCATE(CTHICK(CNL-1),CCON(CNL),CANIS(CNL))
                 CTHICK=CTHICK2
                 CCON=CCON2
                 CANIS=CANIS2
              ELSE
                 CNL=NLAY+1
                 NLL=NLL-1
                 ALLOCATE(CTHICK(CNL-1),CCON(CNL),CANIS(CNL))
                 L=0
                 DO I=1,CNL+1
                    IF (I/=IDUMP) THEN
                       L=L+1
                       IF (L<CNL) CTHICK(L)=CTHICK2(I)
                       CCON(L)=CCON2(I)
                       CANIS(L)=CANIS2(I)
                    END IF
                 END DO
              ENDIF
              DEALLOCATE(CTHICK2,CCON2,CANIS2)
              
              ! Now include layer boundary for Rx
              IF (CZ>0) THEN
                 DDUM=0
                 FOUND=.FALSE.
                 DO I=NLU+1,CNL-1
                    DDUM=DDUM+CTHICK(I)
                    IF (ABS(CZ-DDUM)<ETA) THEN
                       FOUND=.TRUE.
                       NLRX=I+1
                       EXIT
                    ELSEIF(DDUM>CZ) THEN
                       NLRX=I+1
                       ALLOCATE(CTHICK2(CNL-1),CCON2(CNL),CANIS2(CNL))
                       CTHICK2=CTHICK
                       CCON2=CCON
                       CANIS2=CANIS
                       DEALLOCATE(CTHICK,CCON,CANIS)
                       ALLOCATE(CTHICK(CNL),CCON(CNL+1),CANIS(CNL+1))
                       IDUMP=0
                       DO L=1,CNL
                          IF (L==I) THEN
                             IDUMP=1
                             CCON(L)=CCON2(L)
                             CANIS(L)=CANIS2(L)
                             CCON(L+1)=CCON2(L)
                             CANIS(L+1)=CANIS2(L)
                             CTHICK(L+1)=DDUM-CZ
                             CTHICK(L)=CZ-DDUM+CTHICK2(L)
                          ELSE
                             CCON(L+IDUMP)=CCON2(L)
                             CANIS(L+IDUMP)=CANIS2(L)
                             IF (L/=CNL) CTHICK(L+IDUMP)=CTHICK2(L)
                          ENDIF
                       END DO
                       DEALLOCATE(CTHICK2,CCON2,CANIS2)
                       CNL=CNL+1
                       NLL=NLL+1
                       FOUND=.TRUE.
                       EXIT
                    ENDIF
                 ENDDO
                 IF (.NOT.FOUND) THEN
                    NLRX=CNL+1
                    ALLOCATE(CTHICK2(CNL-1),CCON2(CNL),CANIS2(CNL))
                    CTHICK2=CTHICK
                    CCON2=CCON
                    CANIS2=CANIS
                    DEALLOCATE(CTHICK,CCON,CANIS)
                    CNL=CNL+1
                    NLL=NLL+1
                    ALLOCATE(CTHICK(CNL-1),CCON(CNL),CANIS(CNL))
                    CCON(1:CNL-1)=CCON2
                    CANIS(1:CNL-1)=CANIS2
                    CTHICK(1:CNL-2)=CTHICK2
                    CCON(CNL)=CCON2(CNL-1)
                    CANIS(CNL)=CANIS2(CNL-1)
                    CTHICK(CNL-1)=CZ-DDUM
                    DEALLOCATE(CTHICK2,CCON2,CANIS2)
                 ENDIF
              ELSEIF (CZ<0) THEN
                 DDUM=0
                 FOUND=.FALSE.
                 IDUMP2=NLU
                 DO I=IDUMP2,1,-1
                    DDUM=DDUM+CTHICK(I)
                    IF (ABS(-CZ-DDUM)<ETA) THEN
                       FOUND=.TRUE.
                       NLRX=I
                       EXIT
                    ELSEIF(DDUM>-CZ) THEN
                       NLRX=I+1
                       ALLOCATE(CTHICK2(CNL-1),CCON2(CNL),CANIS2(CNL))
                       CTHICK2=CTHICK
                       CCON2=CCON
                       CANIS2=CANIS
                       DEALLOCATE(CTHICK,CCON,CANIS)
                       ALLOCATE(CTHICK(CNL),CCON(CNL+1),CANIS(CNL+1))
                       IDUMP=0
                       DO L=1,CNL
                          IF (L==I) THEN
                             IDUMP=1
                             CCON(L)=CCON2(L)
                             CANIS(L)=CANIS2(L)
                             CCON(L+1)=CCON2(L)
                             CANIS(L+1)=CANIS2(L)
                             CTHICK(L+1)=CTHICK2(L)-DDUM-CZ
                             CTHICK(L)=DDUM+CZ
                          ELSE
                             CCON(L+IDUMP)=CCON2(L)
                             CANIS(L+IDUMP)=CANIS2(L)
                             IF (L/=CNL) CTHICK(L+IDUMP)=CTHICK2(L)
                          ENDIF
                       END DO
                       DEALLOCATE(CTHICK2,CCON2,CANIS2)
                       CNL=CNL+1
                       NLU=NLU+1
                       FOUND=.TRUE.
                       EXIT
                    ENDIF
                 ENDDO
                 IF (.NOT.FOUND) THEN
                    WRITE(6,*) "<I> Measuring in air?"
                    NLRX=1
                    ALLOCATE(CTHICK2(CNL-1),CCON2(CNL),CANIS2(CNL))
                    CTHICK2=CTHICK
                    CCON2=CCON
                    CANIS2=CANIS
                    DEALLOCATE(CTHICK,CCON,CANIS)
                    CNL=CNL+1
                    NLL=NLL+1
                    ALLOCATE(CTHICK(CNL-1),CCON(CNL),CANIS(CNL))
                    CCON(2:CNL)=CCON2
                    CANIS(2:CNL)=CANIS2
                    CTHICK(2:CNL-1)=CTHICK2
                    CCON(1)=1e-8
                    CANIS(1)=1
                    CTHICK(1)=-CZ-DDUM
                    DEALLOCATE(CTHICK2,CCON2,CANIS2)
                 ENDIF
              ELSE
                 NLRX=NLU+1
              ENDIF
              
!              IF (K==1.OR.K==10) THEN
!                 WRITE(6,*) CNL,NLRX,NLL,NLU,CZ
!                 WRITE(6,*) CCON
!                 WRITE(6,*) CTHICK
!                 WRITE(6,*) "COMP:",COMP
!             ENDIF

              F1=CMPLX(0.0_OP,0.0_OP)
              FCN=2
              IF (COMP==3) FCN=1
              DO NN=1,NC(FCN)
                 U=SQROOT(FCN)**(NN-NC0(FCN))/CRHO
                 CALL EZML(U,FKERN)
                 F1=F1+FKERN*HC(NN,FCN)
              END DO

              FRESP(K)=FRESP(K)+F1*(DPL(KK)/NDIP)/DL/CRHO

              IF (K==1.AND.STEP_OFF(IJK)) THEN
                 F2=CMPLX(0.0_OP,0.0_OP)
                 OMEGA=0.0
                 DO NN=1,NC(FCN)
                    U=SQROOT(FCN)**(NN-NC0(FCN))/CRHO
                    CALL EZML(U,FKERN)
                    F2=F2+FKERN*HC(NN,FCN)
                 END DO
                 DC=DC+F2*(DPL(KK)/NDIP)/DL/CRHO
                 OMEGA=OMEGAIR(K)
              ENDIF

              DEALLOCATE(CTHICK,CCON,CANIS)
           END DO
        ENDDO

!        AR=REAL(FRESP(K))
!        AI=AIMAG(FRESP(K))
!        PHASE=ATAN2(AI,AR)*180.0/PI
!        DO WHILE (PHASE>0)
!           PHASE=PHASE-360.0
!        END DO

        IF (.NOT.LFD) THEN
           IF (STEP_OFF(IJK)) THEN
              FRESP(K)=DC-FRESP(K)
           ENDIF
           IF (DERIVATE(IJK)) FRESP(K)=FRESP(K)*CMPLX(0.0,OMEGAIR(K))
           IF (.NOT.SYS_RESP(IJK)) &
                FRESP(K)=FRESP(K)/CMPLX(0.0,OMEGAIR(K))
           IF (IMP(IJK)) &
                FRESP(K)=FRESP(K)*CMPLX(0.0,OMEGAIR(K))
        END IF
     END DO
     
!     DO I=1,NF
!        WRITE(10,*) IJK,OMEGAIR(I),REAL(FRESP(I)),AIMAG(FRESP(I))
!     ENDDO
     N1=NT_F
     IF (LFD) THEN
        DO I=1,NT_F
           EOUT(IJK,I)=REAL(FRESP(I))
           EOUT(IJK,I+NT_F)=AIMAG(FRESP(I))
        END DO
        N1=NT_F*2
     ELSE
        CALL FRTSIN10(FRESP,NF,IJK)
     ENDIF
     IF (CTYPE=="EZEX".OR.CTYPE=="EZBY") THEN
        COST=X(IJK)/CRHO
        DO I=1,N1
           EOUT(IJK,I)=EOUT(IJK,I)*COST
        ENDDO
     ENDIF
     IF (CTYPE=="EZEY".OR.CTYPE=="EZBX") THEN
        SINT=Y(IJK)/CRHO
        DO I=1,N1
           EOUT(IJK,I)=EOUT(IJK,I)*SINT
        ENDDO
     ENDIF
     DEALLOCATE(FRESP,OMEGAIR)
  END DO
  DEALLOCATE(DPIECE1,DPIECE2,DPL,DPZ)
END SUBROUTINE EZ_TX_FHT

SUBROUTINE EZML(ALDA,ANSTM)
  USE calc_mod
  USE mat_funcs,only:CTANH,CSECH
!
  INTEGER J,K
  REAL (OP) :: MU0,AL2,ALDA,FN,SIGMU
  COMPLEX (OP) :: QLOLD,THETA,TNFNEW,ANSTM,IW
  COMPLEX (OP) :: QUOLD,CONU,CONL,THT
  COMPLEX (OP), DIMENSION(:), ALLOCATABLE :: QL,QU

  ALLOCATE(QL(NLL),QU(NLU+1))
  MU0=PI*4.D-07
  AL2=ALDA*ALDA
  IW=CMPLX(0.0,OMEGA)
  FN=CANIS(CNL)
  SIGMU=CCON(CNL)*MU0
  THETA=CSQRT(AL2*FN*FN+IW*SIGMU)
  QL(1)=THETA/CCON(CNL)
  K=1
  DO J=CNL-1,NLU+1,-1
     K=K+1
     QLOLD=QL(K-1)
     FN=CANIS(J)
     SIGMU=CCON(J)*MU0
     THETA=CSQRT(AL2*FN*FN+IW*SIGMU)
     TNFNEW=CTANH(THETA*CTHICK(J))
     QL(K)=THETA*(THETA*TNFNEW+CCON(J)*QLOLD)/(THETA+CCON(J)*TNFNEW*QLOLD)/CCON(J)
  END DO

! Now upper part part:
! This is for the air half-space
!  E=8.854187871D-12
!  THETA=CSQRT(AL2-IW*MU0*OMEGA*E)
!  QU(1)=-THETA/(IW*MU0*E)
! or
  SIGMU=1e-8*MU0
  THETA=CSQRT(AL2+IW*SIGMU)
  QU(1)=-THETA/1e-8

  DO J=1,NLU
     QUOLD=QU(J)
     FN=CANIS(J)
     SIGMU=CCON(J)*MU0
     THETA=CSQRT(AL2*FN*FN+IW*SIGMU)
     TNFNEW=CTANH(THETA*CTHICK(J))
     QU(J+1)=-THETA*(THETA*TNFNEW-CCON(J)*QUOLD)/(THETA-CCON(J)*TNFNEW*QUOLD)/CCON(J)
  END DO

  CONU=CCON(NLU)/CANIS(NLU)**2_OP
  CONL=CCON(NLU+1)/CANIS(NLU+1)**2_OP

! B_phi
  ANSTM = MU0*AL2/4/PI/(QL(NLL)-QU(NLU+1))*(1.0/CONU+1.0/CONL)

  IF (CZ<0) THEN
! Rx higher than Tx        
     DO J=NLU,NLRX,-1
        FN=CANIS(J)
        SIGMU=CCON(J)*MU0
        THETA=CSQRT(AL2*FN*FN+IW*SIGMU)
        THT=THETA*CMPLX(CTHICK(J),0.0)
        ANSTM=ANSTM*(CSECH(THT)*THETA)/(THETA-CCON(J)*QU(J)*CTANH(THT))
     ENDDO
     IF (COMP==2) ANSTM=QU(NLRX)*ANSTM/MU0
     IF (COMP==3) ANSTM=ANSTM/MU0*ALDA/CCON(NLRX)*CANIS(NLRX)**2_OP
  ELSE
! Rx deeper than Tx        
     DO J=NLU+1,NLRX-1
        FN=CANIS(J)
        SIGMU=CCON(J)*MU0
        THETA=CSQRT(AL2*FN*FN+IW*SIGMU)
        THT=THETA*CMPLX(CTHICK(J),0.0)
        ANSTM=ANSTM*(CSECH(THT)*THETA)/(THETA+CCON(J)*QL(CNL-J)*CTANH(THT))
     ENDDO
     IF (COMP==2) ANSTM=QL(CNL-NLRX+1)*ANSTM/MU0
!     IF (ALDA<0.0000011) WRITE(6,*) CNL-NLRX+1,QL(CNL-NLRX+1),ANSTM
     IF (COMP==3) ANSTM=ANSTM/MU0*ALDA/CCON(NLRX-1)*CANIS(NLRX)**2_OP
  ENDIF

  DEALLOCATE(QL,QU)
  RETURN
END SUBROUTINE EZML

