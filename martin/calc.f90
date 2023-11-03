SUBROUTINE CALC
  USE mar_data
  USE calc_mod
  INTEGER I,J,NOD,SIGN,T2,SPL_N,NS,NSPL_N,NPERIOD
  REAL (OP) :: EARLIEST,LATEST,TEST,SPV,DEL,AT_OFFSET,PREV,CAL_AT
  REAL (OP) :: Y2SP,YSP,INTERVALL,HIGH_RES,R,LOW_RES
  LOGICAL ONEMORE,CHANGE
  REAL (OP), DIMENSION(:), ALLOCATABLE :: YAB1,PK1,QK1,AS1,BS1,CS1,DS1,HELP1
  REAL (OP), DIMENSION(:), ALLOCATABLE :: YAB2,PK2,QK2,AS2,BS2,CS2,DS2,HELP2
  REAL (OP), DIMENSION(:), ALLOCATABLE :: SPL_T,SPL_C,HELPT
  REAL (OP) :: INT_CHECK,BUF

!  IF (NDAT<1) THEN
!     WRITE(6,*) '<E> No data set in memory'
!     RETURN
!  ENDIF
!  IF (NLAY<1) THEN
!     WRITE(6,*) '<E> No model in memory'
!     RETURN
!  ENDIF
  NT_F=0
  EARLIEST=1e20
  LATEST=-2
  HIGH_RES=W_RES
  LOW_RES=W_RES
  DO I=1,NLAY
     HIGH_RES=MAX(RES(I),HIGH_RES,maxres) ! added maxres/2 to adjust to data input and maxres can be quite large, but actually the time vector won't need to start early if the maximum resistivity is at greater depth
     LOW_RES=MIN(RES(I),LOW_RES,minres)
  ENDDO
  DO I=1,NDAT
     R=SQRT(X(I)*X(I)+Y(I)*Y(I))
     DO J=1,NT(I)
        IF (SYS_RESP(I)) THEN
           LATEST=MAX(MAX(LATEST,(R/1261)**2/LOW_RES*25),TIME(I,J)*20)
	       ! EARLIEST=MIN(EARLIEST,(R/1261)**2/4/50) ! MAXRES=4 FOR THE NORTHSEA PROBLEM 
            EARLIEST=MIN(EARLIEST,(R/1261)**2/HIGH_RES/50) !ORIGINAL
        ELSE
           LATEST=MAX(LATEST,TIME(I,J))
           EARLIEST=MIN(EARLIEST,TIME(I,J))
        END IF
     END DO
  END DO
  I=1
  DO 
     TEST=EARLIEST*SQRT10**(I-1)
     IF (TEST > LATEST) EXIT
     I=I+1
  END DO

  IF (ALLOCATED(FTIME)) DEALLOCATE(FTIME)
  ALLOCATE(FTIME(I))
  NT_F=MAX(NT_F,I)
  I=1
  DO 
     FTIME(I)=EARLIEST*SQRT10**(I-1)
     IF (FTIME(I) > LATEST) EXIT
     I=I+1
  END DO

  IF (ALLOCATED(EOUT)) DEALLOCATE(EOUT)
  IF (LFD) THEN
     ALLOCATE(EOUT(NDAT,NT_F*2))
  ELSE
     ALLOCATE(EOUT(NDAT,NT_F))
  ENDIF
  CALL INLINE_FD_FHT(EARLIEST)
  CALL HR_FD_FHT(EARLIEST)
  CALL EX_FD_FHT(EARLIEST)
  CALL EZ_TX_FHT(EARLIEST)

  IF (LFD) THEN
     ALLOCATE(YAB1(NT_F),PK1(NT_F-1),QK1(NT_F-1),AS1(NT_F-1),BS1(NT_F-1),&
          CS1(NT_F-1),DS1(NT_F-1),HELP1(NT_F))
     DO NOD=1,NDAT
        DO I=1,NT_F-1
           PK1(I)=10.0
           QK1(I)=10.0
        END DO
        HELP1=EOUT(NOD,1:NT_F)
        YAB1(1)=0
        YAB1(NT_F)=(HELP1(NT_F)-HELP1(NT_F-1))/(FTIME(NT_F)-FTIME(NT_F-1))
        CALL RASPL(NT_F,FTIME,HELP1,PK1,QK1,YAB1,AS1,BS1,CS1,DS1)
        DO I=1,NT(NOD)/2
           CALL YSPLINE(NT_F,FTIME,TIME(NOD,I),AS1,BS1,CS1,DS1,YSP,PK1,QK1)
           CVOLT(NOD,I)=YSP*CALF(NOD)*CUR*DL
        END DO
        HELP1=EOUT(NOD,NT_F+1:NT_F*2)
        YAB1(1)=0
        YAB1(NT_F)=(HELP1(NT_F)-HELP1(NT_F-1))/(FTIME(NT_F)-FTIME(NT_F-1))
        CALL RASPL(NT_F,FTIME,HELP1,PK1,QK1,YAB1,AS1,BS1,CS1,DS1)
        DO I=NT(NOD)/2+1,NT(NOD)
           CALL YSPLINE(NT_F,FTIME,TIME(NOD,I),AS1,BS1,CS1,DS1,YSP,PK1,QK1)
           CVOLT(NOD,I)=-YSP*CALF(NOD)*CUR*DL
        END DO
     END DO
 ELSE
     ALLOCATE(YAB1(NT_F),PK1(NT_F-1),QK1(NT_F-1),AS1(NT_F-1),BS1(NT_F-1),&
          CS1(NT_F-1),DS1(NT_F-1),HELP1(NT_F))
     DO NOD=1,NDAT
! This is for splining the forward data to certain points
        DO I=1,NT_F-1
           PK1(I)=10.0
           QK1(I)=10.0
        END DO
        HELP1=EOUT(NOD,:)
        YAB1(1)=0
        YAB1(NT_F)=(HELP1(NT_F)-HELP1(NT_F-1))/(FTIME(NT_F)-FTIME(NT_F-1))
        CALL RASPL(NT_F,FTIME,HELP1,PK1,QK1,YAB1,AS1,BS1,CS1,DS1)

        IF (SYS_RESP(NOD)) THEN
           NS=SYS_N(NOD)
           ALLOCATE(SPL_T(NT_F+NS),SPL_C(NT_F+NS))

!This is for splining the system response
           ALLOCATE(YAB2(NS),PK2(NS-1),QK2(NS-1))
           ALLOCATE(AS2(NS-1),BS2(NS-1),CS2(NS-1))
           ALLOCATE(DS2(NS-1),HELP2(NS),HELPT(NS))
           DO I=1,NS-1
              PK2(I)=10.0
              QK2(I)=10.0
           END DO
           HELP2=SYS_D(NOD,:)
           HELPT=SYS_T(NOD,:)
           YAB2(1)=0
           YAB2(NS)=(HELP2(NS)-HELP2(NS-1))/(HELPT(NS)-HELPT(NS-1))
           CALL RASPL(NS,HELPT,HELP2,PK2,QK2,YAB2,AS2,BS2,CS2,DS2)
!Ok, now convolve each data point with the system response
           DO I=1,NT(NOD)
              SPV=0.0
              AT_OFFSET=0
              ONEMORE=.TRUE.
              SIGN=1
              T2=0
              NPERIOD=0
              DO WHILE (ONEMORE)
                 NPERIOD=NPERIOD+1
                 IF (NPERIOD.EQ.(NPERIOD/2)*2) PREV=SPV
!First, define the positions, where both functions have to be evaluated
                 SPL_N=0
                 DO J=1,NS
                    CAL_AT=TIME(NOD,I)+AT_OFFSET-SYS_T(NOD,J)
                    IF (CAL_AT<0.0.OR.CAL_AT>FTIME(NT_F)) CYCLE
                    SPL_N=SPL_N+1
                    SPL_T(SPL_N)=CAL_AT
                    SPL_C(SPL_N)=SYS_T(NOD,J)
                 ENDDO
                 DO J=1,NT_F
                    IF (FTIME(J)>=TIME(NOD,I)+AT_OFFSET-SYS_PER(NOD)*0.001.AND. &
                         FTIME(J)<=TIME(NOD,I)+AT_OFFSET) THEN
                       SPL_N=SPL_N+1
                       SPL_T(SPL_N)=FTIME(J)
                       SPL_C(SPL_N)=TIME(NOD,I)+AT_OFFSET-FTIME(J)
                    END IF
                 ENDDO
!Now sort all times
                 CHANGE=.TRUE.
                 DO WHILE (CHANGE)
                    CHANGE=.FALSE.
                    DO J=1,SPL_N-1
                       IF (SPL_C(J)>SPL_C(J+1)) THEN
                          CHANGE=.TRUE.
                          BUF=SPL_C(J)
                          SPL_C(J)=SPL_C(J+1)
                          SPL_C(J+1)=BUF
                          BUF=SPL_T(J)
                          SPL_T(J)=SPL_T(J+1)
                          SPL_T(J+1)=BUF
                       END IF
                    END DO
                 END DO
!Delete double times
                 IF (SPL_N>0) THEN
                    NSPL_N=1
                    DO J=2,SPL_N
                       IF (SPL_C(NSPL_N)/=SPL_C(J)) THEN
                          NSPL_N=NSPL_N+1
                          SPL_C(NSPL_N)=SPL_C(J)
                          SPL_T(NSPL_N)=SPL_T(J)
                       END IF
                    END DO
                    SPL_N=NSPL_N
!                 WRITE(6,*) "Data involved from sys_resp: ",SPL_C(1),SPL_C(SPL_N)
!                 WRITE(6,*) "Data involved from forward:  ",SPL_T(1),SPL_T(SPL_N),FTIME(NT_F)
                 END IF
!Ready, now start to spline
                 IF (NT_F<2) EXIT
                 INT_CHECK=0.0
                 DO J=1,SPL_N
                    IF (SPL_T(J)<FTIME(1)) THEN
                       YSP=HELP1(1)
                    ELSEIF (SPL_T(J)>FTIME(NT_F)) THEN
                       CYCLE
                    ELSE
                       CALL YSPLINE(NT_F,FTIME,SPL_T(J),AS1,BS1,CS1,DS1,YSP,PK1,QK1)
                    ENDIF
                    CALL YSPLINE(NS,HELPT,SPL_C(J),AS2,BS2,CS2,DS2,Y2SP,PK2,QK2)
                    IF (J==1) THEN
                       IF (SPL_N>1) THEN
                          INTERVALL=SPL_C(2)/2.0
                       ELSE
                          INTERVALL=(TIME(NOD,I+1)-TIME(NOD,I))/2.0
                       END IF
                    ELSEIF (J==SPL_N) THEN
                       INTERVALL=MIN(SYS_PER(NOD)*0.001-(SPL_C(J)+ &
                            SPL_C(J-1))/2.0,(SPL_T(J)+SPL_T(J-1))/2.0)
                    ELSE
                       INTERVALL=(SPL_C(J+1)-SPL_C(J-1))/2.0
                    END IF
                    INT_CHECK=INT_CHECK+INTERVALL*1000.0/SYS_PER(NOD)
!                   SPV=SPV+SIGN*YSP*Y2SP/SYS_PER(NOD)*INTERVALL*1000.0
                    SPV=SPV+SIGN*YSP*Y2SP*INTERVALL
!                   SPV=SPV+SIGN*YSP*Y2SP
                 END DO
                 T2=T2+1
                 DEL=ABS((SPV-PREV)/SPV)
!                ONEMORE=(NPERIOD.NE.4)
                 ONEMORE=.NOT.(DEL<0.01.AND.NPERIOD>2.AND. &
                      NPERIOD.EQ.(NPERIOD/2)*2)
                 AT_OFFSET=AT_OFFSET+SYS_PER(NOD)*0.001
                 IF (SYS_TYPE(NOD).EQ.1) SIGN=-1*SIGN
              END DO
              CVOLT(NOD,I)=SPV*CALF(NOD)*CUR*DL
           END DO
           DEALLOCATE(YAB2,PK2,QK2,AS2,BS2,CS2,DS2,HELP2)
           DEALLOCATE(SPL_T,SPL_C,HELPT)
        ELSE
           DO I=1,NT(NOD)
              CALL YSPLINE(NT_F,FTIME,TIME(NOD,I),AS1,BS1,CS1,DS1,YSP,PK1,QK1)
              CVOLT(NOD,I)=YSP*CALF(NOD)*CUR*DL
           END DO
        ENDIF
     END DO
     DEALLOCATE(YAB1,PK1,QK1,AS1,BS1,CS1,DS1,HELP1)
  ENDIF
  DEALLOCATE (EOUT)

  IF (ALLOCATED(FTIME)) DEALLOCATE(FTIME)
END SUBROUTINE CALC

