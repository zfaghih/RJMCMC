SUBROUTINE RASPL(N,X,Y,P,Q,Y1,A,B,C,D)
  USE mar_para
! pc
! Calculate generalized cubic spline coefficients. The coefficients
! A,B,C,D calculated by this subroutine may be used in the spline
! routine YSPLINE.
! INPUT:
!	N	number of input points in X,Y
!	X(N)	arguments of input function
!	Y(N)	function values
!	P(N-1)	(smoothness ??) - tested with P(I)=20.,I=1,N
!	Q(N-1)	(smoothness ??) - tested with Q(I)=20.,I=1,N
!	Y1(N)	first derivative, on input required at first point
!		Y1(1) and at N'th point Y1(N) only
! OUTPUT:
!	A(N-1),B(N-1),C(N-1),D(N-1)	spline coefficients between
!		any two input points
!
! Subroutine taken from "spline-algorithmen, H. Speth, R.Oldenbourg
! Verlag, Muenchen, 1987", pp. 103
  INTEGER N,N1,K,J2,N2,J1
  REAL (OP) :: X(N),Y(N),P(N-1),Q(N-1),Y1(N),A(N-1),B(N-1),C(N-1),D(N-1)
  REAL (OP) :: PP,QQ,PP2,QQ2,P22,Q22,H,H2,R2,HQ,HP,Z,P21,QQ1,H1,R1,P2,Q2
  N1=N-1
  C(1)=0.
  D(1)=0.
  DO K=1,N1
     J2=K+1
     PP=P(K)
     QQ=Q(K)
     PP2=PP*(PP+3.)+3.
     QQ2=QQ*(QQ+3.)+3.
     P22=2.+PP
     Q22=2.+QQ
     A(K)=X(J2)-X(K)
     H=1./A(K)
     B(K)=1./(P22*Q22-1.)
     H2=H*B(K)
     R2=H*H2*(Y(J2)-Y(K))
     IF(K.NE.1) THEN
        HQ=H1*QQ1
        HP=H2*PP2
        Z=1./(HQ*(P21-C(J1))+HP*Q22)
        C(K)=Z*HP
        H=R1*QQ1*(1.+P21)+R2*PP2*(1.+Q22)
        IF(K.EQ.2) H=H-HQ*Y1(1)
        IF(K.EQ.N1) H=H-HP*Y1(N)
        D(K)=Z*(H-HQ*D(J1))
     END IF
     J1=K
     P21=P22
     QQ1=QQ2
     H1=H2
     R1=R2
  END DO
  Y1(N1)=D(N1)
  IF(N1.GT.2) THEN
     N2=N1-1
     DO J1=2,N2
        K=N-J1
        Y1(K)=D(K)-C(K)*Y1(K+1)
     ENDDO
  END IF
  DO K=1,N1
     J2=K+1
     H=B(K)*(Y(J2)-Y(K))
     Z=B(K)*A(K)
     P2=2.+P(K)
     Q2=2.+Q(K)
     C(K)=(1.+Q2)*H-Z*(Y1(J2)+Q2*Y1(K))
     D(K)=-(1.+P2)*H+Z*(P2*Y1(J2)+Y1(K))
     A(K)=Y(K)-C(K)
     B(K)=Y(J2)-D(K)
  ENDDO
  RETURN
END SUBROUTINE RASPL

SUBROUTINE YSPLINE(N,X,XWERT,A,B,C,D,YWERT,P,Q)
  USE mar_para
! pc
! Calculate the function value YWERT at a given variable value
! XWERT by cubic spline interpolation.
! INPUT:
!	N	number of input points in X,Y
!	X(N)	arguments of input function
!	Y(N)	known function values at X(N)
!	P(N-1)	(smoothness ??) - tested with P(I)=20.,I=1,N
!	Q(N-1)	(smoothness ??) - tested with Q(I)=20.,I=1,N
!	A(N-1),B(N-1),C(N-1),D(N-1)	spline coefficients between
!		any two input points as calculated by SUBROUTINE 
!		RASPL1.
!	XWERT	variable value where the function value is desired
! OUTPUT:
!	YWERT	function value at XWERT
! Subroutine coded after from "spline-algorithmen, H. Speth, R.Oldenbourg
! Verlag, Muenchen, 1987", pp. 101, eq. (7.64).
  INTEGER :: N,I
  REAL (OP) :: X(N),P(N-1),Q(N-1),A(N-1),B(N-1),C(N-1),D(N-1)
  REAL (OP) :: T,XWERT,YWERT,U
  
! find proper interval in x:
  I=1
11 IF(XWERT.Lt.X(I+1)) GOTO 10
  I=I+1
  IF(I.EQ.N)          GOTO 8
  GOTO 11
8 I=N-1
! calculate new y-value:
10 T=(XWERT-X(I))/(X(I+1)-X(I))
  U=1.-T
  YWERT=A(I)*U+B(I)*T+C(I)*U*U*U/(P(I)*T+1.)+D(I)*T*T*T/(Q(I)*U+1.)
  RETURN
END SUBROUTINE YSPLINE
