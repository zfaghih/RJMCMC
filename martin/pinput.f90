SUBROUTINE CIN(PROMPT,C)
!------------------------------------------------------------------------
! I/O TASK:		TEXT STRING passed on a CHARACTER*(*) variable
!------------------------------------------------------------------------
! VARIABLES:
!	PROMPT	character string prompting the user for input.
!		The user must NOT enclose the character string
!		in quotes ('xxx').
!	C	character variable, on call of subroutine: default
!		value, on return: default or new string
!------------------------------------------------------------------------
  CHARACTER*(*) PROMPT,C
  CHARACTER ANS*80
  WRITE(6,'(1X,A,1X,A,A,A,$)')PROMPT,'(',C,') '
  READ(5,'(A)')ANS
  IF(ANS.NE.'  ')C=ANS
  RETURN
END SUBROUTINE CIN
!=======================================================================
SUBROUTINE IIN(PROMPT,I)
!------------------------------------------------------------------------
! I/O TASK:		INTEGER passed on an INTEGER variable
!------------------------------------------------------------------------
! VARIABLES:
!	PROMPT	character string prompting the user for input
!	I	integer variable, on call of subroutine: default
!		value, on return: default or new value
!------------------------------------------------------------------------
  CHARACTER*(*) PROMPT
  CHARACTER ANS*80
  INTEGER I
  WRITE(6,'(1X,A,A,I9,A,$)')PROMPT,'(',I,') '
  READ(5,'(A)')ANS
  IF(ANS.NE.'  ')READ(ANS,'(i9)')I
  RETURN
END SUBROUTINE IIN
!=======================================================================
SUBROUTINE RIN(PROMPT,R)
!------------------------------------------------------------------------
! I/O TASK:		REAL*8 passed on an REAL*8 variable
!------------------------------------------------------------------------
! VARIABLES:	same as RIN but for real*8
!------------------------------------------------------------------------
  USE mar_para
  CHARACTER*(*) PROMPT
  REAL (OP):: R
  CHARACTER ANS*80
  WRITE(6,'(1X,A,A,G9.3,A,$)')PROMPT,'(',R,') '
  READ(5,'(A)')ANS
  IF(ANS.NE.'  ')READ(ANS,*)R
  RETURN
END SUBROUTINE RIN
!=======================================================================
SUBROUTINE IINN(PROMPT,I,N)
!-----------------------------------------------------------------------
! PURPOSE:
!	Same as IIN but for an input array.
! PARAMETERS:
!	PROMPT	CHARACTER*(*)	prompt string
!	I(N)	INTEGER		I/O array
!	N	INTEGER	        number of elements in I
!-----------------------------------------------------------------------
  CHARACTER*(*) PROMPT
  CHARACTER ANS*180
  INTEGER N,K,K1,I(N)
  IF (N<1) RETURN
  WRITE(6,'(1X,A,$)')PROMPT
  K1=1
10 CONTINUE
  READ(5,'(A)')ANS
  IF(ANS.NE.' ')READ(ANS,'(i9)',END=20)(I(K),K=K1,N)
  GOTO 30
20 K1=K
  GOTO 10
30 CONTINUE
  RETURN
END SUBROUTINE IINN
!=======================================================================
SUBROUTINE RINN(PROMPT,R,N)
!-----------------------------------------------------------------------
! PURPOSE:
!	Same as RIN but for an input array.
! PARAMETERS:
!	PROMPT	CHARACTER*(*)	prompt string
!	R(N)	REAL		I/O array
!	N	INTEGER	        number of elements in R
!-----------------------------------------------------------------------
  USE mar_para
  CHARACTER*(*) PROMPT
  CHARACTER :: ANS*180
  INTEGER :: K1,K,N
  REAL (OP) :: R(N)
  IF (N<1) RETURN
  WRITE(6,'(1X,A,$)')PROMPT
  K1=1
10 CONTINUE
  READ(5,'(A)')ANS
  IF(ANS.NE.' ')READ(ANS,*,END=20)(R(K),K=K1,N)
  GOTO 30
20 K1=K
  GOTO 10
30 CONTINUE
  RETURN
END SUBROUTINE RINN
!=======================================================================
LOGICAL FUNCTION LIN(PROMPT,LDEF)
!------------------------------------------------------------------------
!  PURPOSE:
!	Function LIN prompts for a YES or NO and returns then a .FALSE. 
!	for NO and a .TRUE. for YES.
!  VARIABLES:
!	PROMPT	CHARACTER*(*)	prompt text
!	LDEF	LOGICAL		default value
!------------------------------------------------------------------------
  CHARACTER*(*) PROMPT
  CHARACTER C*2
  LOGICAL LDEF
! define the string default according to the logical default:
  IF(LDEF)THEN
     C='YE'
  ELSE
     C='NO'
  ENDIF
! prompt for string:
1 CALL CIN(PROMPT,C)
! validity check:
  IF(C.NE.'NO'.AND.C.NE.'no'.AND.C.NE.'YE'.AND.C.NE.'ye') THEN
     WRITE(6,'(1X,A$)')' Answer YES or NO: '
     GOTO 1
  ENDIF
! evaluate function:
  LIN=(C.EQ.'YE'.OR.C.EQ.'ye')
  RETURN
END FUNCTION LIN

SUBROUTINE LCCONV(WORD,LEN)
!	This routine converts lower case in upper case characters.
!
!	parameters: WORD,LEN
!		WORD is a byte string of length LEN. It can be of
!		logical, real, or integer type
  INTEGER LEN,I
  CHARACTER*1 WORD(LEN)
  DO I=1,LEN
     IF(ICHAR(WORD(I)).GE.97.AND.ICHAR(WORD(I)).LE.122) WORD(I)=CHAR(ICHAR(WORD(I))-32)
  ENDDO
END SUBROUTINE LCCONV

