MODULE mat_funcs
  USE mar_para
  PUBLIC :: ANORM,INITM,MULT,TRMULT,DOT,DOT_NON,DIMULT,CTANH
  PUBLIC :: GASDEV,DPROD
CONTAINS 

  REAL (OP) FUNCTION ANORM(N,D)
! RETURNS THE SQUARE OF THE EUCLIDEAN NORM OF A VECTOR
    INTEGER N,I
    REAL (OP) :: D(N)                       
    ANORM = 0.0                                                              
    DO I = 1,N 
!     IF (D(I).LT.1E18) THEN 
       ANORM = ANORM + D(I)*D(I) 
!     ELSE 
!        Y = Y + 1E36 
!     ENDIF
    END DO
    RETURN                                                                   
  END FUNCTION ANORM
!                                                                    
!                                                                      
!****************************************************************************
  SUBROUTINE INITM(MD,M,A)                                               
! ZEROS A 2D MATRIX                                                          
    REAL (OP) :: A(MD,M)
    INTEGER M,MD,I,J                                                        
    DO I = 1,MD                                                            
       DO J = 1,M                                                          
          A(I,J) = 0.0
       END DO
    END DO
    RETURN                                                                   
  END SUBROUTINE INITM
!                                                                           
!                                                                           
!****************************************************************************
  SUBROUTINE MULT(MAD,NA,NB,A,B,C)                                  
! MULTIPLIES TWO 2D MATRICES
    INTEGER MAD,NA,NB,I,J,K
    REAL (OP) :: A(MAD,NA),B(NA,NB),C(MAD,NB),CIJ
    DO I = 1,MAD                                                           
       DO J = 1,NB                                                         
          CIJ = 0.0                                                          
          DO K = 1,NA                                                       
             CIJ = A(I,K)*B(K,J) + CIJ                                       
             C(I,J) = CIJ
          ENDDO
       END DO
    ENDDO
    RETURN                                                                   
  END SUBROUTINE MULT
!                                                                            
!                                                                            
!****************************************************************************
  SUBROUTINE DIMULT(MAD,NA,DIAG,A,B)                                    
    ! MULTIPLIES A 2D MATRIX BY A DIAGONAL MATRIX                          
    INTEGER MAD,NA,I,J
    REAL (OP) :: DIAG(MAD),A(MAD,NA),B(MAD,NA)                               
    DO J = 1,NA                                                           
       DO I = 1,MAD                                                         
          B(I,J) = DIAG(I)*A(I,J)
       END DO
    END DO
    RETURN                                                         
  END SUBROUTINE DIMULT
!                                                                
!                                                                  
!****************************************************************************
  SUBROUTINE TRMULT(MAD,NA,NB,A,B,C)                         
! MULTIPLIES THE TRANSPOSE OF A 2D MATRIX BY ANOTHER MATRIX         
    INTEGER :: MAD,NA,NB,I,J,K              
    REAL (OP) :: A(MAD,NA),B(MAD,NB),C(NA,NB),CIJ
    DO I = 1,NA                                                    
       DO J = 1,NB                                                         
          CIJ = 0.0                                                          
          DO K = 1,MAD                                                       
             CIJ = A(K,I)*B(K,J) + CIJ
          END DO
          C(I,J) = CIJ
       END DO
    END DO
    RETURN                                                         
  END SUBROUTINE TRMULT
!                                                                
!                                                                  
!****************************************************************************
  REAL (OP) FUNCTION DOT(A,B,N)                         
! Computes the dot-product of two vectors         
    INTEGER :: N,I             
    REAL (OP) :: A(N),B(N)
    DOT=0.0
    DO I=1,N
       DOT=DOT+A(I)*B(I)
    END DO
  END FUNCTION DOT

!****************************************************************************
  REAL (OP) FUNCTION DOT_NON(N,X1,INC1,X2,INC2)
!  COMPUTES DOT PRODUCT OF TWO D.P. VECTORS WITH NONUNIT
!  INCREMENTING ALLOWED. REPLACEMENT FOR BLAS SUBROUTINE SDOT.
    INTEGER :: N,K,INC1,INC2,I
    REAL (OP) :: X1(N*INC1),X2(N*INC2)
    IF(INC2.GT.0)THEN
       K=1
    ELSE
       K=N*ABS(INC2)
    ENDIF
    DOT_NON=0.0
    IF(INC1.GT.0)THEN
       DO I=1,N,INC1
          DOT_NON=DOT_NON+X1(I)*X2(K)
          K=K+INC2
       END DO
    ELSE
       DO I=N,1,INC1
          DOT_NON=DOT_NON+X1(I)*X2(K)
          K=K+INC2
       END DO
    ENDIF
    RETURN
  END FUNCTION DOT_NON

!****************************************************************************
  REAL (OP) FUNCTION DPROD (N,N1,N2,B,C)
!
!     DOUBLE PRECISION INNER PRODUCT ROUTINE
!     DPROD = B*C
!       B,C = VECTORS (CAN BE ROWS OR COLS OF ARRAYS)
!         N = LENGTH OF VECTORS
!     N1,N2 = INCREMENT FOR B,C
!           = 1 IF COL OF ARRAY
!           = COL LENGTH (I.E. NO. OF ROWS) IF ROW OF ARRAY
!
!     NOTE ... CHECK IF MACHINE STORES BY ROWS OR COLUMNS
!
!     DPROD MUST BE DECLARED EXTERNAL BY ANY ROUTINE USING
!     IT AS THERE IS A STANDARD INTRINSIC FUNCTION WITH
!     THE SAME NAME. (IF OMITTED YOU GET COMPILATION
!     WARNINGS)
!
!     CALLED BY ESVD,FSATI,NLSQ2,SOLVE2,STANV,STATS
!
    REAL (OP) :: Z1,Z2,Z3
    REAL (OP),DIMENSION(N*N1):: B
    REAL (OP),DIMENSION(N*N2):: C
    INTEGER NA,NB,N,N1,N2,I
!
    Z1=0.0
    IF (N>0) THEN
       NA=1
       NB=1
       DO I=1,N
          Z2=B(NA)
          Z3=C(NB)
          Z1=Z1+Z2*Z3
          NA=NA+N1
          NB=NB+N2
       END DO
    END IF
    DPROD=Z1
    RETURN
!
!     END OF DPROD
!
  END FUNCTION DPROD

!
!
!****************************************************************************

  COMPLEX (OP) FUNCTION CSECH(X)
    COMPLEX (OP) :: X
    
    CSECH=2.*CEXP(-X)/(1+CEXP(-2.*X))

    RETURN
  END FUNCTION CSECH
!
!
!****************************************************************************

  COMPLEX (OP) FUNCTION CTANH(X)
    COMPLEX (OP) :: X

    CTANH=(1.0_OP-EXP(-2.0_OP*X))/(1.0_OP+EXP(-2.0_OP*X))
    RETURN
  END FUNCTION CTANH

!
!
!*********************************************************************
   REAL (OP) FUNCTION GASDEV(GSET)
!------------------------------------------------------------------------------
! FILE:NOISE.FOR
! NAME:GASDEV
!	PURPOSE    : Returns a normally distributed deviate with zero 
!		     mean and unit variace, using RAN(IDUM) as source
!		     of uniform distributed deviates
!
!	USAGE      : deviate=GASDEV(idum)
!
!	       IDUM: INTEGER STARTING VALUE
!
!	   DESCRIPTION -
!			The program uses the Box-Muller transformation
!		It is taken from ' NUMERICAL RECIPES',
!		'The art of scientific computing',p.203
!		William H.Press et al.,Cambridge 1986
!
!------------------------------------------------------------------------------
     INTEGER ISET
     DATA ISET/0/
     REAL GSET,DUM1,DUM2,R,V1,V2,FAC
     SAVE
     IF(ISET.EQ.0) THEN
1       DUM1=RAND(0)
        DUM2=RAND(0)
        V1=2.*DUM1-1.
        V2=2.*DUM2-1.
        R=V1**2+V2**2
        IF(R.GE.1) GOTO 1
        FAC=SQRT(-2.*ALOG(R)/R)
        GSET=V1*FAC                
        GASDEV=V2*FAC
        ISET=1
     ELSE
        GASDEV=GSET
        ISET=0
     ENDIF
     RETURN
   END FUNCTION GASDEV

END MODULE mat_funcs
