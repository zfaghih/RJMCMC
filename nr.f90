MODULE nr
       INTERFACE pythag
          FUNCTION pythag_dp(a,b)
          USE nrtype
          REAL(DP), INTENT(IN) :: a,b
          REAL(DP) :: pythag_dp
          END FUNCTION pythag_dp
!BL
          FUNCTION pythag_sp(a,b)
          USE nrtype
          REAL(SP), INTENT(IN) :: a,b
          REAL(SP) :: pythag_sp
          END FUNCTION pythag_sp
       END INTERFACE
       INTERFACE svbksb
          SUBROUTINE svbksb_dp(u,w,v,b,x)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,v
          REAL(DP), DIMENSION(:), INTENT(IN) :: w,b
          REAL(DP), DIMENSION(:), INTENT(OUT) :: x
          END SUBROUTINE svbksb_dp
!BL
          SUBROUTINE svbksb_sp(u,w,v,b,x)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(IN) :: u,v
          REAL(SP), DIMENSION(:), INTENT(IN) :: w,b
          REAL(SP), DIMENSION(:), INTENT(OUT) :: x
          END SUBROUTINE svbksb_sp
       END INTERFACE
       INTERFACE svdcmp
          SUBROUTINE svdcmp_dp(a,w,v)
          USE nrtype
          REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
          REAL(DP), DIMENSION(:), INTENT(OUT) :: w
          REAL(DP), DIMENSION(:,:), INTENT(OUT) :: v
          END SUBROUTINE svdcmp_dp
!BL
          SUBROUTINE svdcmp_sp(a,w,v)
          USE nrtype
          REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
          REAL(SP), DIMENSION(:), INTENT(OUT) :: w
          REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
          END SUBROUTINE svdcmp_sp
       END INTERFACE
END MODULE nr
