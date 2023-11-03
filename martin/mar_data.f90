MODULE mar_data
  USE mar_para

! Model parameters
  INTEGER, PUBLIC :: NLAY               !Number of layers
  REAL (OP), PUBLIC,DIMENSION(:),ALLOCATABLE :: TNESS,RES,ANIS
  REAL (OP), PUBLIC :: W_DEPTH,W_RES       !Water depth and resistivity
  REAL (OP), PUBLIC :: minres, maxres ! minimum and maximum resistivity
  LOGICAL, PUBLIC, DIMENSION(:), ALLOCATABLE :: LFIXA,LFIXR,LFIXT ! Fixes for the layer parameters
  LOGICAL, PUBLIC :: LFIXWD,LFIXWR ! Fixes for water depth and resistivity
! Data set parameters
  INTEGER, PUBLIC :: NDAT,TNT,MAXNT,ACT_DAT    !Number of data sets, total number of time points, maximum number of data points, active data set
  INTEGER, PUBLIC,DIMENSION(:),ALLOCATABLE :: NT  !Number of time points
  CHARACTER (LEN=15), PUBLIC :: TX_DESIG
  REAL (OP), PUBLIC :: TX_X,TX_Y,TX_Z,DL,CUR,TCHI!Tx-parameters, total misfit
  LOGICAL, PUBLIC :: LFD ! Frequency domain?
  REAL (OP), PUBLIC,DIMENSION(:,:),ALLOCATABLE :: TIME,VOLT,ERROR,CVOLT,WTS 
  REAL (OP), PUBLIC,DIMENSION(:),ALLOCATABLE :: X,Y,Z,ANG1,ANG2 !Receiver geometry
  CHARACTER (LEN=4), PUBLIC,DIMENSION(:),ALLOCATABLE :: TYPE !Receiver type
  REAL (OP), PUBLIC, DIMENSION(:),ALLOCATABLE :: CALF,DELAY,CHI !Calibration factor, delay, misfit
  LOGICAL, PUBLIC, DIMENSION(:), ALLOCATABLE :: LFCAL,SYS_RESP !Fixed factor, system-response
  LOGICAL, PUBLIC, DIMENSION(:), ALLOCATABLE :: IMP,STEP_OFF,DERIVATE !Impulse response, step_off response, derivative of curve
  CHARACTER (LEN=15), PUBLIC,DIMENSION(:),ALLOCATABLE :: SYS_NAME

  REAL (OP), PUBLIC,DIMENSION(:,:), ALLOCATABLE :: SYS_T,SYS_D
  REAL (OP), PUBLIC,DIMENSION(:), ALLOCATABLE :: SYS_PER
  INTEGER, PUBLIC,DIMENSION(:), ALLOCATABLE :: SYS_N,SYS_TYPE
  REAL (OP), PUBLIC :: ERR_MUL ! Normalization if error multiplier

! Unit numbers
  LOGICAL UN_NO(100)

! This is the reference model as used by EQUIVALENCE
  INTEGER :: REF_NL
  REAL (OP), PUBLIC,DIMENSION(:),ALLOCATABLE :: REF_T,REF_A,REF_C,REF_R
  REAL (OP), PUBLIC :: REF_W,REF_D,REF_CHI
  LOGICAL, PUBLIC, DIMENSION(:), ALLOCATABLE :: LRFIXA,LRFIXR,LRFIXT 
END MODULE mar_data
