!==============================================================================
!
!  Reversible Jump MCMC Sampling with parallel tempering for Controlled Source EM
!
!------------------------------------------------------------------------------
!
!  Jan Dettmer, University of Victoria, February 11 2013
!  jand@uvic.ca                       (250) 472 4026
!  http://web.uvic.ca~/jand/
!  Last change: April 13 2013
!
!  Based on Green 1995, Malinverno 2002, Bodin Sambridge 2009, 
!  Agostinetti Malinverno 2010
!
!==============================================================================

MODULE DATA_TYPE
   IMPLICIT NONE
   INTEGER(KIND=4), PARAMETER :: IB=4, RP=KIND(0.0D0), SP=KIND(0.0)
   REAL(KIND=RP),   PARAMETER :: PI2  = 3.141592653589793238462643383279502884197_RP
END MODULE DATA_TYPE

!==============================================================================
MODULE RJMCMC_COM
   USE MPI
   USE DATA_TYPE
   IMPLICIT NONE

!!
!! General switches
!!
      INTEGER(KIND=IB), PARAMETER :: ISETSEED   = 0     ! Fix the random seed 
      INTEGER(KIND=IB), PARAMETER :: ISINPRIOR  = 0     ! Use sine prior with depth
      INTEGER(KIND=IB) :: IMAP           ! WRITE REPLICA AND EXIT
      INTEGER(KIND=IB) :: SAMPLEPPD           ! WRITE tones of REPLICA AND EXIT
!!
!! Sampling move related switches
!!
   INTEGER(KIND=IB), PARAMETER :: ILINROT    = 0     ! 1 = turn on linear rotation proposal
   INTEGER(KIND=IB), PARAMETER :: IDELAY     = 0     ! 1 = turn on delayed rejection
   INTEGER(KIND=IB), PARAMETER :: IEXCHANGE  = 1     ! 1 = turn on exchange moves (parallel tempering)
   INTEGER(KIND=IB), PARAMETER :: ICOOL      = 0     ! Do cooling during burn-in
!!
!!  AUV DATA RJMCMC trial ping
!!
   INTEGER(KIND=IB) :: ICOV                   ! 0 = Sample implicit over sigma
                                              ! 1 = Sample over sigma
                                              ! 2 = Covariance matrix
						    ! 3 = Covariance matrix with variable amplitude
						    ! 4 = Variance only matrix
						    ! 5 = Variance only matrix with factor
                                              ! 6 = fixed sigma
!
   INTEGER(KIND=IB) :: IAR                    ! 1 = Use Autoregressive error model
   INTEGER(KIND=IB) :: ICF                    ! 1 = Use scaling factor on CSEM data
                                              ! 4 = Use scaling factor only on first receiver
   INTEGER(KIND=IB) :: IDT                    ! 1 = invert for time shift
                                              ! 2 = invert for one time shift for all receivers
   INTEGER(KIND=IB) :: IWRES                  ! 1 = invert for water resistivity
   INTEGER(KIND=IB) :: IWDEPTH                ! 0 = invert for water depth
   INTEGER(KIND=IB),PARAMETER :: IdB   = 0    ! 1 = Carry out computation in dB
 
!   INTEGER(KIND=IB),PARAMETER :: NTIME  =64             ! Maximum number of time samples   ! Romina
!   INTEGER(KIND=IB),PARAMETER :: NLMX   =15                ! Max number of layers   
!   INTEGER(KIND=IB),PARAMETER :: NRF=4 ! number of receivers ! Romina     

!  INTEGER(KIND=IB),PARAMETER :: NTIME  =70                ! Maximum number of time samples   ! Romina
!  INTEGER(KIND=IB),PARAMETER :: NLMX   =10                ! Max number of layers   
!  INTEGER(KIND=IB),PARAMETER :: NRF=2 ! number of receivers ! Romina

! Black Sea
!  INTEGER(KIND=IB),PARAMETER :: NTIME  =36                ! Maximum number of time samples   ! Romina
!  INTEGER(KIND=IB),PARAMETER :: NLMX   =15                ! Max number of layers   
!  INTEGER(KIND=IB),PARAMETER :: NRF=3 ! number of receivers ! Romina
! New Zealand
  INTEGER(KIND=IB),PARAMETER :: NTIME  = 38                ! Maximum number of time samples   ! Romina
  INTEGER(KIND=IB),PARAMETER :: NLMX   = 11                ! Max number of layers   
  INTEGER(KIND=IB),PARAMETER :: NRF	   = 4 ! number of receivers ! Romina

!   INTEGER(KIND=IB),PARAMETER :: NTIME  = 45                ! Maximum number of time samples   ! Romina
!   INTEGER(KIND=IB),PARAMETER :: NLMX   = 15                ! Max number of layers   
!   INTEGER(KIND=IB),PARAMETER :: NRF = 3 ! number of receivers ! Romina     
 
   INTEGER(KIND=IB),PARAMETER :: NPL   = 2    ! Number parameters per layer

   CHARACTER(len=64) :: filebasefile      = 'filebase.txt'
!!
!!  Prior variables and good seeding model
!!
   REAL(KIND=RP),DIMENSION(NPL):: minlim   = 0._RP
   REAL(KIND=RP),DIMENSION(NPL):: maxlim   = 0._RP
   !REAL(KIND=RP) :: minres ! defined as public variable in mar_data
   !REAL(KIND=RP) :: maxres
   INTEGER(KIND=IB)            :: kmin         ! Min number of layers
   INTEGER(KIND=IB)            :: kmax         ! Max number of layers
   REAL(KIND=RP)               :: hmin               ! Min allowed layer thickness in m
   REAL(KIND=RP)               :: hmx                ! Max depth in m
   REAL(KIND=RP),PARAMETER     :: fact     = 1.00_RP ! factor for rotated space perturbation
   REAL(KIND=RP),PARAMETER     :: factdelay= 1.50_RP ! shrinking factor for delayed rejection (>1.)
   REAL(KIND=RP),DIMENSION(NPL):: maxpert  = 0._RP
   REAL(KIND=RP),DIMENSION(NPL):: pertsd   = 0._RP 
!   REAL(KIND=RP),DIMENSION(NPL):: pertsdsc = (/ 100._RP,20._RP /)
   REAL(KIND=RP),DIMENSION(NPL):: pertsdsc = 30._RP
   REAL(KIND=RP),DIMENSION(NPL):: w_sl     = 0._RP   ! initial slice width 

   REAL(KIND=RP)               :: area_bn
   REAL(KIND=RP)               :: area_cr

!!
!!  Autoregressive model prior variables:
!!
   INTEGER(KIND=IB),PARAMETER :: NARFP      = 1                ! Order of each process
   REAL(KIND=RP)              :: armx       = 0.2_RP           ! Max AR and ARI model range (R units)
   REAL(KIND=RP), DIMENSION(NARFP*NRF):: minlimar  = 0._RP
   REAL(KIND=RP), DIMENSION(NARFP*NRF):: maxlimar  = 0._RP
   REAL(KIND=RP), DIMENSION(NARFP*NRF):: maxpertar = 0._RP
   REAL(KIND=RP), DIMENSION(NARFP*NRF):: pertarsd  = 0._RP
   REAL(KIND=RP), DIMENSION(NARFP*NRF):: pertarsdsc= 18._RP

!!
!!  Standard deviation prior variables:
!!
   REAL(KIND=RP), DIMENSION(NRF):: minlimsd  = 0._RP
   REAL(KIND=RP), DIMENSION(NRF):: maxlimsd  = 0._RP
   REAL(KIND=RP), DIMENSION(NRF):: maxpertsd = 0._RP
   REAL(KIND=RP), DIMENSION(NRF):: pertsdsd  = 0._RP
   REAL(KIND=RP), DIMENSION(NRF):: pertsdsdsc= 20._RP
!!
!!  C factor prior variables:
!!
   REAL(KIND=RP), DIMENSION(NRF):: minlimcf  = 0._RP
   REAL(KIND=RP), DIMENSION(NRF):: maxlimcf  = 0._RP
   REAL(KIND=RP), DIMENSION(NRF):: maxpertcf = 0._RP
   REAL(KIND=RP), DIMENSION(NRF):: pertcfsd  = 0._RP
   REAL(KIND=RP), DIMENSION(NRF):: pertcfsdsc= 20._RP
!!
!!  dt prior variables:
!!
   REAL(KIND=RP), DIMENSION(NRF,NTIME):: timebuf  = 0._RP
   REAL(KIND=RP)                      :: minlimdt  = 0._RP
   REAL(KIND=RP)                      :: maxlimdt  = 0._RP
   REAL(KIND=RP)                      :: maxpertdt = 0._RP
   REAL(KIND=RP)                      :: pertdtsd  = 0._RP
   REAL(KIND=RP)                      :: pertdtsdsc= 5._RP
!!
!!  Water resistivity and depth prior variables:
!!
   REAL(KIND=RP) :: minlimwres  = 0._RP
   REAL(KIND=RP) :: maxlimwres  = 0._RP
   REAL(KIND=RP) :: maxpertwres = 0._RP
   REAL(KIND=RP) :: pertwressd  = 0._RP
   REAL(KIND=RP) :: pertwressdsc= 20._RP
   REAL(KIND=RP) :: minlimwdepth  = 0._RP
   REAL(KIND=RP) :: maxlimwdepth  = 0._RP
   REAL(KIND=RP) :: maxpertwdepth = 0._RP
   REAL(KIND=RP) :: pertwdepthsd  = 0._RP
   REAL(KIND=RP) :: pertwdepthsdsc= 5._RP   
  
   CHARACTER(len=64)  :: filebase
   CHARACTER(len=64)  :: mapfilein
   INTEGER(KIND=IB)   :: filebaselen
   INTEGER(KIND=IB)   :: mapfilelen
   CHARACTER(LEN=200) :: infile1
   CHARACTER(LEN=64)  :: logfile
   CHARACTER(LEN=64)  :: seedfile
   CHARACTER(LEN=64)  :: mapfile
   CHARACTER(len=64)  :: icovfile
   CHARACTER(LEN=64)  :: obsfile
   CHARACTER(LEN=64)  :: parfile
   CHARACTER(LEN=64)  :: repfile
   CHARACTER(LEN=64)  :: dispfile
   CHARACTER(LEN=64)  :: sdfile
   CHARACTER(LEN=64)  :: samplefile
   CHARACTER(LEN=64)  :: lincovfile
   CHARACTER(LEN=64)  :: modname = 'sample.geom'

!!
!! Parallel Tempering parameters
!!
   INTEGER(KIND=IB),PARAMETER                :: NTPT       =  9_IB    ! # tempering levels (temperatures)
   INTEGER(KIND=IB),PARAMETER                :: NPTCHAINS  =  10_IB    ! # parallel tempering chains
   INTEGER(KIND=IB),DIMENSION(NTPT),PARAMETER:: NCHAINT    = (/2_IB,1_IB,1_IB,1_IB,1_IB,1_IB,1_IB,1_IB,1_IB/) 
!   INTEGER(KIND=IB),PARAMETER               :: NTPT       =  1_IB    ! # tempering levels (temperatures)
!   INTEGER(KIND=IB),PARAMETER               :: NPTCHAINS  =  1_IB   ! # parallel tempering chains
!   INTEGER(KIND=IB),DIMENSION(NTPT),PARAMETER:: NCHAINT   = (/1_IB/) 
   INTEGER(KIND=IB)                        :: ncswap     = 0_IB    ! Temp swap accepted counter
   INTEGER(KIND=IB)                        :: ncswapprop = 0_IB    ! Temp swap proposed counter
   INTEGER(KIND=IB)                        :: nccross    = 0_IB    ! Crossover accepted counter
   INTEGER(KIND=IB)                        :: nccrossprop= 0_IB    ! Crossover accepted counter
   REAL(KIND=RP),PARAMETER                 :: dTlog      = 1.25_RP  ! Temperature increment
   REAL(KIND=RP),DIMENSION(NPTCHAINS)      :: beta_pt              ! Temperature array parallel tempering

!!
!!  Sampling specific parameters
!!
   INTEGER(KIND=IB)           :: NFPMX      = (NLMX * NPL) + (NPL-1)
   INTEGER(KIND=IB)           :: NVV        = ((NLMX*NPL)+NPL-1)*((NLMX*NPL)+NPL-1)*NLMX
   INTEGER(KIND=IB)           :: NVV2       = ((NLMX*NPL)+NPL-1)
   INTEGER(KIND=IB)           :: NSDEVM     = ((NLMX*NPL)+NPL-1)*NLMX
   INTEGER(KIND=IB)           :: NMISC      = 2+NRF ! 2 extra parameters water resistivity and depth and NRF params for dt
   INTEGER(KIND=IB)           :: ioutside   = 0
   INTEGER(KIND=IB)           :: ireject    = 0, iaccept = 0, iaccept_delay = 0, ireject_delay = 0
   INTEGER(KIND=IB)           :: i_bd         ! Birth-Death track (0=MCMC, 1=birth, 2=death)
   INTEGER(KIND=IB)           :: i_sdpert = 0 ! if sigma is perturbed, don't compute forward model

!!
!!  Convergence parameters
!!
   INTEGER(KIND=IB)       :: iconv    = 0       ! Convergence switch slaves
   INTEGER(KIND=IB)       :: iconv2   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iconv3   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iarfail  = 0       ! Tracks failure of AR model when predicted AR series too large

!!
!! RJMCMC parameters
!!
   INTEGER(KIND=IB),PARAMETER    :: NCHAIN     = 1E9_IB  ! # iterations (max # MCMC steps)
   INTEGER(KIND=IB),PARAMETER    :: ICHAINTHIN = 5E0_IB  ! Chain thinning interval !Romina!
   INTEGER(KIND=IB),PARAMETER    :: NKEEP1     = 1E1_IB  ! Number models to keep before writing
   INTEGER(KIND=IB)              :: NKEEP
   INTEGER(KIND=IB),PARAMETER    :: NAP        = 10      ! Misc parameters in sample (for bookeeping)
   INTEGER(KIND=IB),PARAMETER    :: NDM        = 100     ! No. steps in lin rot est

!!
!! Annealing burn-in parameters (sets beta schedule)
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: beta                 ! Inverse Temprerature
   INTEGER(KIND=IB),PARAMETER            :: NTEMP    = 2e3       ! Number of beta values
   INTEGER(KIND=IB),PARAMETER            :: NBAL     = 3e1       ! Number of balanching steps at each T
   REAL(KIND=RP)                         :: beta1gl  = 0.01_RP   ! Global inverse T to define annealing schedule (start)
   REAL(KIND=RP)                         :: beta4gl  = 1.00_RP   ! Global inverse T to define annealing schedule (end for 3 exp legs)

!!
!!  Rotation matrices 
!!
   REAL(KIND=RP),DIMENSION((NLMX*NPL)+NPL-1,(NLMX*NPL)+NPL-1,NLMX):: Cov0   ! Linear parameter cov mat
   REAL(KIND=RP),DIMENSION((NLMX*NPL)+NPL-1,(NLMX*NPL)+NPL-1,NLMX):: VV     ! Linear rotation
   REAL(KIND=RP),DIMENSION((NLMX*NPL)+NPL-1,NLMX)                 :: sdevm  ! Std dev for perturbations
   REAL(KIND=RP),DIMENSION((NLMX*NPL)+NPL-1,NLMX)                 :: dmbest ! Best step-size for derivative est
   INTEGER(KIND=IB),DIMENSION(NLMX)                               :: idmbest! is best step-size already computed in this dimension?
!!
!!  Structures for objects and data 
!!
   TYPE :: objstruc
      SEQUENCE
      REAL(KIND=RP),DIMENSION((NLMX*NPL)+NPL-1):: par     ! Forward parameters, contains depths (z) in meter
      REAL(KIND=RP),DIMENSION(NLMX)            :: z          ! Forward parameters
      REAL(KIND=RP),DIMENSION(NLMX)            :: h          ! Forward parameters
	  REAL(KIND=RP)                          :: wres       ! Water resistivity
  	  REAL(KIND=RP)                          :: wdepth     ! Water depth
      REAL(KIND=RP),DIMENSION(NARFP*NRF)    :: arpar      ! AR model forward parameters
      INTEGER(KIND=IB),DIMENSION(NARFP*NRF) :: idxar      ! AR model forward parameters
      REAL(KIND=RP),DIMENSION(NRF)          :: sdpar      ! Std dev model forward parameters
      REAL(KIND=RP),DIMENSION(NRF)          :: cf         ! Scaling factor
      REAL(KIND=RP),DIMENSION(NRF)          :: dt         ! time shift
      REAL(KIND=RP),DIMENSION(NPL-1)        :: g          ! Acoustic parameters birth-death layer
      REAL(KIND=RP),DIMENSION(NPL-1)        :: gp         ! Acoustic parameters birth-death layer perturbed
      REAL(KIND=RP),DIMENSION(NPL-1,NPL-1)  :: Chat,Chati ! Covariance matrix for perturbing one BD layer
      REAL(KIND=RP)                         :: detChat
      REAL(KIND=RP)                         :: beta
      REAL(KIND=RP)                         :: cpt        ! # chains at this T
      REAL(KIND=RP)                         :: pariwhichrot ! for delayed rejection: rotatet (or not rotated) parameter currently perturbed
      INTEGER(KIND=IB)                      :: k          ! Layer dimension
      INTEGER(KIND=IB)                      :: NFP        ! Number forward parameters
      REAL(KIND=RP)                         :: logL       ! log likelihood
      REAL(KIND=RP)                         :: logPr      ! log Prior probability ratio
      REAL(KIND=RP)                         :: lognorm    ! Data covariance matrices
      INTEGER(KIND=IB)                      :: ireject_bd = 0
      INTEGER(KIND=IB)                      :: iaccept_bd = 0
   END TYPE objstruc

!!
!! Data structure
!!
   TYPE :: datastruc
      SEQUENCE
      REAL(KIND=RP),DIMENSION(NRF,NTIME)     :: Dobs   = 0._RP ! Observed data for one logL eval.
      REAL(KIND=RP),DIMENSION(NRF,NTIME)     :: sdev   = 0._RP ! Observed data for one logL eval.
      REAL(KIND=RP),DIMENSION(NTIME)         :: time   = 0._RP   ! Observed/samples times
      REAL(KIND=RP),DIMENSION(NRF,NTIME)     :: Drep   = 0._RP ! Replica data for trial model
      REAL(KIND=RP),DIMENSION(NRF,NTIME)     :: res    = 0._RP ! Data residuals for trial model
      REAL(KIND=RP),DIMENSION(NRF,NTIME)     :: dar    = 0._RP ! Autoregressive model replica data
      REAL(KIND=RP),DIMENSION(NRF)           :: lognorm=0._RP  ! Data lognorm
      REAL(KIND=RP),DIMENSION(NRF)           :: logdet =0._RP  ! Data log-determinant
      INTEGER(KIND=IB)                       :: NTIME  = NTIME  ! # data/angles (struc copy)
      INTEGER(KIND=IB),DIMENSION(NRF)        :: NDAT   = 0     ! No. data at each receiver
   END TYPE datastruc

      TYPE :: covstruc
!    SEQUENCE
       REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: Cdi   ! Inverse covariance matrix
       REAL(KIND=RP) :: logdet=0._RP   ! determinant of covariance matrix
      END TYPE covstruc
   
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE:: icount

!!
!!  Global variables
!!
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: sample                    ! Posterior sample
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE   :: tmpmap                    ! temporary for reading map
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE   :: buf_save_snd,buf_save_rcv ! Buffers for MPI sending
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE:: buffer1                   !
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE   :: buffer2,buffer3           !

!!
!!  MPI global variables
!!
   INTEGER(KIND=IB)            :: rank,NTHREAD,ncount,ncount2,ierr
   INTEGER(KIND=IB), PARAMETER :: src = 0_IB
   INTEGER                     :: to,from,tag,COMM
   INTEGER                     :: status(MPI_STATUS_SIZE)
   INTEGER(KIND=IB)            :: isize1,isize2,isize3

   INTERFACE
      FUNCTION RANDPERM(num)
         USE data_type, ONLY : IB
         IMPLICIT NONE
         INTEGER(KIND=IB), INTENT(IN) :: num
         INTEGER(KIND=IB), DIMENSION(num) :: RANDPERM
      END FUNCTION RANDPERM
   END INTERFACE

END MODULE RJMCMC_COM
  
!=======================================================================

PROGRAM  RJMCMC_PLANE

!=======================================================================
USE MPI
USE RJMCMC_COM
USE NR
USE mar_data
USE calc_mod
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: i,j,ipar,ilay,imcmc,ithin,isource,ikeep,ik,irf
INTEGER(KIND=IB)  :: idat,it,it2,ic

TYPE (objstruc),DIMENSION(NPTCHAINS) :: obj      ! Objects in likelihood box
TYPE (covstruc),DIMENSION(NRF)        :: cov      ! Structure for inverse data Covariance Matrices

TYPE (objstruc)                      :: objmax   ! Objects in likelihood box
TYPE (datastruc)                     :: dat      ! Data
REAL(KIND=RP)                        :: ran_uni,ran_nor,logLG

INTEGER(KIND=IB),DIMENSION(ICHAINTHIN):: idxchain !! Chain thinning array (use random 
INTEGER(KIND=IB):: NTHIN

!!---------------------------------------------------------------------!
!!     MPI stuff:
!!---------------------------------------------------------------------!
!!
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds

REAL(KIND=RP)               :: tstart, tend              ! Overall time 
REAL(KIND=RP)               :: tstart2, tend2            ! Time for one forward model computation
REAL(KIND=RP)               :: tstartsnd, tendsnd        ! Communication time
REAL(KIND=RP)               :: tstartcmp, tendcmp, tcmp  ! Forward computation time

      CALL MPI_INIT( ierr )
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NTHREAD, ierr )

      OPEN(UNIT=20,FILE=filebasefile,FORM='formatted',STATUS='OLD',ACTION='READ')
       READ(20,*) filebaselen
       READ(20,*) filebase
       READ(20,*) IMAP
       READ(20,*) SAMPLEPPD
       READ(20,*) ICOV
       READ(20,*) IAR
       READ(20,*) ICF
       READ(20,*) IWRES
       READ(20,*) IWDEPTH
       READ(20,*) IDT
       READ(20,*) hmin
       READ(20,*) hmx
       READ(20,*) minres
       READ(20,*) maxres
       READ(20,*) kmin
       READ(20,*) kmax
       READ(20,*) DIPOLE_PARTS
      CLOSE(20)

infile1        = filebase(1:filebaselen) // '.mar'
logfile        = filebase(1:filebaselen) // '_RJMH.log'
seedfile       = filebase(1:filebaselen) // '_seeds.log'
mapfile        = filebase(1:filebaselen) // '_map.dat'
repfile        = filebase(1:filebaselen) // '_rep.dat'
obsfile        = filebase(1:filebaselen) // '_obs.dat'
dispfile       = filebase(1:filebaselen) // '_disp.dat'
icovfile       = 'Cdinv_rec'  !! Inverse data cov mat
sdfile         = filebase(1:filebaselen) // '_sigma.txt'
samplefile     = filebase(1:filebaselen) // '_sample.txt'
parfile        = filebase(1:filebaselen) // '_par.txt'



      IF(rank == src)WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
      IF(rank == src)WRITE(6,*) '~~~                                                        ~~~  '
      IF(rank == src)WRITE(6,*) '~~~             Reversible Jump MCMC Sampling              ~~~  '
      IF(rank == src)WRITE(6,*) '~~~                                                        ~~~  '
      IF(rank == src)WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
      IF(rank == src)WRITE(6,*) '...running on ',NTHREAD,' cores'

NKEEP = NKEEP1*NCHAINT(1)    !! Keep all chains at T=1

ncount = NKEEP*(NFPMX+NAP+2*NRF+(NARFP*NRF)+NMISC) ! max. no. of parameters+ bookkeeping params + no. of cf and sd parameters + no of AR parameters + misc extra params like water depth 
ncount2= NFPMX+NAP+2*NRF+(NARFP*NRF)+NMISC
ALLOCATE( sample(NKEEP,ncount2) )
ALLOCATE( tmpmap(NFPMX+2*NRF+(NARFP*NRF)+1+NMISC) )
ALLOCATE( buffer1(1+MPI_BSEND_OVERHEAD),buffer2(ncount+MPI_BSEND_OVERHEAD),buffer3(1+MPI_BSEND_OVERHEAD) )
ALLOCATE( buf_save_snd(ncount),buf_save_rcv(ncount) )
buf_save_snd = 0._RP
buf_save_rcv = 0._RP
sample       = 0._RP
!------------------------------------------------------------------------
!  Set up tempering schedule
!------------------------------------------------------------------------
it2 = 1_IB
DO it=1,NTPT
DO ic=1,NCHAINT(it)
  beta_pt(it2) = 1._RP/dTlog**REAL(it-1_IB,RP)
  obj(it2)%beta = beta_pt(it2)
  obj(it2)%cpt  = REAL(NCHAINT(it),RP)
  it2 = it2 + 1
ENDDO
ENDDO
IF(rank == src)WRITE(6,*) 'beta_pt = ',beta_pt

!------------------------------------------------------------------------
!  Read in data
!------------------------------------------------------------------------
ETA=EPSILON(1.)
TOL=TINY(1.)/ETA
NLAY=0
NDAT=0
TNT=0
TX_DESIG=''
MAXNT=0
DO ic=1,100
  UN_NO(ic)=.FALSE.
ENDDO
       CALL LOAD(2,infile1)
timebuf = time
dat%Dobs = volt
dat%sdev = abs(cvolt*error/100.)
dat%NDAT = NT
obj%wres = w_res
obj%wdepth = w_depth

!!
IF(rank == src) write(6,*) 'Set Prior bounds ...'
!!
!!           z      rho 
minlim = (/ hmin,   log10(minres) /)
maxlim = (/ hmx,  log10(maxres) /)

!kmin = 1
!kmax = NLMX

!! Maximum perturbation size
maxpert = maxlim-minlim
pertsd = maxpert/pertsdsc
w_sl = maxpert/10._RP

IF(IAR == 1)THEN
   !! Set prior and proposal scaling for AR model:
   minlimar   = -0.6000_RP
   maxlimar   =  0.9999_RP
   pertarsdsc =  10._RP
   maxpertar  = maxlimar-minlimar
   pertarsd   = maxpertar/pertarsdsc
ENDIF

IF(ICOV == 1 .OR. ICOV == 3)THEN
   !! Set prior and proposal scaling for data error standard deviations:
   minlimsd   = LOG(1.e-11_RP)
   maxlimsd   = LOG(1.e-5_RP)  
   maxpertsd  = maxlimsd-minlimsd
   pertsdsd   = maxpertsd/pertsdsdsc
elseif (ICOV == 5) then
    minlimsd   = 1.0_RP
    maxlimsd   = 10.0_RP 
    maxpertsd  = maxlimsd-minlimsd
    pertsdsd   = maxpertsd/pertsdsdsc  
ENDIF

IF(ICOV == 0)THEN
  IF(rank == src)WRITE(6,*) '  ~~~  Implicit std dev sampling                         ~~~  '
ELSEIF(ICOV == 1)THEN
  IF(rank == src)WRITE(6,*) '  ~~~  Implicit std dev sampling non-stat std dev        ~~~  '
      ELSEIF(ICOV == 2)THEN
       IF(rank == src)WRITE(6,*) '  ~~~  Empirical residual cov mat with factors  ~~~  '
      CALL READCOVI(cov,obj,dat) 
      ELSEIF(ICOV == 3 .OR. ICOV ==5)THEN
       IF(rank == src)WRITE(6,*) '  ~~~  Empirical residual cov or main diagonal from stacking mat with magnitude scaling ~~~  '
       IF(rank == src)WRITE(6,*) '  ~~ minlimsd   = ',minlimsd
       IF(rank == src)WRITE(6,*) '  ~~ maxlimsd   = ',maxlimsd
       CALL READCOVI(cov,obj,dat)
      elseif(icov ==4)then
       IF(rank == src)WRITE(6,*) '  ~~~  Stacking standard deviation in cov mat ~~~  '
       CALL READCOVI(cov,obj,dat) 
      ENDIF

IF(IDT >= 1)THEN
   !! Set prior and proposal scaling for data error standard deviations:
   minlimdt   = -0.001_RP
   maxlimdt   = 0.001_RP
   maxpertdt  = maxlimdt-minlimdt
   pertdtsd   = maxpertdt/pertdtsdsc
ENDIF
IF(ICF >= 1)THEN
   !! Set prior and proposal scaling for data error standard deviations:
   minlimcf   = 0.9_RP
   maxlimcf   = 1.15_RP
   maxpertcf  = maxlimcf-minlimcf
   pertcfsd   = maxpertcf/pertcfsdsc
ENDIF
IF(IWRES == 1)THEN ! Romina
   !! Set prior and proposal scaling for water resistivity:
   ! minlimwres   = 0.39_RP
   ! maxlimwres   = 0.42_RP
   minlimwres   = 0.24_RP
   maxlimwres   = 0.28_RP
   maxpertwres  = maxlimwres-minlimwres
   pertwressd   = maxpertwres/pertwressdsc
ENDIF
IF(IWDEPTH == 1)THEN
   !! Set prior and proposal scaling for water depth:
   minlimwdepth   = 30.0_RP
   maxlimwdepth   = 50.0_RP
   maxpertwdepth  = maxlimwdepth-minlimwdepth
   pertwdepthsd   = maxpertwdepth/pertwdepthsdsc
ENDIF
!!------------------------------------------------------------------------
!!
!!  Print sampling parameters to screen for logging
!!
IF(rank == src)THEN
   WRITE(6,210) 'ICOV                 :   ',ICOV
   IF(ICOV == 1)THEN
      WRITE(6,204) '...sampling over standard deviations.                   '
   ELSEIF(ICOV == 2)THEN
      WRITE(6,204) '...sampling over standard deviations with joint logL    '
   ENDIF
   WRITE(6,210) 'NKEEP                :   ',NKEEP
   WRITE(6,210) 'ICHAINTHIN           :   ',ICHAINTHIN
   WRITE(6,210) 'ILINROT              :   ',ILINROT
   WRITE(6,210) 'IDELAY               :   ',IDELAY
   WRITE(6,210) 'ICF                  :   ',ICF
   WRITE(6,210) 'IEXCHANGE            :   ',IEXCHANGE
   WRITE(6,210) 'NTPT                 :   ',NTPT
   WRITE(6,210) 'NPTCHAINS            :   ',NPTCHAINS
   WRITE(6,210) 'NCHAINT              :   ',NCHAINT
   WRITE(6,210) 'IAR                  :   ',IAR
   WRITE(6,210) 'Order                :   ',NARFP
   WRITE(6,210) 'No. Receivers        :   ',NRF
   WRITE(6,210) 'Max no. samples      :   ',NTIME
   WRITE(6,210) 'No. samples per rec. :   ',dat%NDAT
   CALL FLUSH(6)
   WRITE(6,209) 'Sample file:             ',samplefile
   WRITE(6,*) ''
   WRITE(6,*) 'NFPMX         = ',NFPMX
   WRITE(6,*) 'kmin and kmax = ',kmin,kmax
   WRITE(6,*) 'dipole parts  = ',DIPOLE_PARTS
   WRITE(6,*) 'min thickness and resistivity:  '
   WRITE(6,201) minlim
   WRITE(6,*) 'max thickness and resistivity:  '
   WRITE(6,201) maxlim
   WRITE(6,*) 'min/maxlim wres:  '
   WRITE(6,201) minlimwres,maxlimwres
   WRITE(6,*) 'min/maxlim wdepth:  '
   WRITE(6,201) minlimwdepth,maxlimwdepth
   WRITE(6,*) 'min/maxlim sigma:  '
   WRITE(6,201) minlimsd,maxlimsd
   WRITE(6,*) ''
   WRITE(6,*) 'Done reading data.'
   WRITE(6,*) ''
   CALL FLUSH(6)
   WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
   WRITE(6,*) ''
ENDIF

!!
!! Initialize random seeds on each core (Call RANDOM_SEED only once in the whole code. PARALLEL_SEED calls it)
!!
!CALL RANDOM_SEED
CALL PARALLEL_SEED()
!!
!! Make cooling schedule
!!
IF(ICOOL == 1) THEN
   IF(rank == src)THEN
      WRITE(6,*) ''
      WRITE(6,*) '  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
      WRITE(6,*) 'Cooling turned on:'
      WRITE(6,*) 'NTEMP = ',NTEMP
      WRITE(6,*) 'NBAL  = ',NBAL
      WRITE(6,*) 'beta1 = ',beta1gl
      WRITE(6,*) 'beta4 = ',beta4gl
      WRITE(6,*) ''
   ENDIF
   ALLOCATE( beta(NTEMP) )
   CALL MAKE_BETA(beta1gl,beta4gl,NTEMP)
ENDIF

!!
!! some factors for the acceptance ratio in birth/death case:
!!
DO ic = 1,NPTCHAINS
  obj(ic)%Chat    = 0._RP
  obj(ic)%Chati   = 0._RP
  obj(ic)%detChat = 1._RP
  DO ipar = 1,NPL-1
    obj(ic)%Chat(ipar,ipar)  = pertsd(ipar+1)**2._RP
    obj(ic)%Chati(ipar,ipar) = 1._RP/pertsd(ipar+1)**2._RP
    obj(ic)%detChat = obj(ic)%detChat * obj(ic)%Chat(ipar,ipar)
  ENDDO
  obj(ic)%logPr = SUM(LOG(maxlim(2:NPL)-minlim(2:NPL))) 
  obj(ic)%lognorm = LOG(SQRT((2._RP*PI2)**REAL(NPL-1,RP)*obj(ic)%detChat))
ENDDO

!! -------------------------------------------------------------
!! Read mapfile to start:
!!
      IF(SAMPLEPPD == 1)THEN
      IF(rank == src)WRITE(6,*) 'Computing forward models only'
      OPEN(UNIT=21,FILE=mapfile,FORM='formatted',STATUS='OLD',ACTION='READ')
      OPEN(UNIT=51,FILE=repfile,FORM='formatted',STATUS='REPLACE',ACTION='WRITE',RECL=1024)
      tstart2 = MPI_WTIME()
       DO ! Loop over file contents and go to line 520 when reaching end of file
           READ(21,*,end=520) tmpmap
           obj(1)%k   = INT(tmpmap(1),IB)
           obj(1)%par   = 0._RP
           obj(1)%NFP = (obj(1)%k * NPL) + (NPL-1)
           obj(1)%par = 0._RP
           obj(1)%h   = 0._RP
           obj(1)%par(1:obj(1)%NFP) = tmpmap(2:obj(1)%NFP+1)
           obj(1)%sdpar = 0._RP
           obj(1)%arpar = 0._RP
           obj(1)%cf = 1._RP
           obj(1)%dt = 0._RP
           obj(1)%par(1:obj(1)%NFP) = tmpmap(2:obj(1)%NFP+1)
          IF(ICOV == 1 .OR. ICOV == 2 .OR. ICOV ==6) then
           obj(1)%sdpar = LOG(tmpmap(NFPMX+1+1:NFPMX+1+NRF))
          elseif (ICOV==5) then
           obj(1)%sdpar = tmpmap(NFPMX+1+1:NFPMX+1+NRF)
          endif
          IF(ICF >= 1) obj(1)%cf = tmpmap(NFPMX+1+NRF+1:NFPMX+1+2*NRF)
          IF(IAR == 1)THEN
             obj(1)%arpar = tmpmap(NFPMX+1+2*NRF+1:NFPMX+2*NRF+(NARFP*NRF)+1)
             obj(1)%idxar = 1
               DO ipar = 1,NRF
                IF(obj(1)%arpar(ipar) < minlimar(ipar)) obj(1)%idxar(ipar) = 0
               ENDDO!
           ENDIF
            if(iwres==1) obj(1)%wres = tmpmap(NFPMX+1+2*NRF+1+(NARFP*NRF))
              if(iwdepth==1) obj(1)%wdepth=tmpmap(NFPMX+1+2*NRF+1+(NARFP*NRF)+1)
             if(idt>=1) obj(1)%dt = tmpmap(NFPMX+1+2*NRF+1+(NARFP*NRF)+2:NFPMX+3*NRF+1+(NARFP*NRF)+2)
      
              CALL CALCH(obj(1))
     
                CALL CHECKBOUNDS(obj(1))
                IF(ioutside == 1)WRITE(*,*)'FAILED on bounds'

             CALL LOGLHOOD(obj(1),dat,cov)
             IF(rank == src) then
                  WRITE(51,208) obj(1)%logL
                  DO i = 1,NRF
                   WRITE(51,208) dat%Drep(i,:)
                  ENDDO
             endif
208 FORMAT(5000ES20.10)

       ENDDO
        tend2 = MPI_WTIME()
      IF(rank == src)WRITE(6,*) 'time = ',tend2-tstart2 
  520 IF(rank == src)WRITE(6,*) 'Finished.'
       CLOSE(21)
       CLOSE(51)
        STOP
      ENDIF ! SAMPLEPPD
!---------------------------------------------------------------------
write(6,*) mapfile
      OPEN(UNIT=20,FILE=mapfile,FORM='formatted',STATUS='OLD',ACTION='READ')
       READ(20,*) tmpmap
      CLOSE(20)

DO ic = 1,NPTCHAINS
  obj(ic)%k   = INT(tmpmap(1),IB)
  obj(ic)%NFP = (obj(ic)%k * NPL) + (NPL-1)
  !write(6,*) 'reading tmpmap, k and nfp:',obj(ic)%k,obj(ic)%NFP
  obj(ic)%par = 0._RP
  obj(ic)%h   = 0._RP
  obj(ic)%sdpar = 0._RP
  obj(ic)%arpar = 0._RP
  obj(ic)%cf = 1._RP
  obj(ic)%dt = 0._RP
  obj(ic)%par(1:obj(ic)%NFP) = tmpmap(2:obj(ic)%NFP+1)
  !write (6,*) 'parameter:',obj(ic)%par(1:obj(ic)%NFP)
      IF(ICOV == 1 .OR. ICOV == 3 .or. ICOV == 0 .or. ICOV==6) then
         obj(ic)%sdpar = LOG(tmpmap(NFPMX+1+1:NFPMX+1+NRF))
      elseif (ICOV==5) then
          obj(ic)%sdpar = tmpmap(NFPMX+1+1:NFPMX+1+NRF)	
      endif
      IF(ICF >= 1) obj(ic)%cf = tmpmap(NFPMX+1+NRF+1:NFPMX+1+2*NRF)
      IF(IAR == 1)THEN 
          obj(ic)%arpar = tmpmap(NFPMX+1+2*NRF+1:NFPMX+2*NRF+(NARFP*NRF)+1)
          obj(ic)%idxar = 1
       DO ipar = 1,NRF
         IF(obj(ic)%arpar(ipar) < minlimar(ipar)) obj(ic)%idxar(ipar) = 0
       ENDDO!
      ENDIF
      if(iwres==1) obj(ic)%wres = tmpmap(NFPMX+1+2*NRF+1+(NARFP*NRF))
      if(iwdepth==1) obj(ic)%wdepth=tmpmap(NFPMX+1+2*NRF+1+(NARFP*NRF)+1)
      if(idt>=1) obj(ic)%dt = tmpmap(NFPMX+1+2*NRF+1+(NARFP*NRF)+2:NFPMX+3*NRF+1+(NARFP*NRF)+2)

	  
!  IF(ICF == 1) obj(ic)%cf = 1._RP

!  obj(ic)%k   = 4 ! Romina
!  obj(ic)%NFP = (obj(ic)%k * NPL) + (NPL-1)
!  obj(ic)%par(1:obj(ic)%NFP) = (/ 26.0_RP, LOG10(0.91_RP), 208.0_RP, &
!      log10(2.02_RP), 443.0_RP, log10(0.97_RP), 480.6_RP, log10(5.7_RP), &
!      log10(1.17_RP) /)
  
  CALL CALCH(obj(ic))

ENDDO

ALLOCATE( icount(NTHREAD) )
icount = 0
tstart = MPI_WTIME()

      IF(IMAP == 1)THEN
       IF(rank == src)WRITE(6,*) 'IMAP activated, exiting after computing replica for MAP.'
       CALL CHECKBOUNDS(obj(1))
       IF(ioutside == 1)WRITE(*,*)'FAILED on bounds'
   
       tstart2 = MPI_WTIME()
       DO ic = 1,NPTCHAINS
        CALL LOGLHOOD(obj(1),dat,cov)
       ENDDO
      WRITE(6,*) 'save replica'
       tend2 = MPI_WTIME()
      IF(rank == src)CALL SAVEREPLICA(obj(1),dat,repfile,obsfile)
      IF(rank == src)WRITE(6,*) 'time = ',tend2-tstart2
      STOP
ELSE
  IF(rank == src)WRITE(6,*) 'Starting model:'
  IF(rank == src)CALL PRINTPAR(obj(1))
      DO ic = 1,NPTCHAINS
         tstart2 = MPI_WTIME()
       CALL LOGLHOOD(obj(ic),dat,cov)
         tend2 = MPI_WTIME()
      ENDDO
	  
  IF(rank == src)CALL PRINTPAR(obj(1))
  IF(rank == src)WRITE(6,*) 'logL = ',obj(1)%logL
  IF(rank == src)WRITE(6,*) 'time = ',tend2-tstart2
!  DO ic = 1,NPTCHAINS
!    CALL COUPLE_CR(obj(ic))
!  ENDDO

  IF(ioutside == 1)WRITE(*,*)'FAILED in couple'
  DO ic = 1,NPTCHAINS
    CALL CHECKBOUNDS(obj(ic))
  
  ENDDO  

  IF(ioutside == 1)then
  WRITE(*,*)'FAILED on bounds'
  WRITE(6,*)'FAILED on bounds'
  WRITE(6,*) obj(ic)%par, obj(ic)%wres, obj(ic)%wdepth
  STOP
  ENDIF

! this has to be changed because by not assigning VV properly it, the 1._RP are integrated in the sample file! It is like they are floating in space and then randomly settling into the sample file. Therefore you always need to assign in fortran and check bounds in the makefile -fcheck_bounds!
  !VV    = 0._RP
  !DO ik = kmin,kmax
  !DO ipar = 1,NFPMX
  !   VV(ipar,ipar,ik+1) = 1._RP
  !ENDDO
  !ENDDO

  !sdevm = 0._RP
  !     DO ik = kmin,kmax
  !      DO i = 1,ik*NPL+(NPL-1)
  !        ilay = CEILING(REAL(i,RP)/REAL(NPL,RP))
  !       ipar = i-(ilay-1)*NPL
  !      IF(ilay > ik) ipar = ipar + 1
  !      sdevm(i,ik+1) = pertsd(ipar)
  !      ENDDO
  !    ENDDO
!  IF(ILINROT == 1)THEN
!    idmbest = 0
!    dmbest  = 0._RP
!    CALL LINROT(obj(1),dat)
!  ENDIF
      tend2 = MPI_WTIME()
        IF(rank == src)WRITE(6,*) 'time for linrot = ',tend2-tstart2
      ENDIF
      CALL FLUSH(6)

!icount(rank+1) = icount(rank+1) + 1

!------------------------------------------------------------------------
!
!          ************ RJMCMC Sampling ************
!
! -----------------------------------------------------------------------
IF(rank == src)THEN
   WRITE(6,*) 'Starting RJMCMC sampling...'
   !!
   !!  File to save posterior sample
   !!
   OPEN(UNIT=40,FILE=samplefile,FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=1024)
   OPEN(UNIT=44,FILE=logfile,FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=1024)
ENDIF

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!     MASTER PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
IF(rank==src)THEN

      DO imcmc = 1,NCHAIN

         tstartsnd = MPI_WTIME()
       CALL SAVESAMPLE(logLG,isource,tcmp)
         tendsnd = MPI_WTIME()
      IF(MOD(imcmc,30+1)==0)THEN
       WRITE(*,*) ''
       WRITE(6,203) '   imcmc,     logL,      logPr,     chain  T  ,      k,   iacc/irej, iacc_bd,  time(send), source,  time(comp)'
       WRITE(6,203) '--------------------------------------------------------------------------------------------------------------'
      ENDIF

      WRITE(6,202) imcmc,sample(1,(/1,2,3/)),INT(sample(1,4)),sample(1,5+NFPMX+2*NRF+(NARFP*NRF)),&
                INT(sample(1,6+NFPMX+2*NRF+(NARFP*NRF))),tendsnd-tstartsnd,isource,tcmp
      WRITE(44,202) imcmc,sample(1,(/1,2,3/)),INT(sample(1,4)),sample(1,5+NFPMX+2*NRF+(NARFP*NRF)),&
                 INT(sample(1,6+NFPMX+2*NRF+(NARFP*NRF))),tendsnd-tstartsnd,isource,tcmp
   CALL FLUSH(6)
   CALL FLUSH(44)

      ENDDO
      IF(rank == src)THEN
       CLOSE(40)
      ENDIF

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!    SLAVE PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
ELSE
IF(ICOOL == 1)THEN
  !!
  !! Cooling only for one chain, no parallel tempering!
  !!
  objmax = obj(1)
  DO imcmc = 1,NTEMP
    DO ithin = 1,NBAL
	
      CALL EXPLORE_MH(obj(1),dat,cov,beta(imcmc))
      IF(obj(1)%logL > objmax%logL) objmax=obj(1)
    ENDDO
    IF(MOD(imcmc,500) == 0)THEN
      WRITE(*,213) rank,imcmc,'  logL=',obj(1)%logL,'  k=',obj(1)%k,'  objmax_logL=',objmax%logL,'  objmax_k=',objmax%k
!      IF(ILINROT == 1)THEN
!        CALL LINROT(obj(1),dat)
!        WRITE(*,*) rank,'updating Linrot while annealing, k=',obj(1)%k
!      ENDIF
    ENDIF
  ENDDO
  !! Write best model as starting model for all tempering chains
  DO ic = 1,NPTCHAINS
    obj(ic) = objmax
  ENDDO
ENDIF
ikeep   = 1_IB
idxchain = 0
DO imcmc = 1,NCHAIN

   IF(ikeep == 1)tstartcmp = MPI_WTIME()
   NTHIN = ICHAINTHIN
   DO ithin = 1,NTHIN
     DO ic = 1,NPTCHAINS

!write(6,*) 'Before exploring'
!write(6,*) ic, obj(ic)%par, obj(ic)%wres, obj(ic)%wdepth 
       
	   CALL EXPLORE_MH(obj(ic),dat,cov,obj(ic)%beta)
     ENDDO
     ncswap = 0
     ncswapprop = 0
     nccross = 0
     nccrossprop = 0
     IF(IEXCHANGE == 1)CALL TEMPSWP_MH(obj,dat)
   ENDDO
        if (ICOV <=4 .OR. ICOV==6) then
      DO ic = 1,NCHAINT(1)
          sample(ikeep,:) =  (/ obj(ic)%logL, obj(ic)%logPr, 1._RP, REAL(obj(ic)%k,RP), & ! 4 parameters
                           obj(ic)%par,obj(ic)%wres,obj(ic)%wdepth,&
                           obj(ic)%dt,&
                           exp(obj(ic)%sdpar),obj(ic)%cf,obj(ic)%arpar, & 
                           REAL(iaccept,RP)/REAL(iaccept+ireject,RP),REAL(obj(ic)%iaccept_bd,RP),&

                           REAL(obj(ic)%ireject_bd,RP),REAL(i_bd,RP),REAL(ic,RP),REAL(rank,RP) /)
! write(6,*) 'keeping sample'
! WRITE(6,207) sample(ikeep,:)

         ikeep = ikeep + 1_IB
      ENDDO
         else
       DO ic = 1,NCHAINT(1)
         sample(ikeep,:) =  (/ obj(ic)%logL, obj(ic)%logPr, 1._RP, REAL(obj(ic)%k,RP), & ! 4 parameters
                           obj(ic)%par,obj(ic)%wres,obj(ic)%wdepth,&
                           obj(ic)%dt,&
                           obj(ic)%sdpar,obj(ic)%cf,obj(ic)%arpar, & 
                           REAL(iaccept,RP)/REAL(iaccept+ireject,RP),REAL(obj(ic)%iaccept_bd,RP),&
                           REAL(obj(ic)%ireject_bd,RP),REAL(i_bd,RP),REAL(ic,RP),REAL(rank,RP) /)

        ikeep = ikeep + 1_IB
       ENDDO
      endif   
   IF(ikeep == NKEEP+1)THEN
     tendcmp = MPI_WTIME()
     tcmp = tendcmp-tstartcmp
     CALL SAVESAMPLE(logLG,isource,tcmp)
     ikeep  = 1
     sample = 0._RP
     IF(rank == 1)WRITE(6,214)'iaccept_bd         = ',(obj(ic)%iaccept_bd,ic=1,NPTCHAINS)
     IF(IEXCHANGE == 1)THEN 
       IF(rank == 1)WRITE(6,215)'T swap accept      = ',REAL(ncswap,RP)/REAL(ncswapprop,RP)
     ENDIF
     IF(IDELAY == 1) THEN
       IF(rank == 1)WRITE(6,217)'Delayed rej. ratio = ',REAL(iaccept_delay,RP)/REAL(iaccept_delay+ireject_delay,RP)
     ENDIF
     CALL FLUSH(6)
   ENDIF
ENDDO

ENDIF !! MPI ENDIF

201 FORMAT(200F12.4)
202 FORMAT(I8,3F16.6,I8,F16.6,I8,1F13.3,I8,1F13.3)
203 FORMAT(A129)
204 FORMAT(A56)
205 FORMAT(10F10.4)
207 FORMAT(500ES18.8)
209 FORMAT(A26,A40)
210 FORMAT(A26,20I4)
211 FORMAT(20000ES12.4)
212 FORMAT(I4,A8,F12.4,A5,I3)
213 FORMAT(I4,I6,A8,F12.4,A5,I3,A14,F12.4,A11,I3)
214 FORMAT(a21,18I6)
215 FORMAT(a21,F8.4)
216 FORMAT(a21,F8.4)
217 FORMAT(a21,F8.4)
CALL MPI_FINALIZE( ierr ) 

END PROGRAM RJMCMC_PLANE

!=======================================================================
SUBROUTINE LOGLHOOD(obj,dat,cov)
!=======================================================================
USE RJMCMC_COM
USE MPI
USE mar_data
IMPLICIT NONE

INTEGER(KIND=IB)                 :: irf,ipar,id,ilay,iparcur,isamp
TYPE (objstruc)                  :: obj
TYPE (covstruc),DIMENSION(NRF) :: cov
TYPE (datastruc)                 :: dat
REAL(KIND=RP),DIMENSION(NRF)     :: Ltmp
INTEGER(KIND=IB),DIMENSION(obj%k):: idxh
INTEGER(KIND=IB),DIMENSION(NRF)  :: idxrand
REAL(KIND=RP)                    :: tstart, tend, tcmp   ! Overall time 

      CALL CALCH(obj)
!tstart = MPI_WTIME()

       NLAY = obj%k+1

!       write(6,*) 'tness', tness
      IF (ALLOCATED(TNESS)) DEALLOCATE(TNESS,RES,ANIS)
      IF (Allocated(res)) deallocate(res,anis)
      if (obj%k == 0) then
       ALLOCATE(RES(NLAY),ANIS(NLAY))
      else
       ALLOCATE(TNESS(NLAY-1),RES(NLAY),ANIS(NLAY))
       tness = 0.
      endif
      res = 0.
      anis = 1.
      DO ipar=1,obj%k
        tness(ipar) = obj%h(ipar)
        res(ipar)   = 10.**obj%par(ipar*2)
      ENDDO
        res(obj%k+1) = 10.**obj%par(2*obj%k+1)
       
      if(iwres==1) w_res=obj%wres
      if(iwdepth==1) w_depth=obj%wdepth
      if(idt>=1) then
        do irf=1,NRF
           do isamp=1,dat%NDAT(irf)
             time(irf,isamp)=time(irf,isamp)+obj%dt(irf)
           enddo
        enddo
      endif
	  
!	 write(6,*) '! test martin input:'
 !    write(6,*) 'k,tness,res,wres,wdepth', obj%k,tness,res,w_res,w_depth
!	 write(6,*) obj%par, obj%wres, obj%wdepth
!if (w_res==1._RP) then
!	 write(6,*) '! test martin input:'
!     write(6,*) 'k,tness,res,wres,wdepth', obj%k,tness,res,w_res,w_depth
	 
!endif
!! ---------------------------------------------------------
!!  Compute predicted CSEM: CALL CALC (MARTIN main routine)
!!
      CALL CALC
!! ---------------------------------------------------------  
!	 write(6,*) obj%par, obj%wres, obj%wdepth
!	 write(6,*) obj%cf, obj%sdpar
      if(idt>=1) time=timebuf
	  
       dat%Drep = cvolt
      IF(ICF >= 1)THEN
       DO ipar = 1,NRF
  !        write(6,*) obj%cf(ipar)
          dat%Drep(ipar,:) = obj%cf(ipar)*dat%Drep(ipar,:)
       ENDDO
      ENDIF
       dat%res = (dat%Dobs-dat%Drep)

       iarfail = 0
      IF(IAR == 1)THEN
        dat%dar  = 0._RP
   !!
   !!  Compute autoregressive model
   !!
        DO irf = 1,NRF
         CALL ARFB(obj,dat,irf,1,dat%NDAT(irf))
        ENDDO
   
   !! Recompute replica data as ith autoregressive model
           dat%res = dat%res-dat%dar
   !! Check if predicted AR model data are outside max allowed bounds
        CALL CHECKBOUNDS_ARMX(dat)

      ENDIF
!!
!!  Compute log likelihood
!!
      IF(iarfail == 0)THEN
       Ltmp = 0._RP
         IF(ICOV == 0)THEN
          !!
          !! Implicitly sample over sigma
          !!
	!	write(6,*) 'Implicitly sample over sigma:'
     !      DO irf = 1,NRF
     !        Ltmp(irf) = -REAL(dat%NDAT(irf),RP)/2._RP*LOG(SUM(dat%res(irf,:)*dat%res(irf,:)))
	 !  write(6,*) 'L',irf,Ltmp(irf)
     !      ENDDO
		  !!
		  !! fixed sigma
		  !!
	!	  write(6,*) 'fixed sigma:'
           DO irf = 1,NRF
       Ltmp(irf) = -REAL(dat%NDAT(irf),RP)/2._RP*log(2._RP*PI2*exp(obj%sdpar(irf)))-&
       SUM((dat%res(irf,:)*dat%res(irf,:))/(2._RP*(exp(obj%sdpar(irf))**2._RP)))

           ENDDO
	
         ELSEIF(ICOV == 1 .OR. ICOV == 6)THEN
           !!
           !! Sample over sigma
           !! SIGMA IS IN LOG SPACE!
           !!
           DO irf = 1,NRF
           ! Gaussian
           Ltmp(irf) = LOG(1._RP/(2._RP*PI2)**(REAL(dat%NDAT(irf),RP)/2._RP)) &
                     -(SUM(dat%res(irf,1:dat%NDAT(irf))**2._RP)/ &
                      (2._RP*exp(obj%sdpar(irf))**2._RP)&
                      +REAL(dat%NDAT(irf),RP)*obj%sdpar(irf))
      
    ! Laplace ---------------------------------------------
    !   Ltmp(irf) = log(1._RP/2._RP**REAL(dat%NDAT(irf),RP))&
    !               -(SUM(abs(dat%res(irf,1:dat%NDAT(irf))))/ &
	!		       exp(obj%sdpar(irf))&
    !               +REAL(dat%NDAT(irf),RP)*obj%sdpar(irf))
	! ------------------------------------------------------
           ENDDO
        ELSEIF(ICOV == 2 .OR. ICOV == 4)THEN
              !! Empirical cov mat estimate
          DO irf = 1,NRF
           Ltmp(irf) = -DOT_PRODUCT(MATMUL(dat%res(irf,1:dat%NDAT(irf)),&
                       cov(irf)%Cdi),&
                       dat%res(irf,1:dat%NDAT(irf)))/2._RP
          ENDDO
		    ELSEIF(ICOV == 5)THEN
    !! Empirical cov mat estimate and magnitude scaling (scaling is done via std
    !! dev array...
       DO irf = 1,NRF
            Ltmp(irf) = -DOT_PRODUCT(MATMUL(dat%res(irf,1:dat%NDAT(irf)),&
            cov(irf)%Cdi),dat%res(irf,1:dat%NDAT(irf)))&
           /(2._RP*obj%sdpar(irf)**2._RP) - &
            REAL(dat%NDAT(irf),RP)*log(obj%sdpar(irf))
			! cov(irf)%logdet/2._RP
	!		write(6,*) 'irf,Ltmp(irf),log(obj%sdpar(irf))',irf,Ltmp(irf),log(obj%sdpar(irf))
       ENDDO
          ELSEIF(ICOV == 3)THEN
    !! Empirical cov mat estimate and magnitude scaling (scaling is done via std
    !! dev array...
       DO irf = 1,NRF
         Ltmp(irf)= -DOT_PRODUCT(MATMUL(dat%res(irf,1:dat%NDAT(irf)),&
                    cov(irf)%Cdi),dat%res(irf,1:dat%NDAT(irf)))&
                    /(2._RP*exp(obj%sdpar(irf))) - &
                    REAL(dat%NDAT(irf),RP)/2._RP*obj%sdpar(irf)
       ENDDO		  
        ENDIF
      obj%logL = SUM(Ltmp)
!	 write(6,*) 'Likelihood'
!     write(6,*) obj%logL
!	 write(6,*) 'Parameters'
!	 write(6,*) obj%par, obj%wres, obj%wdepth

      ELSE
       obj%logL = -HUGE(1._RP)
       iarfail = 0
      ENDIF

	  !	  if (w_res==1._RP) then
 !     write(6,*) 'wres and logl', obj%wres, obj%logL
!	  endif

RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD
!==============================================================================

SUBROUTINE TEMPSWP_MH(obj,dat)
!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                     :: ic,ic1,ic2,ik
INTEGER(KIND=IB),DIMENSION(NPTCHAINS):: idxrand
TYPE(objstruc),DIMENSION(NPTCHAINS)  :: obj
TYPE(objstruc)                       :: objtmp1,objtmp2
TYPE(datastruc)                      :: dat
REAL(KIND=RP)                        :: ran_uni1
REAL(KIND=RP)                        :: vol1,vol2   ! Volumes for each model
REAL(KIND=RP)                        :: logP1,logP2,logfactk !
REAL(KIND=RP)                        :: logratio    !
REAL(KIND=RP)                        :: betaratio   !

DO ic = 1,NPTCHAINS*NPTCHAINS
  !! Random permutation of chain indeces:
  idxrand = RANDPERM(NPTCHAINS)
  ic1 = idxrand(1)
  ic2 = idxrand(2)

  CALL RANDOM_NUMBER(ran_uni1)
  betaratio = obj(ic2)%beta-obj(ic1)%beta
  logratio  = betaratio*(obj(ic1)%logL-obj(ic2)%logL) !&

  IF(ran_uni1 <= EXP(logratio))THEN
    !! ACCEPT SWAP
    objtmp1       = obj(ic1)
    objtmp2       = obj(ic2)
    obj(ic1)   = objtmp2
    obj(ic2)   = objtmp1
    !! Temperature does not swap
    obj(ic1)%beta = objtmp1%beta
    obj(ic2)%beta = objtmp2%beta
    !! Acceptance counters do not swap
    obj(ic1)%iaccept_bd = objtmp1%iaccept_bd
    obj(ic2)%iaccept_bd = objtmp2%iaccept_bd
    obj(ic1)%ireject_bd = objtmp1%ireject_bd
    obj(ic2)%ireject_bd = objtmp2%ireject_bd
    ncswap    = ncswap+1_IB
  ENDIF
  ncswapprop    = ncswapprop+1_IB
ENDDO
END SUBROUTINE TEMPSWP_MH
!==============================================================================

SUBROUTINE EXPLORE_MH(obj,dat,cov,beta_mh)
!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,NFPnew,idel,ipar,ipar2,ilay,idxz
TYPE(objstruc)                              :: obj,objnew1,objnew2
TYPE (covstruc),DIMENSION(NRF) :: cov
TYPE(datastruc)                             :: dat
REAL(KIND=RP)                               :: logPLratio,logPQ,logy,logq1_1,logq1_2
REAL(KIND=RP)                               :: ran_uni,ran_uni_BD, ran_unik,ran_uni_ar
REAL(KIND=RP)                               :: znew,beta_mh,Lr1,Lr2
REAL(KIND=RP)                               :: logarp
INTEGER(KIND=IB)                            :: arptype
INTEGER(KIND=IB),DIMENSION(NFPMX)           :: idxrand
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
REAL(KIND=RP),DIMENSION(obj%k)              :: ztmp

objnew1 = obj
!! Draw uniform Birth-Death probability
      CALL RANDOM_NUMBER(ran_uni_BD)

!! Do normal MCMC with 0.5 probability
      IF(ran_uni_BD <= 0.5_RP)THEN
!! Do nothing to k and z
         i_bd = 0
      ENDIF

!! Do BIRTH-DEATH MCMC with 0.5 probability -----------------
      IF(kmin /= kmax)THEN
        IF(ran_uni_BD > 0.5_RP)THEN 
          !! Perturbing k:
         CALL RANDOM_NUMBER(ran_unik)
         i_bd = 0
         IF(obj%k == kmax)THEN  !! If k == kmax, no birth allowed, but 50% death and 50% perturb
              IF(ran_unik>=0.5_RP) i_bd = 2
         ELSEIF(obj%k == kmin)THEN  !! If k == kmin, no death allowed, but 50% birth and 50% perturb
              IF(ran_unik>=0.5_RP) i_bd = 1
         ELSE !! in all other cases it's 33% birth 33% death 33% perturb
             IF((ran_unik <= 1._RP/3._RP)) i_bd = 1
             IF((ran_unik > 1._RP/3._RP) .AND. (ran_unik <= 2._RP/3._RP)) i_bd = 2
         ENDIF

         IF(i_bd == 1)THEN
           CALL BIRTH(obj,objnew1)
if (objnew1%wdepth == 1._RP) then
write(6,*) 'call birth'
 write(6,*) 'old', obj%par, obj%wres, obj%wdepth
 write(6,*) 'new', objnew1%par, objnew1%wres, objnew1%wdepth
endif

         ELSEIF(i_bd == 2)THEN
           CALL DEATH(obj,objnew1)
		   
if (objnew1%wdepth == 1._RP)then
write(6,*) 'call death'
 write(6,*) 'old', obj%par, obj%wres, obj%wdepth
 write(6,*) 'new', objnew1%par, objnew1%wres, objnew1%wdepth
endif
         ENDIF
        ENDIF
      ENDIF
! ---------------------------------------------------------------------
!!
!! Do Metropolis-Hastings
!!
      IF(obj%k /= objnew1%k)THEN
   !!
   !! If k changed, check BD acceptance
   !!
      CALL CHECKBOUNDS(objnew1)
      IF(ioutside == 0)THEN
        CALL LOGLHOOD(objnew1,dat,cov)
        IF(i_bd == 1)THEN
         !! BIRTH ACCEPTANCE RATIO
         logPQ = -objnew1%logPr + objnew1%lognorm + &
           (0.5_RP*DOT_PRODUCT(MATMUL((objnew1%gp-objnew1%g),&
           objnew1%Chati),(objnew1%gp-objnew1%g)))
         logPLratio =  logPQ + (objnew1%logL - obj%logL)*beta_mh
		 
!		 write(6,*) 'birth acceptance ratio', logPLratio
                ELSEIF(i_bd == 2)THEN
                     !! DEATH ACCEPTANCE RATIO
                   logPQ = objnew1%logPr - objnew1%lognorm - &
                       (0.5_RP*DOT_PRODUCT(MATMUL((objnew1%gp-objnew1%g),objnew1%Chati),(objnew1%gp-objnew1%g)))
                   logPLratio = logPQ + (objnew1%logL - obj%logL)*beta_mh  
    !    if (objnew1%k==0) then
    !       write(6,*) '--------------'
    !       write(6,*) 'Death acceptance ratio'
    !       write(6,*) ' objnew1%logPr, objnew1%lognorm,loggp', objnew1%logPr, objnew1%lognorm,-1.*(logPQ-objnew1%logPr+objnew1%lognorm)
    !       write(6,*) ' objnew1%logL, obj%logL,beta_mh', objnew1%logL, obj%logL,beta_mh
    !       write(6,*) 'PQ, PLratio', exp(logPQ), exp(logPLratio) 
    !       write(6,*) '-------------------------------------'
    !    endif
                ENDIF
          CALL RANDOM_NUMBER(ran_uni)
                IF(ran_uni >= EXP(logPLratio))THEN
                    objnew1 = obj
                    obj%ireject_bd = obj%ireject_bd + 1
                ELSE
       if (objnew1%k==0) then
    !  write(6,*) 'b/d happened, rejectence rate:', obj%ireject_bd
    !  write(6,*) 'old parameters:'
    !  write(6,*) obj%par
    !  write(6,*) 'new parameters:'
    !  write(6,*) objnew1%par
    !  write(6,*) 'Death acceptance ratio'
    !  write(6,*) ' objnew1%logPr, objnew1%lognorm,loggp', objnew1%logPr, objnew1%lognorm,-1.*(logPQ-objnew1%logPr+objnew1%lognorm)
    !   write(6,*) ' objnew1%logL, obj%logL,beta_mh', objnew1%logL, obj%logL,beta_mh
    !   write(6,*) 'PQ, PLratio', exp(logPQ), exp(logPLratio) 
    !  write(6,*) '----------------------------------------------------'
       endif
if (objnew1%wres == 1._RP) then
 write(6,*) 'BD accepted'
 write(6,*) 'new,old L', objnew1%logL, obj%logL
 write(6,*) 'old', obj%par, obj%wres
 write(6,*) 'new wres', objnew1%par, objnew1%wres
endif

if (objnew1%wdepth == 1._RP) then
 write(6,*) 'BD accepted'
 write(6,*) 'new,old L', objnew1%logL, obj%logL
 write(6,*) 'old', obj%par, obj%wdepth
 write(6,*) 'new wdepth', objnew1%par, objnew1%wdepth
endif

                    obj = objnew1
	!		write(6,*) 'model accepted'	
	!		write(6,*) obj%par, obj%wres		
					
                    obj%iaccept_bd = obj%iaccept_bd + 1
               !		 stop
                ENDIF 
          ELSE
             objnew1 = obj
             obj%ireject_bd = obj%ireject_bd + 1
             ioutside = 0
          ENDIF
      ELSE  ! k-change if
   !!
   !! If k did not change, carry out MH sweep
   !!
          idxrand = 0
          idxrand(1:obj%NFP) = RANDPERM(obj%NFP)
      !!
      !! Do Metropolis-Hastings update on rho and depth
      !!
      DO ipar = 1,obj%NFP ! sweeps over all parameters
        iwhich = idxrand(ipar) ! selects a parameter out of all possibilities but mixed so there is no bias which parameter is always perturbed first
        CALL PROPOSAL(obj,objnew1,iwhich,1._RP)
        CALL CHECKBOUNDS(objnew1)
        IF(ioutside == 0)THEN
          CALL LOGLHOOD(objnew1,dat,cov)
          logPLratio = (objnew1%logL - obj%logL)*beta_mh
          CALL RANDOM_NUMBER(ran_uni)
          IF(ran_uni >= EXP(logPLratio))THEN
            !!
            !! Try delayed rejection
            !!
            IF(IDELAY == 1)THEN
              objnew2 = obj
              CALL PROPOSAL(obj,objnew2,iwhich,factdelay)
if (objnew1%wdepth == 1._RP)then
write(6,*) 'perturb'
 write(6,*) 'old', obj%par, obj%wres, obj%wdepth
 write(6,*) 'new', objnew2%par, objnew2%wres, objnew2%wdepth
endif
              CALL CHECKBOUNDS(objnew2)
              IF(ioutside == 0)THEN
                ilay  = CEILING(REAL(iwhich,RP)/REAL(NPL,RP))
                ipar2  = iwhich-(ilay-1)*NPL
                IF(ilay > obj%k) ipar2 = ipar2 + 1
                !! Cauchy:
                logq1_1 = LOG((fact*pertsd(ipar2))/((objnew1%pariwhichrot-obj%pariwhichrot)**2._RP    +(fact*pertsd(ipar2))**2._RP))
                logq1_2 = LOG((fact*pertsd(ipar2))/((objnew1%pariwhichrot-objnew2%pariwhichrot)**2._RP+(fact*pertsd(ipar2))**2._RP))

                CALL LOGLHOOD(objnew2,dat,cov)
                Lr1 = EXP(objnew1%logL-objnew2%logL)
                IF(Lr1 > 1._RP)Lr1 = 1._RP
                Lr2 = EXP(objnew1%logL-obj%logL)
                IF(Lr2 > 1._RP)Lr2 = 1._RP
                logPLratio = (objnew2%logL - obj%logL + logq1_2 - logq1_1 &
                             + LOG(1._RP - Lr1) &
                             - LOG(1._RP - Lr2))*beta_mh
                CALL RANDOM_NUMBER(ran_uni)
                IF(ran_uni >= EXP(logPLratio))THEN
                  objnew1 = obj
                  ireject = ireject + 1
                  ireject_delay = ireject_delay + 1
                ELSE
if (objnew2%wres == 1._RP) then
 write(6,*) 'parameter change accepted'
 write(6,*) 'new,old L', objnew1%logL, obj%logL
 write(6,*) 'old', obj%par, obj%wres
 write(6,*) 'new wres', objnew2%par, objnew2%wres
endif
if (objnew2%wdepth == 1._RP) then
 write(6,*) 'parameter change accepted'
 write(6,*) 'new,old L', objnew1%logL, obj%logL
 write(6,*) 'old', obj%par, obj%wres
 write(6,*) 'new wdepth', objnew2%par, objnew2%wdepth
endif
                  obj = objnew2
                  iaccept = iaccept + 1
                  iaccept_delay = iaccept_delay + 1
                ENDIF
              ELSE
                objnew1 = obj
                ireject = ireject + 1
                ioutside = 0
              ENDIF
            ELSE  !! No delayed rejection, just reject
              objnew1 = obj
              ireject = ireject + 1
            ENDIF
          ELSE
            obj = objnew1
	!		write(6,*) 'accepted'
	!             write(6,*) objnew1%par,objnew1%wres 		
            iaccept = iaccept + 1
          ENDIF
        ELSE !! outside before delayed rejection
          objnew1 = obj
          ireject = ireject + 1
          ioutside = 0
        ENDIF
      ENDDO
   !!
   !! Do Metropolis-Hastings on water resistivity and water depth
   !!
      IF(IWRES == 1)THEN
      !! Perturb water resistivity
      !! CAUCHY proposal
       CALL RANDOM_NUMBER(ran_uni)
       objnew1%wres = obj%wres + &
         pertwressd*TAN(PI2*(ran_uni-0.5_RP))
       IF(((objnew1%wres - minlimwres) < 0._RP)&
	     .OR.((maxlimwres - objnew1%wres) < 0._RP)) ioutside = 1
       IF(ioutside == 0)THEN
!	   write(6,*) 'new water resistivity', objnew1%wres
if (objnew1%wdepth == 1._RP)then
write(6,*) 'call death'
 write(6,*) 'old', obj%par, obj%wres, obj%wdepth
 write(6,*) 'new', objnew1%par, objnew1%wres, objnew1%wdepth
endif
          CALL LOGLHOOD(objnew1,dat,cov)
!	   write(6,*) 'likelihood new, old', objnew1%logL, obj%logL
            logPLratio = (objnew1%logL - obj%logL)*beta_mh
          CALL RANDOM_NUMBER(ran_uni)
          IF(ran_uni >= EXP(logPLratio))THEN
             objnew1 = obj
             ireject = ireject + 1
          ELSE
		  
if (objnew1%wres == 1._RP) then
 write(6,*) 'water resistivity accepted'
 write(6,*) 'new,old L, oldwres', objnew1%logL, obj%logL, obj%wres
endif
if (objnew1%wdepth == 1._RP) then
 write(6,*) 'water resistivity accepted'
 write(6,*) 'new,old L, oldwdepth', objnew1%logL, obj%logL, obj%wdepth
endif
             obj = objnew1
             iaccept = iaccept + 1
	!		 write(6,*) 'water resistivity accepted'
     !        write(6,*) objnew1%par,objnew1%wres 
          ENDIF
       else
          objnew1 = obj
          ireject = ireject + 1
          ioutside = 0
       endif
       endif ! iwres
	   
      IF(IWDEPTH == 1)THEN
      !! Perturb water depth
      !! CAUCHY proposal
       CALL RANDOM_NUMBER(ran_uni)
       objnew1%wdepth = obj%wdepth + &
         pertwdepthsd*TAN(PI2*(ran_uni-0.5_RP))
       IF(((objnew1%wdepth - minlimwdepth) < 0._RP)&
	     .OR.((maxlimwdepth - objnew1%wdepth) < 0._RP)) ioutside = 1
       IF(ioutside == 0)THEN
          CALL LOGLHOOD(objnew1,dat,cov)
             logPLratio = (objnew1%logL - obj%logL)*beta_mh
          CALL RANDOM_NUMBER(ran_uni)
          IF(ran_uni >= EXP(logPLratio))THEN
             objnew1 = obj
             ireject = ireject + 1
          ELSE
             obj = objnew1
             iaccept = iaccept + 1
          ENDIF
       else
          objnew1 = obj
          ireject = ireject + 1
          ioutside = 0
       endif
       endif ! iwdepth
   
      IF(IDT == 1)THEN
       DO ipar = 1,NRF
      !! Perturb time shift
      !! CAUCHY proposal
         CALL RANDOM_NUMBER(ran_uni)
          objnew1%dt(ipar) = obj%dt(ipar) + &
            pertdtsd*TAN(PI2*(ran_uni-0.5_RP))
         IF(((timebuf(ipar,1)+objnew1%dt(ipar)) < 0._RP)&
           .OR.((minlimdt-objnew1%dt(ipar)) > 0._RP)&
	       .OR.((maxlimdt - objnew1%dt(ipar)) < 0._RP)) ioutside = 1
       
         IF(ioutside == 0)THEN
          CALL LOGLHOOD(objnew1,dat,cov)
             logPLratio = (objnew1%logL - obj%logL)*beta_mh
          CALL RANDOM_NUMBER(ran_uni)
          IF(ran_uni >= EXP(logPLratio))THEN
             objnew1 = obj
             ireject = ireject + 1
          ELSE
             obj = objnew1
             iaccept = iaccept + 1
          ENDIF
         else
          objnew1 = obj
          ireject = ireject + 1
          ioutside = 0
         endif
       enddo
      elseif(IDT == 2)THEN
        CALL RANDOM_NUMBER(ran_uni)
           objnew1%dt(1) = obj%dt(1) + &
            pertdtsd*TAN(PI2*(ran_uni-0.5_RP))
         IF(((timebuf(1,1)+objnew1%dt(1)) < 0._RP)&
           .OR.((minlimdt-objnew1%dt(1)) > 0._RP)&
	       .OR.((maxlimdt - objnew1%dt(1)) < 0._RP)) ioutside = 1
         if (NRF >= 2) then
           DO ipar = 2,NRF
           objnew1%dt(ipar)=objnew1%dt(1)
           IF(((timebuf(ipar,1)+objnew1%dt(ipar)) < 0._RP)&
             .OR.((minlimdt-objnew1%dt(ipar)) > 0._RP)&
	            .OR.((maxlimdt - objnew1%dt(ipar)) < 0._RP)) ioutside = 1
           enddo
         endif
         IF(ioutside == 0)THEN
          CALL LOGLHOOD(objnew1,dat,cov)
             logPLratio = (objnew1%logL - obj%logL)*beta_mh
          CALL RANDOM_NUMBER(ran_uni)
          IF(ran_uni >= EXP(logPLratio))THEN
             objnew1 = obj
             ireject = ireject + 1
          ELSE
             obj = objnew1
             iaccept = iaccept + 1
          ENDIF
         else
          objnew1 = obj
          ireject = ireject + 1
          ioutside = 0
         endif
	  
      endif ! idt   
	   
   !!
   !! Do Metropolis-Hastings on data-error standard deviations
   !!
      IF(ICOV == 1 .OR. ICOV == 5)THEN
      !! Perturb std devs with .25 probability
       CALL RANDOM_NUMBER(ran_uni_ar)
       IF(ran_uni_ar>=0.25_RP)THEN
         DO ipar = 1,NRF
            CALL PROPOSAL_SD(obj,objnew1,ipar)
            IF(ioutside == 0)THEN
               CALL LOGLHOOD(objnew1,dat,cov)
               logPLratio = (objnew1%logL - obj%logL)*beta_mh
               CALL RANDOM_NUMBER(ran_uni)
               IF(ran_uni >= EXP(logPLratio))THEN
                  objnew1 = obj
                  ireject = ireject + 1
               ELSE
                  obj = objnew1
                  iaccept = iaccept + 1
               ENDIF
            ELSE
               objnew1 = obj
               ireject = ireject + 1
               ioutside = 0
            ENDIF
            i_sdpert = 0
         ENDDO
       ENDIF
      ENDIF ! ICOV if
   !!
   !! Do Metropolis-Hastings on cf
   !!
   IF(ICF >= 1)THEN
      !! Perturb std devs with .25 probability
      CALL RANDOM_NUMBER(ran_uni_ar)
      IF(ran_uni_ar>=0.25_RP)THEN
         DO ipar = 1,NRF-(ICF-1)
            CALL PROPOSAL_CF(obj,objnew1,ipar)
            IF(ioutside == 0)THEN
               CALL LOGLHOOD(objnew1,dat,cov)
               logPLratio = (objnew1%logL - obj%logL)*beta_mh
               CALL RANDOM_NUMBER(ran_uni)
               IF(ran_uni >= EXP(logPLratio))THEN
                  objnew1 = obj
                  ireject = ireject + 1
               ELSE
                  obj = objnew1
                  iaccept = iaccept + 1
               ENDIF
            ELSE
               objnew1 = obj
               ireject = ireject + 1
               ioutside = 0
            ENDIF
            i_sdpert = 0
         ENDDO
      ENDIF
   ENDIF ! ICOV if

   !!
   !! Do Metropolis-Hastings on autoregressive model
   !!
   IF(IAR == 1)THEN
      !! Perturb AR model with .25 probability
      CALL RANDOM_NUMBER(ran_uni_ar)
      IF(ran_uni_ar>=0.25_RP)THEN
         DO ipar = 1,NARFP*NRF
            IF(obj%idxar(ipar) == 0)THEN
              !! Propose birth
              arptype = 1
              logarp = LOG(0.5_RP)
            ELSE
              CALL RANDOM_NUMBER(ran_uni_ar)
              IF(ran_uni_ar>=0.5_RP)THEN
                !! Propose death
                arptype = 2
                logarp = LOG(2._RP)
              ELSE
                !! Propose perturb
                arptype = 3
                logarp = 0._RP
              ENDIF
            ENDIF
            CALL PROPOSAL_AR(obj,objnew1,ipar,arptype)
            IF(ioutside == 0)THEN
               CALL LOGLHOOD(objnew1,dat,cov)
               !! logPLratio = (objnew1%logL - obj%logL)*beta_mh
               !! Input Birth Death AR here:
               logPLratio = logarp + (objnew1%logL - obj%logL)*beta_mh
               CALL RANDOM_NUMBER(ran_uni)
               IF(ran_uni >= EXP(logPLratio))THEN
                  objnew1 = obj
                  ireject = ireject + 1
               ELSE
                  obj = objnew1
                  iaccept = iaccept + 1
               ENDIF
            ELSE
               objnew1 = obj
               ireject = ireject + 1
               ioutside = 0
            ENDIF
            i_sdpert = 0
         ENDDO
      ENDIF
   ENDIF ! AR if
ENDIF ! k-change if


END SUBROUTINE EXPLORE_MH
!=======================================================================

SUBROUTINE CHECKBOUNDS(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                            :: ip,iwhich,ipar,ilay,ncra
INTEGER(KIND=IB)                            :: ih
TYPE(objstruc)                              :: obj
REAL(KIND=RP)                               :: vspmin,vspmax,zhere
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
REAL(KIND=RP),DIMENSION(obj%k)              :: ztmp

      if (obj%k==0) then
             IF(obj%par((obj%k-1)*NPL+1) > maxlim(2))ioutside = 1
             IF(obj%par(1) < minlim(2)) ioutside = 1
      else
             IF(obj%par((obj%k-1)*NPL+1) > maxlim(1))ioutside = 1
             IF(obj%par(1) < minlim(1)) ioutside = 1
      endif

      DO ilay = 1,obj%k
       IF(hmin > obj%h(ilay))ioutside = 1
       IF(maxlim(1) < obj%h(ilay))ioutside = 1
      ENDDO

      CALL GETIDXCRA(obj,idxcra,ncra)
      DO ip = 1,ncra
       iwhich = idxcra(ip)
       ilay = CEILING(REAL(iwhich,RP)/REAL(NPL,RP))
       ipar = iwhich-(ilay-1)*NPL
      IF(ilay > obj%k) ipar = ipar + 1
        IF(((obj%par(iwhich) - minlim(ipar)) < 0._RP).OR.((maxlim(ipar) - obj%par(iwhich)) < 0._RP))THEN
        ioutside = 1
      !PRINT*,iwhich,obj%par(iwhich), minlim(ipar),maxlim(ipar)
        ENDIF
      ENDDO 

!DO ip = 1,NMISC
!  IF(((obj%misc(ip) - minlimmisc(ip)) < 0._RP).OR.((maxlimmisc(ip) - obj%misc(ip)) < 0._RP))THEN
!    ioutside = 1
!    !PRINT*,obj%misc(ip), minlimmisc(ip), maxlimmisc(ip)
!  ENDIF
!ENDDO
RETURN
END SUBROUTINE CHECKBOUNDS
!=======================================================================

SUBROUTINE PROPOSAL(obj,objnew,iwhich,factor)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i,iwhich,ipar,ipar2,iloop,ilay,ilay2
TYPE(objstruc) :: obj,objnew,objtmp
REAL(KIND=RP)  :: ran_uni, ran_nor, factor
!REAL(KIND=RP),DIMENSION(obj%NFP) :: m2p,m1
!REAL(KIND=RP),DIMENSION(obj%NFP,obj%NFP) :: VVtmp

!VVtmp = 0._RP
!VVtmp = VV(1:obj%NFP,1:obj%NFP,obj%k)

      objnew = obj
      ilay  = CEILING(REAL(iwhich,RP)/REAL(NPL,RP))
      ipar  = iwhich-(ilay-1)*NPL
      IF(ilay > obj%k) ipar = ipar + 1
!!
!! Gaussian proposal
!!
      CALL GASDEVJ(ran_nor)
!objnew%par(iwhich) = obj%par(iwhich) + pertsd(ipar)*ran_nor
!!
!! CAUCHY proposal
!!
       CALL RANDOM_NUMBER(ran_uni)

!IF(ILINROT == 0)THEN
  !! Gaussian:
!  objnew%par(iwhich) = obj%par(iwhich) + fact/factor*pertsd(ipar)*ran_nor
  !! Cauchy:
      objnew%par(iwhich) = obj%par(iwhich) + &
                 fact/factor*pertsd(ipar)*TAN(PI2*(ran_uni-0.5_RP))
  !! Save the parameters for delayed rejection
      obj%pariwhichrot = obj%par(iwhich)
      objnew%pariwhichrot = objnew%par(iwhich)

      CALL CALCH(objnew)

      RETURN
      END SUBROUTINE PROPOSAL
!=======================================================================

      SUBROUTINE PROPOSAL_AR(obj,objnew,iwhich,arptype)   
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,arptype
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor,ran_uni

!! Birth: sample uniform from prior
      IF(arptype == 1)THEN
       CALL RANDOM_NUMBER(ran_uni)
       objnew%arpar(iwhich) = ran_uni*(maxlimar(iwhich)-minlimar(iwhich))+minlimar(iwhich)
       objnew%idxar(iwhich) = 1
      IF(((objnew%arpar(iwhich) - minlimar(iwhich)) < 0._RP).OR. &
         ((maxlimar(iwhich) - objnew%arpar(iwhich)) < 0._RP))ioutside = 1
      ENDIF
!! Death
      IF(arptype == 2)THEN
       objnew%arpar(iwhich) = minlimar(iwhich)-1._RP
       objnew%idxar(iwhich) = 0
      ENDIF
!! Perturb
      IF(arptype == 3)THEN
       CALL GASDEVJ(ran_nor)
       objnew%arpar(iwhich) = obj%arpar(iwhich) + pertarsd(iwhich)*ran_nor
       IF(((objnew%arpar(iwhich) - minlimar(iwhich)) < 0._RP).OR. &
         ((maxlimar(iwhich) - objnew%arpar(iwhich)) < 0._RP))ioutside = 1
      ENDIF

      RETURN
      END SUBROUTINE PROPOSAL_AR
!!=======================================================================
!
!SUBROUTINE PROPOSAL_MISC(obj,objnew,iwhich)
!!=======================================================================
!USE DATA_TYPE
!USE RJMCMC_COM
!IMPLICIT NONE
!INTEGER(KIND=IB) :: iwhich
!TYPE(objstruc) :: obj,objnew
!REAL(KIND=RP)  :: ran_nor
!
!!
!! Gaussian proposal
!!
!CALL GASDEVJ(ran_nor)
!objnew%misc(iwhich) = obj%misc(iwhich) + pertmiscsd(iwhich)*ran_nor
!IF(((objnew%misc(iwhich) - minlimmisc(iwhich)) < 0._RP).OR.((maxlimmisc(iwhich) - objnew%misc(iwhich)) < 0._RP))ioutside = 1
!
!RETURN
!END SUBROUTINE PROPOSAL_MISC
!=======================================================================

SUBROUTINE CHECKBOUNDS_ARMX(dat)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ifr
TYPE(datastruc):: dat

DO ifr = 1,NRF
!   IF(MAXVAL(dat%dar(ifr,:)) > MAXVAL(dat%Dobs(ifr,:)))THEN
!      iarfail = 1
!   ENDIF
!   IF(MINVAL(dat%dar(ifr,:)) < MINVAL(dat%Dobs(ifr,:)))THEN
!      iarfail = 1
!   ENDIF
   IF(MAXVAL(dat%dar(ifr,:)) > armx)THEN
      iarfail = 1
   ENDIF
   IF(MINVAL(dat%dar(ifr,:)) < -armx)THEN
      iarfail = 1
   ENDIF
ENDDO

RETURN
END SUBROUTINE CHECKBOUNDS_ARMX
!=======================================================================

SUBROUTINE PROPOSAL_SD(obj,objnew,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
      CALL GASDEVJ(ran_nor)
objnew%sdpar(iwhich) = obj%sdpar(iwhich) + pertsdsd(iwhich)*ran_nor
      IF(((objnew%sdpar(iwhich) - minlimsd(iwhich)) < 0._RP)&
	     .OR.((maxlimsd(iwhich) - objnew%sdpar(iwhich)) < 0._RP)) ioutside = 1

RETURN
END SUBROUTINE PROPOSAL_SD
!=======================================================================

SUBROUTINE PROPOSAL_CF(obj,objnew,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
      CALL GASDEVJ(ran_nor)
        objnew%cf(iwhich) = obj%cf(iwhich) + pertcfsd(iwhich)*ran_nor
      IF(((objnew%cf(iwhich) - minlimcf(iwhich)) < 0._RP).OR.((maxlimcf(iwhich) - objnew%cf(iwhich)) < 0._RP))ioutside = 1

      RETURN
      END SUBROUTINE PROPOSAL_CF
!=======================================================================

SUBROUTINE GETIDXCRA(obj,idxcra,ncra)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ipar,ilay,ipl,isd,ncra
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)) :: idxcra

IF(obj%k == 0)THEN
  !idxcra = (/1:NPL-1:1/)
  DO ipar = 1,NPL-1
     idxcra(ipar) = ipar
  ENDDO
  ncra = NPL-1
ELSE
  !! Layers:
  idxcra = 0
  ipar = 1
  DO ilay = 1,obj%k
    DO ipl = 2,NPL
      idxcra(ipar) = (ilay-1)*NPL+ipl
      ipar = ipar + 1
    ENDDO
  ENDDO
  !! Half-space:
  DO ipl = 1,NPL-1
    idxcra(ipar) = obj%k*NPL+ipl
    ipar = ipar + 1
  ENDDO
  ncra = (obj%k+1)*(NPL-1)
ENDIF

END SUBROUTINE GETIDXCRA
!=======================================================================

SUBROUTINE CALCH(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i,j
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION(obj%k) :: idxz
REAL(KIND=IB),DIMENSION(obj%k) :: h

IF(obj%k == 0)THEN
  obj%h = 0._RP
  obj%z = 0._RP
ELSE
  idxz = 0
  !idxz = (/1:obj%k*NPL:NPL/)
  j = 1
  DO i = 1,obj%k*NPL,NPL
    idxz(j) = i
    j = j+1
  ENDDO
  obj%z = 0._RP
  obj%z(1:obj%k) = obj%par(idxz)
  obj%h = 0._RP
  obj%h(1) = obj%par(1)
  DO i = 2,obj%k
    obj%h(i) = obj%par(idxz(i))-obj%par(idxz(i-1))
  ENDDO
ENDIF
END SUBROUTINE CALCH
!=======================================================================

SUBROUTINE CALCM(obj,mtmp)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: i,j
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION(obj%k) :: idxh
REAL(KIND=IB),DIMENSION(NFPMX) :: mtmp

mtmp = 0._RP
mtmp = obj%par
idxh = 0
!idxh = (/1:obj%k*NPL:NPL/)
j = 1
DO i = 1,obj%k*NPL,NPL
  idxh(j) = i
  j = j+1
ENDDO
mtmp(idxh) = obj%h(1:obj%k)

END SUBROUTINE CALCM
!=======================================================================
SUBROUTINE DEATH(obj,objnew)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                 :: idel,idel2
TYPE(objstruc)                   :: obj,objnew
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxdeath
REAL(KIND=RP)                    :: ran_uni
REAL(KIND=RP),DIMENSION(NPL-1)   :: cra

objnew%k   = obj%k - 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)
!! Pick random layer:
idxdeath = 0
idxdeath(1:obj%k) = RANDPERM(obj%k)
idel = idxdeath(1)

!! Delete interface (setting new layer properties to one of old layers
!! properties at random ):
      CALL RANDOM_NUMBER(ran_uni)
      idel2 = idel
      IF(ran_uni>=0.5_RP) idel2 = idel+1

      objnew%par = 0._RP  ! Ensure no "left overs" from other dimensions exist
      cra = 0._RP
      IF(idel == 1)THEN
   !! Killed interface is first interface:
         !! no interface left
         if (objnew%k==0) then
         IF(idel2 == 1) objnew%par(1:NPL-1)=obj%par(2:NPL)
         IF(idel2 == 2) objnew%par(1:NPL-1)=obj%par(1+NPL:obj%NFP)  ! parameter 2 or 3 = resistivity 1 or 2    
            objnew%g  = obj%par(2:NPL)
            objnew%gp = obj%par(1+NPL:obj%NFP)
		!	write(6,*) 'k==0, idel2,obj%par,objnew%par', idel2,obj%par,objnew%par
         else
            cra = obj%par((idel2-1)*NPL+2:(idel2-1)*NPL+NPL)
            objnew%par(1:objnew%NFP) = obj%par(NPL+1:obj%NFP)
            objnew%par(2:NPL) = cra ! new resistivity is either the resistivity of the 1st or 2nd layer that now combine
   !! This records the perturbation for the bd acceptance ratio
       objnew%g  = obj%par(2:NPL)
       objnew%gp = obj%par(2+NPL:2*NPL)
	     endif
      ELSEIF(idel == obj%k)THEN
   !! Killed interface is last interface:
         IF(idel2 == obj%k) cra = obj%par((idel2-1)*NPL+2:(idel2-1)*NPL+NPL)
         IF(idel2 > obj%k) cra = obj%par(obj%NFP-(NPL-2):obj%NFP)
         objnew%par(1:objnew%NFP-(NPL-1)) = (/ obj%par(1:obj%NFP-(NPL*2-1)) /)
         objnew%par(objnew%NFP-(NPL-2):objnew%NFP) = cra
   !! This records the perturbation for the bd acceptance ratio
         objnew%g  = obj%par((idel-1)*NPL+2:(idel-1)*NPL+NPL)
         objnew%gp = obj%par(obj%NFP-(NPL-2):obj%NFP) 
      ELSE
   !! Killed interface is in normal layer stack:
        cra = obj%par((idel2-1)*NPL+2:(idel2-1)*NPL+NPL)
        objnew%par(1:objnew%NFP) = (/ obj%par(1:(idel-1)*NPL),obj%par(idel*NPL+1:obj%NFP) /)
        objnew%par((idel-1)*NPL+2:(idel-1)*NPL+NPL) = cra
   !! This records the perturbation for the bd acceptance ratio
        objnew%g  = obj%par((idel-1)*NPL+2:(idel-1)*NPL+NPL)
        objnew%gp = obj%par((idel)*NPL+2:(idel)*NPL+NPL)
       ENDIF
      CALL CALCH(objnew)

      END SUBROUTINE DEATH
!=======================================================================

SUBROUTINE BIRTH(obj,objnew)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)              :: i, iznew, ipert, iloop, iwhich, ipar
TYPE(objstruc)                :: obj,objnew,objnew2
REAL(KIND=RP)                 :: znew,hnew1,hnew2,ran_uni,ran_nor
REAL(KIND=RP),DIMENSION(obj%k):: ztmp

objnew%k   = obj%k + 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)
!!
!! Draw new z until new layer > hmin
!!
      CALL RANDOM_NUMBER(ran_uni)
         znew = maxpert(1)*ran_uni
       if (obj%k==0) then
         ztmp = hmx - znew
!		 write(6,*) 'k and znew', obj%k,znew
       else
         ztmp = obj%z(1:obj%k) - znew
       endif
      IF(hmin > MINVAL(ABS(ztmp))) ioutside = 1
      IF(hmin > znew) ioutside = 1
!!
!! Insert new interface
!!
      if (obj%k==0) then
      iznew=1 
!	  write(6,*) 'birth of interface'
      else
        iznew = 0
        DO i = 1,obj%k
         IF(ztmp(i) > 0._RP) EXIT
          iznew = i
        ENDDO
        iznew = iznew + 1
      endif
        objnew%z   = 0._RP
        objnew%par = 0._RP
      IF(iznew == 1)THEN
!!   New layer is first layer:
!   hnew1 = znew
!   hnew2 = obj%z(iznew)-znew
       if (objnew%k==1) then
         objnew%par(1)=znew
         objnew%par(2:NPL*2-1)=(/ obj%par(1:NPL-1),obj%par(1:NPL-1) /)
         !write(6,*) 'new model k and par', objnew%k,objnew%par
		  
       else
         objnew%par(1:objnew%NFP) = (/ obj%par(1:NPL),obj%par /)
         objnew%par(1)            = znew
       endif
!   objnew%par(NPL+1)        = hnew2
      ELSEIF(iznew > obj%k)THEN
   !! New layer is created in half-space:
!   hnew1 = znew - obj%z(iznew-1)
      objnew%par(1:objnew%NFP) = (/ obj%par(1:obj%NFP-(NPL-1)),hnew1, &
                                 obj%par(obj%NFP-(NPL-2):obj%NFP), &
                                 obj%par(obj%NFP-(NPL-2):obj%NFP) /)
      objnew%par(objnew%NFP-((NPL*2)-2)) = znew
      ELSE
   !! New layer is within layer stack:
!   hnew1 = znew - obj%z(iznew-1)
!   hnew2 = obj%z(iznew)-znew
      objnew%par(1:objnew%NFP) = (/ obj%par(1:(iznew)*NPL), &
                                 obj%par(((iznew-1)*NPL)+1:iznew*NPL), &
                                 obj%par(iznew*NPL+1:obj%NFP) /)
      objnew%par(((iznew-1)*NPL)+1) = znew
!   objnew%par((iznew*NPL)+1)     = hnew2
      ENDIF
      CALL CALCH(objnew)
!!
!! Pick one of the two new layers at random and perturb:
!!
    CALL RANDOM_NUMBER(ran_uni)
IF(ran_uni <= 0.5_RP)THEN
   ipert = iznew
ELSE
   ipert = iznew+1
ENDIF

!! Perturb only one parameter (determined in EXPLORE_MH call)
objnew2 = objnew
      DO ipar = 2,NPL
   iloop = 0
   iwhich = (ipert-1)*NPL+ipar
       IF(ipert > objnew%k) iwhich = iwhich - 1
   !!
   !! Gaussian proposal
   !!
      CALL GASDEVJ(ran_nor)
       objnew2%par(iwhich) = objnew%par(iwhich) + pertsd(ipar)*ran_nor

   !! This records the perturbation for the bd acceptance ratio
      objnew2%g(ipar-1)  = objnew%par(iwhich)
      objnew2%gp(ipar-1) = objnew2%par(iwhich)

      IF(((objnew2%par(iwhich) - minlim(ipar)) < 0._RP).OR. &
        ((maxlim(ipar) - objnew%par(iwhich)) < 0._RP)) ioutside = 1
      ENDDO
      objnew = objnew2
END SUBROUTINE BIRTH
!=======================================================================

SUBROUTINE ARFB(obj,dat,irf,idata,idatb)
!=======================================================================
!!
!! Autoregressive model to model data error correlations.
!! AR process is computed forward and backward and average is used.
!!
! irf=receiver no, idatb=no of data for this receiver, idata=1
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj
TYPE (datastruc) :: dat
INTEGER          :: i,j,k,irf,idata,idatb
REAL(KIND=RP),DIMENSION(idatb-idata+1)::dres1,dres2,dar1,dar2,dar3

      IF(obj%idxar(irf) == 1)THEN
        k = NARFP                         ! order==1
        dat%dar(irf,idata)=0._RP          ! Matlab sets first point to zero...

       dres1 = 0._RP
       dres1 = dat%res(irf,idata:idatb)

       dar1(1)=0._RP          ! Matlab sets first point to zero...
       DO i=2,idatb-idata+1
        dar1(i) = 0
        IF(k >= i)THEN
          DO j=1,i-1
           dar1(i) = dar1(i) + obj%arpar(NARFP*(irf-1)+j) * dres1(i-j)
          ENDDO
        ELSE
          DO j=1,k
           dar1(i) = dar1(i) + obj%arpar(NARFP*(irf-1)+j) * dres1(i-j)
          ENDDO
        ENDIF
       ENDDO
           dat%dar(irf,idata:idatb) = dar1
           dat%dar(irf,idata) = 0._RP
           dat%dar(irf,idatb) = 0._RP
      ENDIF

END SUBROUTINE ARFB
!=======================================================================

SUBROUTINE PARALLEL_SEED()
!!
!!  Ensure unique random seed for each CPU
!!
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: i
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed1
REAL(KIND=RP) :: ran_uni


      CALL RANDOM_SEED
      CALL RANDOM_SEED(SIZE=iseedsize)
ALLOCATE( iseed1(iseedsize) )
IF(ISETSEED == 1)THEN
   iseed1 = (/2303055,     2435432,     5604058,     4289794,     3472290, &
      7717070,      141180,     3783525,     3087889,     4812786,     3028075, &
      3712062,     6316731,      436800,     7957708,     2055697,     1944360, &
      1222992,     7537775,     7769874,     5588112,     7590383,     1426393, &
      1753301,     7681841,     2842400,     4411488,     7304010,      497639, &
      4978920,     5345495,      754842,     7360599,     5776102/)
   CALL RANDOM_SEED(PUT=iseed1)
ELSE
   CALL RANDOM_SEED(GET=iseed1)
!   IF(rank==src)WRITE(6,*) 'Master seed:',iseed1
ENDIF

ALLOCATE( iseed(iseedsize), rseeds(iseedsize,NTHREAD), iseeds(iseedsize,NTHREAD) )
iseed = 0
rseeds = 0._RP
iseeds = 0
IF(rank == src)THEN
   CALL RANDOM_NUMBER(rseeds)
   iseeds = -NINT(rseeds*1000000._RP)
ENDIF
DO i = 1,iseedsize
   CALL MPI_BCAST( iseeds(i,:), NTHREAD, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )
ENDDO
iseed = iseeds(:,rank+1)

!!
!! Write seeds to seed logfile:
!!
IF(rank == src)THEN
   OPEN(UNIT=50,FILE=seedfile,FORM='formatted',STATUS='UNKNOWN', &
   ACTION='WRITE',POSITION='REWIND',RECL=1024)
   WRITE(50,*) 'Rank: ',rank
   WRITE(50,201) iseed
   WRITE(50,*) ''
   CLOSE(50)
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
DO i = 1,NTHREAD-1
   IF(rank == i)THEN
      OPEN(UNIT=50,FILE=seedfile,FORM='formatted',STATUS='UNKNOWN', &
      ACTION='WRITE',POSITION='APPEND',RECL=1024)
      WRITE(50,*) 'Rank: ',rank
      WRITE(50,201) iseed
      WRITE(50,*) ''
      CLOSE(50)
   ENDIF
   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
ENDDO
CALL RANDOM_SEED(PUT=iseed)

DO i = 1,100
   CALL RANDOM_NUMBER(ran_uni)
ENDDO
DO i = 1,CEILING(ran_uni*10000._RP)
   CALL RANDOM_NUMBER(ran_uni)
ENDDO

201   FORMAT(50I10)
END SUBROUTINE PARALLEL_SEED
!=======================================================================

SUBROUTINE SAVESAMPLE(logLG,isource,tcmp)
!=======================================================================
!!
!! Exchanging and saving posterior samples
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: i,j,ikeep,isource
REAL(KIND=RP)   :: logLG,tcmp
REAL(KIND=RP)   :: t1,t2

!!
!!  Sending samples to master
!!
IF(rank == src)THEN
   ikeep = NKEEP
   IF(iconv == 1)THEN
      iconv2 = 1
      iconv3 = iconv3 + 1
   ENDIF
   tag   = MPI_ANY_TAG
   buf_save_rcv = 0._RP

   CALL MPI_RECV(buf_save_rcv, ncount, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,&
                 tag, MPI_COMM_WORLD, status, ierr )
   isource = status(MPI_SOURCE)
   CALL MPI_RECV(tcmp, 1, MPI_DOUBLE_PRECISION, isource,&
                 tag, MPI_COMM_WORLD, status, ierr )

   isize1 = SIZE(buffer1,1)*IB
   CALL MPI_BUFFER_ATTACH(buffer1,isize1,ierr)
   CALL MPI_BSEND(iconv,1,MPI_INTEGER,isource,rank,MPI_COMM_WORLD,ierr)
   CALL MPI_BUFFER_DETACH(buffer1,isize1,ierr)

   sample = RESHAPE(buf_save_rcv,(/ NKEEP,ncount2 /))


   !!
   !! Master writes sample to file
   !!
   DO j=1,NKEEP
! write(*,*) 'ikeep',j
! write(*,*) j,sample(j,:)
      WRITE(40,207) sample(j,:)
   ENDDO
   CALL FLUSH(40)
   ikeep = ikeep + NKEEP
   logLG = MAXVAL(sample(:,1))
ELSE
   buf_save_snd = 0._RP
   buf_save_snd = PACK(sample(1:NKEEP,:),.true.)

   isize2 = SIZE(buffer2,1)*RP
   CALL MPI_BUFFER_ATTACH(buffer2,isize2,ierr)
   CALL MPI_BSEND( buf_save_snd, ncount, MPI_DOUBLE_PRECISION, src, &
                   rank, MPI_COMM_WORLD, ierr )
   CALL MPI_BUFFER_DETACH(buffer2,isize2,ierr)

   isize3 = size(buffer3,1)*RP
   CALL MPI_BUFFER_ATTACH(buffer3,isize3,ierr)
   CALL MPI_BSEND( tcmp, 1, MPI_DOUBLE_PRECISION, src, &
                   rank, MPI_COMM_WORLD, ierr )
   CALL MPI_BUFFER_DETACH(buffer3,isize3,ierr)

   CALL MPI_RECV(iconv,1,MPI_INTEGER,src,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
ENDIF !  MPI
CALL FLUSH(6)
207   FORMAT(500ES18.8)
RETURN
END SUBROUTINE SAVESAMPLE
!=======================================================================

SUBROUTINE READCOVI(cov,obj,dat)
!!==============================================================================
!!
!! Reads covariance matrix
!!
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB) :: idat,istn,irow,irf
TYPE(covstruc),DIMENSION(NRF) :: cov  ! Structure for inverse data cov mat
TYPE (objstruc) :: obj
TYPE (datastruc) :: dat
REAL(KIND=RP),DIMENSION(NRF)   :: factsdpar
character(len=20) :: x1,covfile

      if (icov == 2) then
       IF(rank == src) WRITE(6,*) ' Reading factors for standard deviation from file '
         OPEN(UNIT=20,FILE='fort.44',FORM='formatted',STATUS='OLD',ACTION='READ')
           do irf=1,NRF
             READ(20,203) factsdpar(irf)
    IF(rank == src) write(6,*) 'stdev factor ',irf, factsdpar(irf)
           enddo
         close(20)
      endif
		
      DO istn = 1,NRF 
       write(x1,686) istn
        covfile=TRIM(ADJUSTL(icovfile))//trim(x1)//'.dat'

       OPEN(UNIT=20,FILE=covfile,FORM='formatted',STATUS='OLD',ACTION='READ')
        ALLOCATE( cov(istn)%Cdi(dat%NDAT(istn),dat%NDAT(istn)) )
         cov(istn)%Cdi = 0._RP
	!     cov(istn)%logdet=0._RP
       
        DO irow = 1,dat%NDAT(istn)
          if (icov>=4) then
          READ(20,676) cov(istn)%Cdi(irow,irow)
!  		 cov(istn)%Cdi(irow,irow)=cov(istn)%Cdi(1,1)
     !      cov(istn)%logdet=cov(istn)%logdet-log(cov(istn)%Cdi(irow,irow))
          else
             DO idat = 1,dat%NDAT(istn)
               READ(20,676) cov(istn)%Cdi(irow,idat)
               if (icov==2)  cov(istn)%Cdi(irow,idat)=cov(istn)%Cdi(irow,idat)/(factsdpar(istn)**2)
             ENDDO
          endif
        ENDDO
		
      CLOSE(20)
	!  write(6,*) 'Covariance estimate:',istn,cov(istn)%Cdi(1:5,1:5)
	 
      ENDDO
686    format(I1)
676	   format(e24.16)
203    FORMAT(24F12.4)
RETURN
END SUBROUTINE READCOVI
!=======================================================================

SUBROUTINE COUPLE_CR(obj)
!=======================================================================
!!
!! Constrain to Hamilton bounds THIS IS FOR v in km/s!
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                            :: i,j,ncra
TYPE (objstruc)                             :: obj
INTEGER(KIND=IB),DIMENSION((kmax+1)*(NPL-1)):: idxcra
INTEGER(KIND=IB),DIMENSION(kmax+1)          :: idxc,idxr,idxa
REAL(KIND=RP),DIMENSION(kmax+1)             :: c,r,a
REAL(KIND=RP)                               :: cl,ch,al,ah

idxcra= 0
idxc  = 0
idxr  = 0
CALL GETIDXCRA(obj,idxcra,ncra)
idxc = idxcra(1:ncra:(NPL-1))
idxr = idxcra(2:ncra:(NPL-1))
idxa = idxcra(3:ncra:(NPL-1))

c = 0
r = 0
c(1:obj%k+1) = obj%par(idxc(1:obj%k+1))
r(1:obj%k+1) = obj%par(idxr(1:obj%k+1))
a(1:obj%k+1) = obj%par(idxa(1:obj%k+1))

!!
!! ioutside must not be set to 0 here!! Would destroy the 
!! perturbation setting from earlier (in PROPOSAL)
!!
DO i=1,obj%k+1
   cl=(1.54_RP-0.907_RP*r(i)+0.3695_RP*r(i)**1.88_RP)*1.5004_RP
   ch=(1.6_RP-0.907_RP*r(i)+0.3695_RP*r(i)**2.01_RP)*1.5014_RP
   IF(c(i) > ch)THEN
      ioutside = 1
      IF(rank == src)WRITE(*,201) i,'ch',ioutside,c(i),ch
   ENDIF
   IF(c(i) < cl)THEN
      ioutside = 1
      IF(rank == src)WRITE(*,201) i,'cl',ioutside,c(i),cl
   ENDIF
ENDDO
!DO i=1,obj%k+1
!   al = (0.00002_RP*EXP(0.0019_RP*r(i)*1000._RP))*1.5294_RP
!   ah = (0.00332_RP*EXP(0.00265_RP*r(i)*1000._RP)+0.1_RP)*1.5294_RP
!   IF(a(i) > ah)THEN
!      ioutside = 1
!      IF(rank == src)WRITE(*,201) i,'ah:',ioutside,a(i),ah
!   ENDIF
!   IF(a(i) < al)THEN
!      ioutside = 1
!      IF(rank == src)WRITE(*,201) i,'al:',ioutside,a(i),al
!   ENDIF
!ENDDO

201 FORMAT(I3,a,I2,2F10.4)
END SUBROUTINE COUPLE_CR
!=======================================================================

SUBROUTINE SAVEREPLICA(obj,dat,filename,filename2)
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: i
TYPE (objstruc)  :: obj      ! Best object
TYPE (datastruc) :: dat      ! Best object
CHARACTER(len=64):: filename,filename2 ! replica file name

WRITE(6,*) 'Global best model:'
CALL PRINTPAR(obj)
WRITE(6,*) 'Global best logL = ',obj%logL
OPEN(UNIT=50,FILE=filename,FORM='formatted',STATUS='REPLACE', &
ACTION='WRITE',RECL=1024)
DO i = 1,NRF
   WRITE(50,208) dat%Drep(i,:)
ENDDO
CLOSE(50)

OPEN(UNIT=50,FILE=filename2,FORM='formatted',STATUS='REPLACE', &
ACTION='WRITE',RECL=1024)
DO i = 1,NRF
   WRITE(50,208) dat%Dobs(i,:)
ENDDO
CLOSE(50)

208 FORMAT(5000ES20.10)
RETURN
END SUBROUTINE SAVEREPLICA
!==============================================================================

SUBROUTINE MAKE_BETA(beta1,beta4,NTEMP1)
!!
!! Make cooling schedule (three geometrically spaced legs)
!!
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: itemp,NTEMP1
REAL(KIND=RP)                  :: beta1,beta2,beta3,beta4
REAL(KIND=RP),DIMENSION(NTEMP1):: logbeta

beta2 = beta1+0.2_RP*ABS(beta1-beta4)
beta3 = beta1+0.75_RP*ABS(beta1-beta4)
logbeta = 0._RP
beta = 0._RP
logbeta(1) = LOG(beta1)
DO itemp=2,NTEMP1/20
   logbeta(itemp) = logbeta(itemp-1)+(LOG(beta2)-LOG(beta1))/(REAL(NTEMP1/20,RP)-1)
ENDDO
DO itemp=NTEMP1/20+1,NTEMP1/3
   logbeta(itemp)=logbeta(itemp-1)+(LOG(beta3)-LOG(beta2))/(REAL(NTEMP1/3,RP)-REAL(NTEMP1/20,RP))
ENDDO
DO itemp=NTEMP1/3+1,NTEMP1
   logbeta(itemp)=logbeta(itemp-1)+(LOG(beta4)-LOG(beta3))/(REAL(NTEMP1,RP)-REAL(NTEMP1/3,RP))
ENDDO
beta = EXP(logbeta)

209 FORMAT(50ES16.8)
RETURN
END SUBROUTINE MAKE_BETA
!=======================================================================

SUBROUTINE PRINTPAR(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER :: i
TYPE(objstruc) :: obj

DO i=1,obj%k
   WRITE(6,201) obj%par((i-1)*NPL+1:i*NPL)
ENDDO
! write(6,*) 'write out parameter',obj%k,NPL,obj%k*NPL+1,obj%k*NPL+(NPL-1),obj%par(1:1)

WRITE(6,202) '            ',obj%par(obj%k*NPL+1:obj%k*NPL+(NPL-1))
WRITE(6,*) 'Partition:'
WRITE(6,203) obj%z
WRITE(6,*) 'Layers:'
WRITE(6,203) obj%h
WRITE(6,*) 'Water resistivity and depth:'
WRITE(6,203) obj%wres, obj%wdepth
IF(ICOV == 1 .OR. ICOV == 3)THEN
   WRITE(6,*) 'SD parameters:'
   WRITE(6,206) exp(obj%sdpar)
elseif(icov == 5) then
   WRITE(6,*) 'SD factors:'
   WRITE(6,206) obj%sdpar
ENDIF
IF(IAR == 1)THEN
   WRITE(6,*) 'AR parameters:'
   WRITE(6,203) obj%arpar
   WRITE(6,*) 'AR switch:'
   WRITE(6,205) obj%idxar
ENDIF
IF(ICF >= 1)THEN
   WRITE(6,*) 'CF parameters:'
   WRITE(6,203) obj%cf
ENDIF

201 FORMAT(5F12.4)
202 FORMAT(A12,8F12.4)
203 FORMAT(24F12.4)
205 FORMAT(24I8)
204 FORMAT(4F16.4)
206 FORMAT(14F14.10)
END SUBROUTINE PRINTPAR
!=======================================================================

SUBROUTINE PRINTVV()
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER :: i,j,l

WRITE(6,*) 'Rotation matrices:'
DO l=3,7
WRITE(6,*) 'k=',l
DO i=1,l*NPL+3
   WRITE(6,201) VV(i,:,l)
ENDDO
WRITE(6,*) 'sdevm:'
WRITE(6,201) sdevm(:,l)
WRITE(6,*) ''
ENDDO

201 FORMAT(100F8.3)
END SUBROUTINE PRINTVV

!=======================================================================

SUBROUTINE PRINTCOV0()
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER :: i,j,l

WRITE(6,*) 'Cov0 matrices:'
DO l=1,2
WRITE(6,*) 'k=',l
DO i=1,l*NPL+3
   WRITE(6,201) Cov0(i,1:l*NPL+3,l)
ENDDO
WRITE(6,*) ''
ENDDO

201 FORMAT(100F8.3)
END SUBROUTINE PRINTCOV0

!=======================================================================

FUNCTION RANDPERM(num)
!==============================================================================
USE data_type, ONLY : IB, RP
IMPLICIT NONE
INTEGER(KIND=IB), INTENT(IN) :: num
INTEGER(KIND=IB) :: numb, i, j, k
INTEGER(KIND=IB), DIMENSION(num) :: randperm
REAL(KIND=RP), DIMENSION(num) :: rand2
INTRINSIC RANDOM_NUMBER
      CALL RANDOM_NUMBER(rand2)
      DO i=1,num
        numb=1
        DO j=1,num
          IF (rand2(i) > rand2(j)) THEN
              numb=numb+1
          END IF
        END DO
      DO k=1,i-1
         IF (rand2(i) <= rand2(k) .AND. rand2(i) >= rand2(k)) THEN
           numb=numb+1
         END IF
      END DO
       randperm(i)=numb
      END DO
RETURN
END FUNCTION RANDPERM
!====================================================================

SUBROUTINE GASDEVJ(harvest)
!====================================================================
USE DATA_TYPE
!USE nrtype
USE nr
IMPLICIT NONE
REAL(RP), INTENT(OUT) :: harvest
REAL(RP) :: rsq,v1,v2
REAL(RP), SAVE :: g
LOGICAL, SAVE :: gaus_stored=.FALSE.
IF (gaus_stored) THEN
   harvest=g
   gaus_stored=.FALSE.
ELSE
   DO
      CALL RANDOM_NUMBER(v1)
      CALL RANDOM_NUMBER(v2)
      v1=2.0_RP*v1-1.0_RP
      v2=2.0_RP*v2-1.0_RP
      rsq=v1**2+v2**2
      IF (rsq > 0.0_RP .AND. rsq < 1.0_RP) EXIT
   END DO
   rsq=SQRT(-2.0_RP*LOG(rsq)/rsq)
   harvest=v1*rsq
   g=v2*rsq
   gaus_stored=.TRUE.
END IF
END SUBROUTINE GASDEVJ
!====================================================================

Function ASINC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,ii,asinc
ii    = cmplx(0.,1.)
asinc = -ii*LOG(ii*z+SQRT(1.-z**2))
RETURN
END FUNCTION
!=======================================================================

Function COSC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,cosc
REAL(KIND=RP)    :: x,y
x    = REAL(z)
y    = AIMAG(z)
cosc = CMPLX(COS(x)*COSH(y),-SIN(x)*SINH(y),RP)
RETURN
END FUNCTION
!=======================================================================

Function SINC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,sinc
REAL(KIND=RP)    :: x,y
x    = REAL(z)
y    = AIMAG(z)
sinc = CMPLX(SIN(x)*COSH(y),COS(x)*SINH(y),RP)
RETURN
END FUNCTION
!==============================================================================

FUNCTION CACOS(z)
!==============================================================================

USE DATA_TYPE
COMPLEX(KIND=RP) :: CACOS
COMPLEX(KIND=RP) :: z
REAL(KIND=RP) :: zrp1,zrm1,zi,zizi,a1,a2,a,b

!CACOS = -CMPLX(0._RP,1._RP,RP)*LOG(z+CMPLX(0._RP,1._RP,RP)*SQRT(1._RP-z*z))
!!
!! This version from IDL; much faster than above
!!
zrp1 = REAL(z,RP)+1._RP
zrm1 = zrp1-2._RP
zi = AIMAG(z)
zizi = zi*zi
a1 = 0.5_RP*SQRT(zrp1*zrp1 + zizi)
a2 = 0.5_RP*SQRT(zrm1*zrm1 + zizi)
a = a1+a2
b = a1- a2
IF(zi >= 0._RP)THEN
   CACOS = ACOS(b) - CMPLX(0._RP,1._RP,RP)*LOG(a + SQRT(a*a - 1))
ELSE
   CACOS = ACOS(b) + CMPLX(0._RP,1._RP,RP)*LOG(a + SQRT(a*a - 1))
ENDIF

RETURN
END FUNCTION CACOS
!==============================================================================

FUNCTION ASINH(x)
!==============================================================================

USE DATA_TYPE
REAL(KIND=RP) :: ASINH
REAL(KIND=RP) :: x

ASINH = LOG(x+SQRT(x**2._RP+1))

RETURN
END FUNCTION ASINH
!==============================================================================

FUNCTION CSIN(z)
!==============================================================================
!! Complex sine (Jan's version)

USE DATA_TYPE
COMPLEX(KIND=RP) :: CSIN
COMPLEX(KIND=RP) :: z

CSIN =  (EXP( CMPLX(0._RP,1._RP,RP)*z) -EXP(-CMPLX(0._RP,1._RP,RP)*z)) &
                            /CMPLX(0._RP,2._RP,RP)
RETURN
END FUNCTION CSIN
!==============================================================================

FUNCTION CTAN(z)
!==============================================================================
!! Complex TAN

USE DATA_TYPE
COMPLEX(KIND=RP) :: CTAN
COMPLEX(KIND=RP) :: z

CTAN =  -CMPLX(0._RP,1._RP,RP)*(EXP( CMPLX(0._RP,1._RP,RP)*z) -EXP(-CMPLX(0._RP,1._RP,RP)*z)) &
                            /(EXP( CMPLX(0._RP,1._RP,RP)*z)+EXP(-CMPLX(0._RP,1._RP,RP)*z))
RETURN
END FUNCTION CTAN
!=======================================================================

Function LOGPLUS(x,y)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
REAL(KIND=RP)    :: logplus,x,y

IF(x > y)THEN
   LOGPLUS = x+LOG(1._RP+EXP(y-x))
ELSE
   LOGPLUS = y+log(1._RP+EXP(x-y))
ENDIF

RETURN
END FUNCTION LOGPLUS
!=======================================================================

      Subroutine CHOLESKY(A,n,np,L)
!=======================================================================
!  Cholesky decomposition of symmetric, positive-definite matix A
!  into lower-triangular matrix L. Modified from Numerical Recipes.
!-----------------------------------------------------------------------
USE DATA_TYPE
IMPLICIT NONE
INTEGER      :: n,np,i,j,k
REAL(KIND=RP) :: A(np,np),A2(np,np),L(np,np),p(n),summ

A2 = A
DO i=1,n
   DO j=1,n
      IF (i > j) A2(i,j) = 0._RP
   ENDDO
ENDDO

DO i=1,n
   DO j=i,n
      summ = A2(i,j)
      DO k=i-1,1,-1
         summ = summ-A2(i,k)*A2(j,k)
      ENDDO
      IF (i == j) THEN
         IF (summ <= 0.) THEN
            WRITE(6,*) 'Cholesky Decomp Failed ',summ
            STOP
         ENDIF
         p(i) = SQRT(summ)
      ELSE
         A2(j,i) = summ/p(i)
      ENDIF
   ENDDO
ENDDO

DO i=1,n
   DO j=1,n
      IF (i > j) L(i,j) = A2(i,j)
      IF (i == j) L(i,i) = p(i)
   ENDDO
ENDDO

RETURN
END SUBROUTINE
!=======================================================================
! This is the end my friend...
! EOF
