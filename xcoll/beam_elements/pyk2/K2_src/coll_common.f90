

! ================================================================================================ !
!  Collimation Common Variables
! ================================================================================================ !
module coll_common

  use floatPrecision
  use numerical_constants, only : zero

  implicit none

  ! Logical Flags
  logical, save :: coll_debug         = .true.
  logical, save :: dowrite_impact     = .false.
  logical, save :: dowrite_dist       = .false.
  logical, save :: dowrite_secondary  = .false.
  logical, save :: dowrite_amplitude  = .false.
  logical, save :: dowrite_tracks     = .false.
  logical, save :: dowrite_efficiency = .false.
  logical, save :: dowrite_crycoord   = .false.
  logical, save :: coll_hasCrystal    = .false.

  ! Various Variables
  integer, save :: rnd_seed   = 0

  ! Collimation Particle Arrays
  real(kind=fPrec), allocatable, save :: rcx(:)
  real(kind=fPrec), allocatable, save :: rcxp(:)
  real(kind=fPrec), allocatable, save :: rcy(:)
  real(kind=fPrec), allocatable, save :: rcyp(:)
  real(kind=fPrec), allocatable, save :: rcp(:)
  real(kind=fPrec), allocatable, save :: rcs(:)

  ! Process index for interaction with crystals
  integer, allocatable, save :: cry_proc(:)
  integer, allocatable, save :: cry_proc_prev(:)
  integer, allocatable, save :: cry_proc_tmp(:)

  ! Pencil Beam
  integer,          save :: ipencil       = 0
  integer,          save :: pencil_distr  = 0
  real(kind=fPrec), save :: pencil_offset = zero
  real(kind=fPrec), save :: pencil_rmsx   = zero
  real(kind=fPrec), save :: pencil_rmsy   = zero
  real(kind=fPrec), save :: xp_pencil0    = zero
  real(kind=fPrec), save :: yp_pencil0    = zero
  real(kind=fPrec), allocatable, save :: x_pencil(:)
  real(kind=fPrec), allocatable, save :: y_pencil(:)
  real(kind=fPrec), allocatable, save :: xp_pencil(:)
  real(kind=fPrec), allocatable, save :: yp_pencil(:)
  real(kind=fPrec), allocatable, save :: pencil_dx(:)

  ! Other Arrays
  integer,          allocatable, save :: cn_impact(:)
  integer,          allocatable, save :: cn_absorbed(:)
  real(kind=fPrec), allocatable, save :: caverage(:)
  real(kind=fPrec), allocatable, save :: csigma(:)
  real(kind=fPrec), allocatable, save :: gap_rms_error(:)
  real(kind=fPrec), allocatable, save :: csum(:)
  real(kind=fPrec), allocatable, save :: csqsum(:)

  ! Output File Names
  character(len=12), parameter :: coll_survivalFile   = "survival.dat"
  character(len=12), parameter :: coll_gapsFile       = "collgaps.dat"
  character(len=10), parameter :: coll_impactFile     = "impact.dat"
  character(len=11), parameter :: coll_tracksFile     = "tracks2.dat"
  character(len=17), parameter :: coll_positionsFile  = "CollPositions.dat"
  character(len=20), parameter :: coll_pencilFile     = "pencilbeam_distr.dat"
  character(len=16), parameter :: coll_ellipseFile    = "coll_ellipse.dat"
  character(len=15), parameter :: coll_allImpactFile  = "all_impacts.dat"
  character(len=19), parameter :: coll_allAbsorbFile  = "all_absorptions.dat"
  character(len=16), parameter :: coll_scatterFile    = "Coll_Scatter.dat"
  character(len=16), parameter :: coll_fstImpactFile  = "FirstImpacts.dat"
  character(len=17), parameter :: coll_flukImpFile    = "FLUKA_impacts.dat"
  character(len=21), parameter :: coll_flukImpAllFile = "FLUKA_impacts_all.dat"
  character(len=17), parameter :: coll_sigmaSetFile   = "sigmasettings.out"
  character(len=16), parameter :: coll_settingsFile   = "collsettings.dat"
  character(len=16), parameter :: coll_jawProfileFile = "jaw_profiles.dat"
  character(len=13), parameter :: coll_ampFile        = "amplitude.dat"
  character(len=17), parameter :: coll_orbitCheckFile = "orbitchecking.dat"
  character(len=16), parameter :: coll_summaryFile    = "coll_summary.dat"
  character(len=14), parameter :: coll_efficFile      = "efficiency.dat"
  character(len=19), parameter :: coll_efficDPFile    = "efficiency_dpop.dat"
  character(len=17), parameter :: coll_effic2DFile    = "efficiency_2d.dat"
  character(len=16), parameter :: coll_cryEntFile     = "cry_entrance.dat"
  character(len=12), parameter :: coll_cryExitFile    = "cry_exit.dat"
  character(len=19), parameter :: coll_cryInterFile   = "cry_interaction.dat"

  ! Output File Units
  integer, save :: outlun              = -1
  integer, save :: coll_survivalUnit   = -1
  integer, save :: coll_gapsUnit       = -1
  integer, save :: coll_impactUnit     = -1
  integer, save :: coll_tracksUnit     = -1
  integer, save :: coll_positionsUnit  = -1
  integer, save :: coll_pencilUnit     = -1
  integer, save :: coll_ellipseUnit    = -1
  integer, save :: coll_allImpactUnit  = -1
  integer, save :: coll_allAbsorbUnit  = -1
  integer, save :: coll_scatterUnit    = -1
  integer, save :: coll_fstImpactUnit  = -1
  integer, save :: coll_flukImpUnit    = -1
  integer, save :: coll_flukImpAllUnit = -1
  integer, save :: coll_sigmaSetUnit   = -1
  integer, save :: coll_settingsUnit   = -1
  integer, save :: coll_jawProfileUnit = -1
  integer, save :: coll_ampUnit        = -1
  integer, save :: coll_orbitCheckUnit = -1
  integer, save :: coll_summaryUnit    = -1
  integer, save :: coll_efficUnit      = -1
  integer, save :: coll_efficDPUnit    = -1
  integer, save :: coll_effic2DUnit    = -1
  integer, save :: coll_cryEntUnit     = -1
  integer, save :: coll_cryExitUnit    = -1
  integer, save :: coll_cryInterUnit   = -1



contains

subroutine coll_expandArrays(npart_new)

  use mod_alloc
  use numerical_constants

  integer, intent(in) :: npart_new

  call alloc(rcx,  npart_new, zero, "rcx")
  call alloc(rcxp, npart_new, zero, "rcxp")
  call alloc(rcy,  npart_new, zero, "rcy")
  call alloc(rcyp, npart_new, zero, "rcyp")
  call alloc(rcp,  npart_new, zero, "rcp")
  call alloc(rcs,  npart_new, zero, "rcs")

  call alloc(cry_proc, npart_new, -1, "cry_proc")
  call alloc(cry_proc_prev, npart_new, -1, "cry_proc_prev")
  call alloc(cry_proc_tmp, npart_new, -1, "cry_proc_tmp")

end subroutine coll_expandArrays

subroutine coll_expandNColl(nColl)

  use mod_alloc
  use numerical_constants

  integer, intent(in) :: nColl

  call alloc(cn_impact,     nColl, 0,    "cn_impact")
  call alloc(cn_absorbed,   nColl, 0,    "cn_absorbed")
  call alloc(caverage,      nColl, zero, "caverage")
  call alloc(csigma,        nColl, zero, "csigma")
  call alloc(gap_rms_error, nColl, zero, "gap_rms_error")
  call alloc(csum,          nColl, zero, "csum")
  call alloc(csqsum,        nColl, zero, "csqsum")
  call alloc(x_pencil,      nColl, zero, "x_pencil")
  call alloc(y_pencil,      nColl, zero, "y_pencil")
  call alloc(xp_pencil,     nColl, zero, "xp_pencil")
  call alloc(yp_pencil,     nColl, zero, "yp_pencil")
  call alloc(pencil_dx,     nColl, zero, "pencil_dx")

end subroutine coll_expandNColl

end module coll_common

! ================================================================================================ !
!  Collimator Material Variables
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Cross section inputs and material property database
! ================================================================================================ !
module coll_materials

  use floatPrecision
  use numerical_constants, only : zero, one, c1e10

  implicit none

  !integer, parameter :: nmat  = 16 ! Total number of materials
  !integer, parameter :: nrmat = 14 ! Number of real materials

  ! pp cross-sections and parameters for energy dependence
  real(kind=fPrec), parameter :: pptref = 0.04_fPrec
  real(kind=fPrec), parameter :: freeco = 1.618_fPrec

  ! Collimator Material Arrays
  !character(4),     public, save :: colmats(nmat)   ! Material Names
  !real(kind=fPrec), public, save :: exenergy(nmat)  ! Mean excitation energy [GeV]
  !real(kind=fPrec), public, save :: anuc(nmat)      ! Atomic mass
  !real(kind=fPrec), public, save :: zatom(nmat)     ! Atomic Z
  !real(kind=fPrec), public, save :: rho(nmat)       ! Density
  !real(kind=fPrec), public, save :: emr(nmat)       ! Nuclear radius
  !real(kind=fPrec), public, save :: hcut(nmat)      ! T cut (upper)
  !real(kind=fPrec), public, save :: radl(nmat)      ! Remaining length
  !real(kind=fPrec), public, save :: bnref(nmat)     ! Nuclear elastic slope from Schiz et al., PRD 21(3010)1980
  real(kind=fPrec), public, save :: csect(0:5) ! Cross section
  real(kind=fPrec), public, save :: xintl     ! Interaction length
  real(kind=fPrec), public, save :: bn        ! Nuclear elastic related
  real(kind=fPrec), public, save :: freep     ! Number of nucleons involved in single scattering
  real(kind=fPrec), public, save :: cgen(200)  ! Used by FUNLUX / Rutherford routine

  ! All cross-sections are in barns. Nuclear values from RPP at 20 GeV
  ! Coulomb is integerated above t=tLcut[Gev2] (+-1% out Gauss mcs)

  ! In Cs and CsRef,1st index: Cross-sections for processes
  ! 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
  ! 4:Single Diffractive pp or pn, 5:Coulomb for t above mcs

  ! Claudia 2013: updated cross section values. Unit: Barn. New 2013:
  !real(kind=fPrec), public, parameter :: csref(0:5,nmat) = reshape([ &
  !  [0.271_fPrec, 0.192_fPrec, zero, zero, zero, 0.0035e-2_fPrec], & ! BE
  !  [0.643_fPrec, 0.418_fPrec, zero, zero, zero, 0.0340e-2_fPrec], & ! AL
  !  [1.253_fPrec, 0.769_fPrec, zero, zero, zero, 0.1530e-2_fPrec], & ! CU
  !  [2.765_fPrec, 1.591_fPrec, zero, zero, zero, 0.7680e-2_fPrec], & ! W
  !  [3.016_fPrec, 1.724_fPrec, zero, zero, zero, 0.9070e-2_fPrec], & ! PB
  !  [0.337_fPrec, 0.232_fPrec, zero, zero, zero, 0.0076e-2_fPrec], & ! C
  !  [0.337_fPrec, 0.232_fPrec, zero, zero, zero, 0.0076e-2_fPrec], & ! C2
  !  [0.664_fPrec, 0.430_fPrec, zero, zero, zero, 0.0390e-2_fPrec], & ! Si
  !  [1.388_fPrec, 0.844_fPrec, zero, zero, zero, 0.1860e-2_fPrec], & ! Ge
  !  [0.362_fPrec, 0.247_fPrec, zero, zero, zero, 0.0094e-2_fPrec], & ! MoGR
  !  [0.572_fPrec, 0.370_fPrec, zero, zero, zero, 0.0279e-2_fPrec], & ! CuCD
  !  [1.713_fPrec, 1.023_fPrec, zero, zero, zero, 0.2650e-2_fPrec], & ! Mo
  !  [1.246_fPrec, 0.765_fPrec, zero, zero, zero, 0.1390e-2_fPrec], & ! Glid
  !  [2.548_fPrec, 1.473_fPrec, zero, zero, zero, 0.5740e-2_fPrec], & ! Iner
  !  [       zero,        zero, zero, zero, zero,            zero], & ! VA
  !  [       zero,        zero, zero, zero, zero,            zero]  & ! BL
  !], shape=[6,nmat])

  ! Cprob to choose an interaction in iChoix
  real(kind=fPrec), public, save :: cprob(0:5) = [zero, zero, zero, zero, zero, one]

  ! Electron density and plasma energy
  real(kind=fPrec), public, save :: edens = zero
  real(kind=fPrec), public, save :: pleng = zero

contains


! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-16
!  Updated: 2019-09-16
!  Get collimator material number from name (case sensitive)
! ================================================================================================ !
!integer function collmat_getCollMatID(matName)
!
!  character(len=*), intent(in) :: matName
!  integer i, matID
!
!  matID = -1
!  do i=1,nmat
!    if(colmats(i) == matName) then
!      matID = i
!      exit
!    end if
!  end do
!  collmat_getCollMatID = matID
!
!end function collmat_getCollMatID

end module coll_materials
