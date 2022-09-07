subroutine pyk2_init(random_generator_seed)
  use mod_ranlux ,       only : rluxgo     ! for ranlux init
  use coll_common ,      only : rnd_seed !, rcx, rcxp, rcy, rcyp, rcp, rcs, &

  implicit none

  integer, intent(in)          :: random_generator_seed

  rnd_seed = random_generator_seed

  ! Initialize random number generator
  !if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  call rluxgo(3, rnd_seed, 0, 0)

end subroutine


subroutine make_ruth_dist(cgen, zatom, emr, hcut)
  use mod_funlux ,       only : funlxp
  use coll_k2,           only : zatom_curr, emr_curr, k2coll_ruth

  implicit none

  real(kind=8), intent(inout) :: cgen(200)
  real(kind=8), intent(in)    :: zatom
  real(kind=8), intent(in)    :: emr
  real(kind=8), intent(in)    :: hcut
  real(kind=8), parameter     :: tlcut = 0.0009982

  zatom_curr = zatom
  emr_curr = emr
  call funlxp(k2coll_ruth, cgen(1), tlcut, hcut)

end subroutine


real(kind=8) function pyk2_rand() 
  use mod_ranlux, only: coll_rand

  implicit none

  pyk2_rand = coll_rand()

end function pyk2_rand


subroutine pyk2_funlux(array,xran,len)
  use mod_funlux, only: funlux

  implicit none

  real(kind=8), intent(in)  :: array(200)
  integer, intent(in)       :: len
  real(kind=8), intent(inout) :: xran(len)

  call funlux(array,xran,len)

end subroutine


real(kind=8) function pyk2_rand_gauss(cut)
  use mod_ranlux, only: ran_gauss

  implicit none
  real(kind=8), intent(in)    :: cut

  pyk2_rand_gauss = ran_gauss(cut)

end function pyk2_rand_gauss




subroutine pyk2_cryinteract(ci_x,xp,y,yp,pc,length,s_P,x_P,ci_exenergy,ci_rho,ci_anuc,ci_zatom,ci_emr,&
                      ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt,ci_hcut,ci_csref0,ci_csref1,ci_csref4,&
                      ci_csref5,ci_bnref,ci_csect,ci_tilt,ci_rcurv,ci_alayer,ci_xmax,ci_ymax,ci_orient,&
                      ci_miscut,ci_bend,ci_cBend,ci_sBend,ci_cpTilt,ci_spTilt,ci_cnTilt,ci_snTilt,ci_iProc,&
                      dest,c_v1,c_v2,nam,zn,mom,enr,gammar,betar,bgr,mep,tmax,plen,const_dech,&
                      s,s_length,L_chan,s_K,x_K,s_M,x_M,r,A_F,B_F,C_F,x_F,s_F,&
                      tdefl,xp_rel,ymin,ymax)
    use coll_crystal
    use mod_ranlux
    use mod_funlux
    use mod_common_main
    use floatPrecision
    ! use coll_materials, only : zatom, exenergy, rho, anuc
    use mathlib_bouncer
    use physical_constants
    implicit none
  
    ! integer,          intent(in)    :: is  ! Material number
    real(kind=8), intent(inout) :: ci_x
    real(kind=8), intent(inout) :: xp
    real(kind=8), intent(inout) :: y
    real(kind=8), intent(inout) :: yp
    real(kind=8), intent(inout) :: pc
    real(kind=8), intent(in)    :: length
    real(kind=8), intent(in)    :: s_P
    real(kind=8), intent(in)    :: x_P
  
    real(kind=8), intent(in)    :: ci_exenergy
    real(kind=8), intent(in)    :: ci_rho
    real(kind=8), intent(in)    :: ci_anuc
    real(kind=8), intent(in)    :: ci_zatom
    real(kind=8), intent(in)    :: ci_emr
    
    real(kind=8), intent(in)    :: ci_dlri
    real(kind=8), intent(in)    :: ci_dlyi
    real(kind=8), intent(in)    :: ci_ai
    real(kind=8), intent(in)    :: ci_eUm
    real(kind=8), intent(in)    :: ci_collnt
    real(kind=8), intent(in)    :: ci_hcut
    real(kind=8), intent(in)    :: ci_csref0
    real(kind=8), intent(in)    :: ci_csref1
    real(kind=8), intent(in)    :: ci_csref4
    real(kind=8), intent(in)    :: ci_csref5
  
    real(kind=8), intent(in)    :: ci_bnref
    real(kind=8), intent(in)    :: ci_csect

    real(kind=8), intent(in)    :: ci_tilt
    real(kind=8), intent(in)    :: ci_rcurv
    real(kind=8), intent(in)    :: ci_alayer
    real(kind=8), intent(in)    :: ci_xmax
    real(kind=8), intent(in)    :: ci_ymax
    integer,      intent(in)    :: ci_orient
    real(kind=8), intent(in)    :: ci_miscut
    real(kind=8), intent(in)    :: ci_bend
    real(kind=8), intent(in)    :: ci_cBend
    real(kind=8), intent(in)    :: ci_sBend
    real(kind=8), intent(in)    :: ci_cpTilt
    real(kind=8), intent(in)    :: ci_spTilt
    real(kind=8), intent(in)    :: ci_cnTilt
    real(kind=8), intent(in)    :: ci_snTilt
    integer,      intent(inout) :: ci_iProc

    real(kind=8), intent(inout) :: dest
    real(kind=8), intent(in)    :: c_v1
    real(kind=8), intent(in)    :: c_v2
    integer,      intent(in)    :: nam
    integer,      intent(in)    :: zn
    real(kind=8), intent(in)    :: mom
    real(kind=8), intent(in)    :: enr
    real(kind=8), intent(in)    :: gammar
    real(kind=8), intent(in)    :: betar
    real(kind=8), intent(in)    :: bgr
    real(kind=8), intent(in)    :: mep
    real(kind=8), intent(in)    :: tmax
    real(kind=8), intent(in)    :: plen
    real(kind=8), intent(in)    :: const_dech
                      
    real(kind=8), intent(inout)    :: s
    real(kind=8), intent(inout)    :: s_length
    real(kind=8), intent(inout)    :: L_chan
    real(kind=8), intent(inout)    :: s_K
    real(kind=8), intent(inout)    :: x_K
    real(kind=8), intent(inout)    :: s_M
    real(kind=8), intent(inout)    :: x_M
    real(kind=8), intent(inout)    :: r
    real(kind=8), intent(inout)    :: A_F
    real(kind=8), intent(inout)    :: B_F
    real(kind=8), intent(inout)    :: C_F
    real(kind=8), intent(inout)    :: x_F
    real(kind=8), intent(inout)    :: s_F

    real(kind=8), intent(in)    :: tdefl
    real(kind=8), intent(in)    :: xp_rel
    real(kind=8), intent(in)    :: ymin
    real(kind=8), intent(in)    :: ymax

    ! integer nam,zn                        ! Switch on/off the nuclear interaction (NAM) and the MCS (ZN)
    ! real(kind=8) ymax,ymin            ! Crystal geometrical parameters
    ! real(kind=8) s_length             ! Element length along s
    ! real(kind=8) DESt                 ! Changed energy loss by ionization now calculated and not tabulated
    real(kind=8) x0,y0                ! Coordinates of the particle [m]
    ! real(kind=8) s                    ! Long coordinates of the particle [m]
    real(kind=8) a_eq,b_eq,c_eq,delta ! Second order equation param.
    real(kind=8) Ang_rms, Ang_avr     ! Volume reflection mean angle [rad]
    real(kind=8) Dechan               ! Probability for dechanneling
    real(kind=8) Lrefl,Srefl          ! Distance of the reflection point [m]
    real(kind=8) Vcapt                ! Volume capture probability
    real(kind=8) Chann                ! Channeling probability
    real(kind=8) N_atom               ! Probability for entering channeling near atomic planes
    real(kind=8) Dxp                  ! Variation in angle
    real(kind=8) xpcrit               ! Critical angle for curved crystal[rad]
    real(kind=8) xpcrit0              ! Critical angle for str. crystal [rad]
    real(kind=8) Rcrit                ! Critical curvature radius [m]
    real(kind=8) ratio                ! X=c_rcurv/Rcrit
    real(kind=8) TLdech1,TLdech2      ! Typical dechanneling length(2) [m]
    real(kind=8) tdech,Ldech,Sdech    ! Angle, lenght, and S coordinate of dechanneling point
    real(kind=8) Rlength,Red_S        ! Reduced length/s coordinate (in case of dechanneling)
    real(kind=8) am_len               ! Amorphous length
    real(kind=8) len_xs,len_ys        ! Amorphous length
    ! real(kind=8) xp_rel               ! Xp-c_miscut angle in mrad
    real(kind=8) alpha                ! Par for new chann prob
    real(kind=8) Pvr                  ! Prob for VR->AM transition
  
    ! Quantities for length and deflection calculation
    ! real(kind=8) const_dech,xpin,ypin,tchan,tdefl,tP,L_chan,mep
    ! real(kind=8) s_K,x_K,s_M,x_M,s_F,x_F,r,a
    ! real(kind=8) A_F,B_F,C_F,alpha_F,beta_F
    real(kind=8) alpha_F,beta_F,a,tP,xpin,ypin,tchan
  
    ! real(kind=8), parameter :: c_v1 =  1.7 ! Fitting coefficient
    ! real(kind=8), parameter :: c_v2 = -1.5 ! Fitting coefficient
  
    ! real(kind=8) enr
    ! real(kind=8) mom
    ! real(kind=8) betar
    ! real(kind=8) bgr
    ! real(kind=8) gammar
    ! real(kind=8) tmax
    ! real(kind=8) plen
  
    !   nam = 1 ! Switch on/off the nuclear interaction (NAM) and the MCS (ZN)
    !   zn  = 1
  
    !   ! dE/dX and dechanneling length calculation
    !   mom    = pc*c1e3                ! [GeV]
    !   enr    = sqrt(mom**2 + pmap**2) ! [MeV]
    !   gammar = enr/pmap
    !   betar  = mom/enr
    !   bgr = betar*gammar
    !   mep    = pmae/pmap  ! Electron/proton
  
    !   tmax = (two*pmae*bgr**2)/(one + two*gammar*mep + mep**2)  ! [MeV]
    !   plen = sqrt((ci_rho*ci_zatom)/ci_anuc)*28.816e-6 ! [MeV]
  
    !   const_dech = ((256.0/(nine*pi**2)) * &
    !     (one/(log_mb(((two*pmae)*gammar)/(ci_exenergy*c1e3)) - one))) * ((aTF*dP)/(crade*pmae)) ! [m/MeV]
    !   const_dech = const_dech*c1e3 ! [m/GeV]
  
    !   s        = zero
    !   s_length = ci_rcurv*sin_mb(length/ci_rcurv)
    !   L_chan   = length
  
    ! ! MISCUT second step: fundamental coordinates (crystal edges and plane curvature radius)
    !   s_K = ci_rcurv*sin_mb(length/ci_rcurv)
    !   x_K = ci_rcurv*(1-cos_mb(length/ci_rcurv))
    !   s_M = (ci_rcurv-ci_xmax)*sin_mb(length/ci_rcurv)
    !   x_M = ci_xmax + (ci_rcurv-ci_xmax)*(1-cos_mb(length/ci_rcurv))
    !   r   = sqrt(s_P**2 + (ci_x-x_P)**2)
  
    !   ! MISCUT third step: F coordinates (exit point) on crystal exit face
    !   A_F = (tan_mb(length/ci_rcurv))**2 + one
    !   B_F = ((-two)*(tan_mb(length/ci_rcurv))**2)*ci_rcurv + (two*tan_mb(length/ci_rcurv))*s_P - two*x_P
    !   C_F = ((tan_mb(length/ci_rcurv))**2)*(ci_rcurv**2) - ((two*tan_mb(length/ci_rcurv))*s_P)*ci_rcurv + s_P**2 + x_P**2 - r**2
  
    !   x_F = (-B_F-sqrt(B_F**2-four*(A_F*C_F)))/(two*A_F)
    !   s_F = (-tan_mb(length/ci_rcurv))*(x_F-ci_rcurv)
  
    ! if(x_F >= x_K .and. x_F <= x_M .and. s_F >= s_M .and. s_F <= s_K) then
    !   ! No additional steps required for miscut
    ! else if (ci_miscut == 0 .and. abs(x_F-x_K) <= c1m13 .and. abs(s_F-s_K) <= c1m13) then
    !   ! no miscut, entrance from below: exit point is K (lower edge)
    !   x_F = x_K
    !   s_F = s_K
    ! else if (ci_miscut == 0 .and. abs(x_F-x_M) <= c1m13 .and. abs(s_F-s_M) <= c1m13) then
    !   ! no miscut, entrance from above: exit point is M (upper edge)
    !   x_F = x_M
    !   s_F = s_M
    ! else
  
    !   ! MISCUT Third step (bis): F coordinates (exit point)  on bent side
    !   if(ci_miscut < 0) then
    !     ! Intersect with bottom side
    !     alpha_F = (ci_rcurv-x_P)/x_P
    !     beta_F = -(r**2-s_P**2-x_P**2)/(two*s_P)
    !     A_F = alpha_F**2 + one
    !     B_F = two*(alpha_F*beta_F) - two*ci_rcurv
    !     C_F = beta_F**2
    !   else
    !     ! Intersect with top side
    !     alpha_F = (ci_rcurv-x_P)/s_P
    !     beta_F = -(r**2+ci_xmax*(ci_xmax-(two*ci_rcurv))-s_P**2-x_P**2)/(two*s_P)
    !     A_F = alpha_F**2 + one
    !     B_F = two*(alpha_F*beta_F) - two*ci_rcurv
    !     C_F = beta_F**2 - ci_xmax*(ci_xmax-two*ci_rcurv)
    !   endif
      
    !   x_F = (-B_F-sqrt(B_F**2-four*(A_F*C_F)))/(two*A_F)
    !   s_F = alpha_F*x_F + beta_F
      
    ! endif
  
    ! ! MISCUT fourth step: deflection and length calculation
    ! a = sqrt(s_F**2+(ci_x-x_F)**2)
    ! tP = acos_mb((2*(r**2)-a**2)/(2*(r**2)))
    ! tdefl = asin_mb((s_f-s_P)/r)
    ! L_chan = r*tP
  
    ! xp_rel = xp - ci_miscut
  
    ! ymin = -ci_ymax/two
    ! ymax =  ci_ymax/two
  
    ! ! FIRST CASE: p don't interact with crystal
    ! if(y < ymin .or. y > ymax .or. ci_x > ci_xmax) then
    !   ci_x  = ci_x + xp*s_length
    !   y     = y + yp*s_length
    !   ci_iProc = proc_out
    !   return
  
    ! ! SECOND CASE: p hits the amorphous layer
    ! else if(ci_x < ci_alayer .or. y-ymin < ci_alayer .or. ymax-y < ci_alayer) then
    !   x0    = ci_x
    !   y0    = y
    !   a_eq  = one + xp**2
    !   b_eq  = (two*ci_x)*xp - (two*xp)*ci_rcurv
    !   c_eq  = ci_x**2 - (two*ci_x)*ci_rcurv
    !   delta = b_eq**2 - (four*a_eq)*c_eq
    !   s     = (-b_eq+sqrt(delta))/(two*a_eq)
    !   if(s >= s_length) then
    !     s = s_length
    !   end if
    !   ci_x   =  xp*s + x0
    !   len_xs = sqrt((ci_x-x0)**2 + s**2)
    !   if(yp >= zero .and. y + yp*s <= ymax) then
    !     len_ys = yp*len_xs
    !   else if(yp < zero .and. y + yp*s >= ymin) then
    !     len_ys = yp*len_xs
    !   else
    !     s      = (ymax-y)/yp
    !     len_ys = sqrt((ymax-y)**2 + s**2)
    !     ci_x   = x0 + xp*s
    !     len_xs = sqrt((ci_x-x0)**2 + s**2)
    !   end if
    !   am_len = sqrt(len_xs**2 + len_ys**2)
    !   s     = s/two
    !   ci_x  = x0 + xp*s
    !   y     = y0 + yp*s
    !   ci_iProc = proc_AM

    !   PRINT *, "BEFORE crycalcionloss1"
    !   PRINT *, pc
    !   PRINT *, am_len
    !   PRINT *, dest
    !   PRINT *, betar
    !   PRINT *, bgr
    !   PRINT *, gammar
    !   PRINT *, tmax
    !   PRINT *, plen
    !   PRINT *, ci_exenergy
    !   PRINT *, ci_zatom
    !   PRINT *, ci_rho
    !   PRINT *, ci_anuc
    !   PRINT *, ci_dlri
    !   PRINT *, ci_dlyi
    !   PRINT *, ci_ai
    !   PRINT *, ci_eUm
    !   PRINT *, ci_collnt
    !   PRINT *
    !   call pyk2_crycalcionloss(pc,am_len,dest,betar,bgr,gammar,tmax,plen,&
    !                         ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !   PRINT *, "AFTER crycalcionloss1"
    !   PRINT *, pc
    !   PRINT *, am_len
    !   PRINT *, dest
    !   PRINT *, betar
    !   PRINT *, bgr
    !   PRINT *, gammar
    !   PRINT *, tmax
    !   PRINT *, plen
    !   PRINT *, ci_exenergy
    !   PRINT *, ci_zatom
    !   PRINT *, ci_rho
    !   PRINT *, ci_anuc
    !   PRINT *, ci_dlri
    !   PRINT *, ci_dlyi
    !   PRINT *, ci_ai
    !   PRINT *, ci_eUm
    !   PRINT *, ci_collnt
    !   PRINT *
      
    !   PRINT *, "BEFORE crymoveam1"
    !   PRINT *, xp
    !   PRINT *, yp
    !   PRINT *, pc
    !   PRINT *, ci_iProc
    !   PRINT *, nam
    !   PRINT *, s_length
    !   PRINT *, dest
    !   PRINT *, ci_dlyi
    !   PRINT *, ci_dlri
    !   PRINT *, xp
    !   PRINT *, yp
    !   PRINT *, pc
    !   PRINT *, ci_anuc
    !   PRINT *, ci_zatom
    !   PRINT *, ci_emr
    !   PRINT *, ci_hcut
    !   PRINT *, ci_bnref
    !   PRINT *, ci_csref0      
    !   PRINT *, ci_csref1
    !   PRINT *, ci_csref4
    !   PRINT *, ci_csref5
    !   PRINT *, ci_collnt
    !   PRINT *, ci_iProc
    !   PRINT *
    !   call pyk2_crymoveam(nam,am_len,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csref0,&
    !                   ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
    !   PRINT *, "AFTER crymoveam1"
    !   PRINT *, xp
    !   PRINT *, yp
    !   PRINT *, pc
    !   PRINT *, ci_iProc
    !   PRINT *, nam
    !   PRINT *, s_length
    !   PRINT *, dest
    !   PRINT *, ci_dlyi
    !   PRINT *, ci_dlri
    !   PRINT *, xp
    !   PRINT *, yp
    !   PRINT *, pc
    !   PRINT *, ci_anuc
    !   PRINT *, ci_zatom
    !   PRINT *, ci_emr
    !   PRINT *, ci_hcut
    !   PRINT *, ci_bnref
    !   PRINT *, ci_csref0      
    !   PRINT *, ci_csref1
    !   PRINT *, ci_csref4
    !   PRINT *, ci_csref5
    !   PRINT *, ci_collnt
    !   PRINT *, ci_iProc
    !   PRINT *
    !   ci_x = ci_x + xp*(s_length-s)
    !   y = y + yp*(s_length-s)
      
    !   return
  
    ! else if(ci_x > ci_xmax-ci_alayer .and. ci_x < ci_xmax) then
    !   ci_iProc = proc_AM

    !   PRINT *, "BEFORE crycalcionloss2"
    !   PRINT *, pc
    !   PRINT *, am_len
    !   PRINT *, dest
    !   PRINT *, betar
    !   PRINT *, bgr
    !   PRINT *, gammar
    !   PRINT *, tmax
    !   PRINT *, plen
    !   PRINT *, ci_exenergy
    !   PRINT *, ci_zatom
    !   PRINT *, ci_rho
    !   PRINT *, ci_anuc
    !   PRINT *, ci_dlri
    !   PRINT *, ci_dlyi
    !   PRINT *, ci_ai
    !   PRINT *, ci_eUm
    !   PRINT *, ci_collnt
    !   PRINT *
    !   call pyk2_crycalcionloss(pc,s_length,dest,betar,bgr,gammar,tmax,plen,&
    !                       ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !   PRINT *, "AFTER crycalcionloss2"
    !   PRINT *, pc
    !   PRINT *, am_len
    !   PRINT *, dest
    !   PRINT *, betar
    !   PRINT *, bgr
    !   PRINT *, gammar
    !   PRINT *, tmax
    !   PRINT *, plen
    !   PRINT *, ci_exenergy
    !   PRINT *, ci_zatom
    !   PRINT *, ci_rho
    !   PRINT *, ci_anuc
    !   PRINT *, ci_dlri
    !   PRINT *, ci_dlyi
    !   PRINT *, ci_ai
    !   PRINT *, ci_eUm
    !   PRINT *, ci_collnt
    !   PRINT *

    !   PRINT *, "BEFORE crymoveam2"
    !   PRINT *, xp
    !   PRINT *, yp
    !   PRINT *, pc
    !   PRINT *, ci_iProc
    !   PRINT *, nam
    !   PRINT *, s_length
    !   PRINT *, dest
    !   PRINT *, ci_dlyi
    !   PRINT *, ci_dlri
    !   PRINT *, xp
    !   PRINT *, yp
    !   PRINT *, pc
    !   PRINT *, ci_anuc
    !   PRINT *, ci_zatom
    !   PRINT *, ci_emr
    !   PRINT *, ci_hcut
    !   PRINT *, ci_bnref
    !   PRINT *, ci_csref0      
    !   PRINT *, ci_csref1
    !   PRINT *, ci_csref4
    !   PRINT *, ci_csref5
    !   PRINT *, ci_collnt
    !   PRINT *, ci_iProc
    !   PRINT *
    !   call pyk2_crymoveam(nam,s_length,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csref0,&
    !                   ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
    !   PRINT *, "AFTER crymoveam2"
    !   PRINT *, xp
    !   PRINT *, yp
    !   PRINT *, pc
    !   PRINT *, ci_iProc
    !   PRINT *, nam
    !   PRINT *, s_length
    !   PRINT *, dest
    !   PRINT *, ci_dlyi
    !   PRINT *, ci_dlri
    !   PRINT *, xp
    !   PRINT *, yp
    !   PRINT *, pc
    !   PRINT *, ci_anuc
    !   PRINT *, ci_zatom
    !   PRINT *, ci_emr
    !   PRINT *, ci_hcut
    !   PRINT *, ci_bnref
    !   PRINT *, ci_csref0      
    !   PRINT *, ci_csref1
    !   PRINT *, ci_csref4
    !   PRINT *, ci_csref5
    !   PRINT *, ci_collnt
    !   PRINT *, ci_iProc
    !   PRINT *
      
    !   return
  
    ! end if


    ! ! THIRD CASE: the p interacts with the crystal.
    ! ! Define typical angles/probabilities for orientation 110
    ! xpcrit0 = sqrt((c2m9*ci_eUm)/pc)    ! Critical angle (rad) for straight crystals
    ! Rcrit   = (pc/(c2m6*ci_eUm))*ci_ai ! Critical curvature radius [m]
  
    ! ! If R>Rcritical=>no channeling is possible (ratio<1)
    ! ratio  = ci_rcurv/Rcrit
    ! xpcrit = (xpcrit0*(ci_rcurv-Rcrit))/ci_rcurv ! Critical angle for curved crystal
  
    ! if(ratio <= one) then ! no possibile channeling
    !   Ang_rms = ((c_v1*0.42)*xpcrit0)*sin_mb(1.4*ratio) ! RMS scattering
    !   Ang_avr = ((c_v2*xpcrit0)*c5m2)*ratio                         ! Average angle reflection
    !   Vcapt   = zero                                                ! Probability of VC
  
    ! else if(ratio <= three) then ! Strongly bent crystal
    !   Ang_rms = ((c_v1*0.42)*xpcrit0)*sin_mb(0.4713*ratio + 0.85) ! RMS scattering
    !   Ang_avr = (c_v2*xpcrit0)*(0.1972*ratio - 0.1472)                  ! Average angle reflection
    !   Vcapt   = 7.0e-4*(ratio - 0.7)/pc**c2m1                           ! Correction by sasha drozdin/armen
    !   ! K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)
  
    ! else ! Rcry >> Rcrit
    !   Ang_rms = (c_v1*xpcrit0)*(one/ratio)                ! RMS scattering
    !   Ang_avr = (c_v2*xpcrit0)*(one - 1.6667/ratio) ! Average angle for VR
    !   Vcapt   = 7.0e-4*(ratio - 0.7)/pc**c2m1 ! Probability for VC correction by sasha drozdin/armen
    !   ! K=0.0007 is taken based on simulations using CATCH.f (V.Biryukov)
  
    ! end if
  
    ! if(ci_orient == 2) then
    !   Ang_avr = Ang_avr*0.93
    !   Ang_rms = Ang_rms*1.05
    !   xpcrit  = xpcrit*0.98
    ! end if
  
    ! if(abs(xp_rel) < xpcrit) then
    !   alpha  = xp_rel/xpcrit
    !   Chann  = sqrt(0.9*(one - alpha**2))*sqrt(one-(one/ratio)) ! Saturation at 95%
    !   N_atom = c1m1
  
    !   ! if they can channel: 2 options
    !   if(coll_rand() <= chann) then ! option 1:channeling
  
    !     TLdech1 = (const_dech*pc)*(one-one/ratio)**2 ! Updated calculate typical dech. length(m)
    !     if(coll_rand() <= n_atom) then
    !       TLdech1 = ((const_dech/c2e2)*pc)*(one-one/ratio)**2  ! Updated dechanneling length (m)
    !     end if
  
    !     Dechan = -log_mb(coll_rand()) ! Probability of dechanneling
    !     Ldech  = TLdech1*Dechan   ! Actual dechan. length
  
    !     ! careful: the dechanneling lentgh is along the trajectory
    !     ! of the particle -not along the longitudinal coordinate...
    !     if(ldech < l_chan) then
    !       ci_iProc = proc_DC
    !       Dxp   = Ldech/r ! Change angle from channeling [mrad]
    !       Sdech = Ldech*cos_mb(ci_miscut + half*Dxp)
    !       ci_x  = ci_x  + Ldech*(sin_mb(half*Dxp+ci_miscut)) ! Trajectory at channeling exit
    !       xp    = xp + Dxp + (two*(coll_rand()-half))*xpcrit
    !       y     = y  + yp * Sdech
  
    !       call pyk2_crycalcionloss(pc,ldech,dest,betar,bgr,gammar,tmax,plen,&
    !                             ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !       pc = pc - half*dest*Ldech ! Energy loss to ionization while in CH [GeV]
    !       ci_x  = ci_x  + (half*(s_length-Sdech))*xp
    !       y  = y  + (half*(s_length-Sdech))*yp
  
    !       call pyk2_crycalcionloss(pc,s_length-sdech,dest,betar,bgr,gammar,tmax,plen,&
    !                             ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !       call pyk2_crymoveam(nam,s_length-sdech,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,&
    !                       ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
    !       ci_x = ci_x + (half*(s_length-Sdech))*xp
    !       y = y + (half*(s_length-Sdech))*yp
    !     else
    !       ci_iProc = proc_CH
    !       xpin  = XP
    !       ypin  = YP
  
    !       ! check if a nuclear interaction happen while in CH
    !       call pyk2_crymovech(nam,l_chan,ci_x,xp,yp,pc,ci_rcurv,rcrit,ci_rho,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csect,&
    !                       ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_eUm,ci_collnt,ci_iProc)
    !       if(ci_iProc /= proc_CH) then
    !         ! if an nuclear interaction happened, move until the middle with initial xp,yp then
    !         ! propagate until the "crystal exit" with the new xp,yp accordingly with the rest
    !         ! of the code in "thin lens approx"
    !         ci_x = ci_x + (half*L_chan)*xpin
    !         y = y + (half*L_chan)*ypin
    !         ci_x = ci_x + (half*L_chan)*XP
    !         y = y + (half*L_chan)*YP
  
    !         call pyk2_crycalcionloss(pc,length,dest,betar,bgr,gammar,tmax,plen,&
    !                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !         pc = pc - dest*length ! energy loss to ionization [GeV]
    !       else
    !         Dxp = tdefl + (half*ran_gauss(zero))*xpcrit ! Change angle[rad]
          
    !         xp  = Dxp
    !         ci_x = ci_x + L_chan*(sin_mb(half*Dxp)) ! Trajectory at channeling exit
    !         y   = y + s_length * yp
  
    !         call pyk2_crycalcionloss(pc,length,dest,betar,bgr,gammar,tmax,plen,&
    !                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !         pc = pc - (half*dest)*length ! energy loss to ionization [GeV]
    !       end if
    !     end if
  
    !   else ! Option 2: VR
  
    !     ! good for channeling but don't channel (1-2)
    !     ci_iProc = proc_VR
  
    !     xp = xp + (0.45*(xp_rel/xpcrit + one))*Ang_avr
    !     ci_x  = ci_x  + (half*s_length)*xp
    !     y  = y  + (half*s_length)*yp
  
    !     call pyk2_crycalcionloss(pc,s_length,dest,betar,bgr,gammar,tmax,plen,&
    !                           ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !     call pyk2_crymoveam(nam,s_length,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csref0,&
    !                     ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
  
    !     ci_x = ci_x + (half*s_length)*xp
    !     y = y + (half*s_length)*yp
  
    !   end if
  
    ! else ! case 3-2: no good for channeling. check if the  can VR
  
    !   Lrefl = xp_rel*r ! Distance of refl. point [m]
    !   Srefl = sin_mb(xp_rel/two + ci_miscut)*Lrefl
  
    !   if(Lrefl > zero .and. Lrefl < L_chan) then ! VR point inside
  
    !     ! 2 options: volume capture and volume reflection
  
    !     if(coll_rand() > Vcapt .or. ZN == zero) then ! Option 1: VR
  
    !       ci_iProc = proc_VR
    !       ci_x  = ci_x + xp*Srefl
    !       y     = y + yp*Srefl
    !       Dxp   = Ang_avr
    !       xp    = xp + Dxp + Ang_rms*ran_gauss(zero)
    !       ci_x  = ci_x  + (half*xp)*(s_length - Srefl)
    !       y     = y  + (half*yp)*(s_length - Srefl)
  
    !       call pyk2_crycalcionloss(pc,s_length-srefl,dest,betar,bgr,gammar,tmax,plen,&
    !                             ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !       call pyk2_crymoveam(nam,s_length-srefl,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,&
    !                       ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
    !       ci_x = ci_x + (half*xp)*(s_length - Srefl)
    !       y = y + (half*yp)*(s_length - Srefl)
  
    !     else ! Option 2: VC
  
    !       ci_x = ci_x + xp*Srefl
    !       y = y + yp*Srefl
  
    !       TLdech2 = (const_dech/c1e1)*pc*(one-one/ratio)**2          ! Updated typical dechanneling length(m)
    !       Ldech   = TLdech2*(sqrt(c1m2 - log_mb(coll_rand())) - c1m1)**2 ! Updated DC length
    !       tdech   = Ldech/ci_rcurv
    !       Sdech   = Ldech*cos_mb(xp + half*tdech)
  
    !       if(Ldech < Length-Lrefl) then
  
    !         ci_iProc = proc_DC
    !         Dxp   = Ldech/ci_rcurv + (half*ran_gauss(zero))*xpcrit
    !         ci_x  = ci_x + Ldech*(sin_mb(half*Dxp+xp)) ! Trajectory at channeling exit
    !         y     = y + Sdech*yp
    !         xp    =  Dxp
    !         Red_S = (s_length - Srefl) - Sdech
    !         ci_x  = ci_x + (half*xp)*Red_S
    !         y     = y + (half*yp)*Red_S
  
    !         call pyk2_crycalcionloss(pc,srefl,dest,betar,bgr,gammar,tmax,plen,&
    !                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !         pc = pc - dest*Srefl ! "added" energy loss before capture
  
    !         call pyk2_crycalcionloss(pc,sdech,dest,betar,bgr,gammar,tmax,plen,&
    !                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !         pc = pc - (half*dest)*Sdech ! "added" energy loss while captured
  
    !         call pyk2_crycalcionloss(pc,red_s,dest,betar,bgr,gammar,tmax,plen,&
    !                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !         call pyk2_crymoveam(nam,red_s,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csref0,&
    !                         ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
    !         ci_x = ci_x + (half*xp)*Red_S
    !         y = y + (half*yp)*Red_S
  
    !       else
  
    !         ci_iProc   = proc_VC
    !         Rlength = Length - Lrefl
    !         tchan   = Rlength/ci_rcurv
    !         Red_S   = Rlength*cos_mb(xp + half*tchan)
  
    !         call pyk2_crycalcionloss(pc,lrefl,dest,betar,bgr,gammar,tmax,plen,&
    !                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !         pc   = pc - dest*Lrefl ! "added" energy loss before capture
    !         xpin = xp
    !         ypin = yp
  
    !         ! Check if a nuclear interaction happen while in ch
    !         call pyk2_crymovech(nam,rlength,ci_x,xp,yp,pc,ci_rcurv,rcrit,ci_rho,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csect,&
    !                         ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_eUm,ci_collnt,ci_iProc)
                            
    !         if(ci_iProc /= proc_VC) then
    !           ! if an nuclear interaction happened, move until the middle with initial xp,yp then propagate until
    !           ! the "crystal exit" with the new xp,yp aciordingly with the rest of the code in "thin lens approx"
    !           ci_x = ci_x + (half*Rlength)*xpin
    !           y = y + (half*Rlength)*ypin
    !           ci_x = ci_x + (half*Rlength)*XP
    !           y = y + (half*Rlength)*YP
  
    !           call pyk2_crycalcionloss(pc,rlength,dest,betar,bgr,gammar,tmax,plen,&
    !                                 ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !           pc = pc - dest*Rlength
    !         else
    !           Dxp = (Length-Lrefl)/ci_rcurv
    !           ci_x = ci_x + sin_mb(half*Dxp+xp)*Rlength ! Trajectory at channeling exit
    !           y   = y + red_S*yp
    !           xp  = tdefl + (half*ran_gauss(zero))*xpcrit ! [mrad]
  
    !           call pyk2_crycalcionloss(pc,rlength,dest,betar,bgr,gammar,tmax,plen,&
    !                                 ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !           pc = pc - (half*dest)*Rlength  ! "added" energy loss once captured
    !         end if
    !       end if
    !     end if
  
    !   else
  
    !     ! Case 3-3: move in amorphous substance (big input angles)
    !     ! Modified for transition vram daniele
    !     if(xp_rel > tdefl-ci_miscut + two*xpcrit .or. xp_rel < -xpcrit) then
    !       ci_iProc = proc_AM
    !       ci_x  = ci_x + (half*s_length)*xp
    !       y     = y + (half*s_length)*yp
    !       if(zn > zero) then
    !         call pyk2_crycalcionloss(pc,s_length,dest,betar,bgr,gammar,tmax,plen,&
    !                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !         call pyk2_crymoveam(nam,s_length,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,ci_csref0,&
    !                         ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
    !       end if
    !       ci_x = ci_x + (half*s_length)*xp
    !       y = y + (half*s_length)*yp
    !     else
    !       Pvr = (xp_rel-(tdefl-ci_miscut))/(two*xpcrit)
    !       if(coll_rand() > Pvr) then
    !         ci_iProc = proc_TRVR
    !         ci_x  = ci_x + xp*Srefl
    !         y     = y + yp*Srefl
  
    !         Dxp = (((-three*Ang_rms)*xp_rel)/(two*xpcrit) + Ang_avr) + ((three*Ang_rms)*(tdefl-ci_miscut))/(two*xpcrit)
    !         xp  = xp + Dxp
    !         ci_x = ci_x + (half*xp)*(s_length-Srefl)
    !         y   = y + (half*yp)*(s_length-Srefl)
  
    !         call pyk2_crycalcionloss(pc,s_length-srefl,dest,betar,bgr,gammar,tmax,plen,&
    !                               ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !         call pyk2_crymoveam(nam,s_length-srefl,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,&
    !                         ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
    !         ci_x = ci_x + (half*xp)*(s_length - Srefl)
    !         y = y + (half*yp)*(s_length - Srefl)
    !       else
    !         ci_iProc = proc_TRAM
    !         ci_x = ci_x + xp*Srefl
    !         y = y + yp*Srefl
    !         Dxp = ((((-one*(13.6/pc))*sqrt(s_length/ci_dlri))*c1m3)*xp_rel)/(two*xpcrit) + &
    !           (((13.6/pc)*sqrt(s_length/ci_dlri))*c1m3)*(one+(tdefl-ci_miscut)/(two*xpcrit))
    !         xp = xp+Dxp
    !         ci_x  = ci_x + (half*xp)*(s_length-Srefl)
    !         y  = y + (half*yp)*(s_length-Srefl)
  
    !         call pyk2_crycalcionloss(pc,s_length-srefl,dest,betar,bgr,gammar,tmax,plen,&
    !                             ci_exenergy,ci_zatom,ci_rho,ci_anuc,ci_dlri,ci_dlyi,ci_ai,ci_eUm,ci_collnt)
    !         call pyk2_crymoveam(nam,s_length-srefl,dest,ci_dlyi,ci_dlri,xp,yp,pc,ci_anuc,ci_zatom,ci_emr,ci_hcut,ci_bnref,&
    !                         ci_csref0,ci_csref1,ci_csref4,ci_csref5,ci_collnt,ci_iProc)
    !         ci_x = ci_x + (half*xp)*(s_length - Srefl)
    !         y = y + (half*yp)*(s_length - Srefl)
    !       end if
    !     end if
    !   end if
    ! end if
                    

end subroutine pyk2_cryinteract


subroutine pyk2_crycalcionloss(pc,dz,EnLo,cc_betar,cc_bgr,cc_gammar,cc_tmax,cc_plen,cc_exenergy,cc_zatom,&
                           cc_rho,cc_anuc,cc_dlri,cc_dlyi,cc_ai,cc_eUm,cc_collnt)

  use mod_ranlux
  use mod_funlux
  use floatPrecision
  ! use coll_materials, only : zatom, exenergy, rho, anuc
  use mathlib_bouncer
  use physical_constants

  ! integer,          intent(in)  :: is
  real(kind=8), intent(in)  :: pc
  real(kind=8), intent(in)  :: dz
  real(kind=8), intent(inout) :: EnLo
  real(kind=8), intent(in)  :: cc_betar
  real(kind=8), intent(in)  :: cc_bgr
  real(kind=8), intent(in)  :: cc_gammar
  real(kind=8), intent(in)  :: cc_tmax
  real(kind=8), intent(in)  :: cc_plen

  real(kind=8), intent(in)  :: cc_exenergy
  real(kind=8), intent(in)  :: cc_zatom
  real(kind=8), intent(in)  :: cc_rho
  real(kind=8), intent(in)  :: cc_anuc
  
  real(kind=8), intent(in)  :: cc_dlri
  real(kind=8), intent(in)  :: cc_dlyi
  real(kind=8), intent(in)  :: cc_ai
  real(kind=8), intent(in)  :: cc_eUm
  real(kind=8), intent(in)  :: cc_collnt

  
  real(kind=fPrec) thl,tt,cs_tail,prob_tail
  real(kind=fPrec), parameter :: k = 0.307075_fPrec ! Constant in front bethe-bloch [mev g^-1 cm^2]

  thl       = (((((four*k)*cc_zatom)*dz)*c1e2)*cc_rho)/(cc_anuc*cc_betar**2) ! [MeV]
  EnLo      = ((k*cc_zatom)/(cc_anuc*cc_betar**2)) * ( &
    half*log_mb(((((two*pmae)*cc_bgr)*cc_bgr)*cc_tmax)/(c1e6*cc_exenergy**2)) - &
    cc_betar**2 - log_mb(cc_plen/(cc_exenergy*c1e3)) - log_mb(cc_bgr) + half    &
  )
  EnLo      = ((EnLo*cc_rho)*c1m1)*dz ! [GeV]
  Tt        = (EnLo*c1e3)+thl          ! [MeV]

  cs_tail   = ((k*cc_zatom)/(cc_anuc*cc_betar**2)) * ((half*((one/Tt)-(one/cc_tmax))) - &
    (log_mb(cc_tmax/Tt)*(cc_betar**2)/(two*cc_tmax)) + ((cc_tmax-Tt)/((four*(cc_gammar**2))*(pmap**2))))
  prob_tail = ((cs_tail*cc_rho)*dz)*c1e2

  if(coll_rand() < prob_tail) then
    EnLo = ((k*cc_zatom)/(cc_anuc*cc_betar**2)) * ( &
      half*log_mb((two*pmae*cc_bgr*cc_bgr*cc_tmax)/(c1e6*cc_exenergy**2)) -      &
      cc_betar**2 - log_mb(cc_plen/(cc_exenergy*c1e3)) - log_mb(cc_bgr) + half + &
      cc_tMax**2/(eight*(cc_gammar**2)*(pmap**2)) &
    )
    EnLo = (EnLo*cc_rho)*c1m1 ! [GeV/m]
  else
    EnLo = EnLo/dz ! [GeV/m]
  end if

end subroutine pyk2_crycalcionloss



subroutine pyk2_crymoveam(nam,dz,dei,dly,dlr,xp,yp,pc,cm_anuc,cm_zatom,cm_emr,cm_hcut,cm_bnref,&
                    cm_csref0,cm_csref1,cm_csref4,cm_csref5,cm_collnt,cm_iProc)

  use coll_crystal, only: cry_moveAM

  ! integer,          intent(in)    :: is
  integer,      intent(in)    :: nam
  real(kind=8), intent(in)    :: dz
  real(kind=8), intent(in)    :: dei
  real(kind=8), intent(in)    :: dly
  real(kind=8), intent(in)    :: dlr
  real(kind=8), intent(inout) :: xp
  real(kind=8), intent(inout) :: yp
  real(kind=8), intent(inout) :: pc

  real(kind=8), intent(in)    :: cm_csref0
  real(kind=8), intent(in)    :: cm_csref1
  real(kind=8), intent(in)    :: cm_csref4
  real(kind=8), intent(in)    :: cm_csref5

  real(kind=8), intent(in)    :: cm_anuc
  real(kind=8), intent(in)    :: cm_zatom
  real(kind=8), intent(in)    :: cm_emr
  real(kind=8), intent(in)    :: cm_hcut
  real(kind=8), intent(in)    :: cm_bnref
  real(kind=8), intent(in)    :: cm_collnt
  integer,      intent(inout) :: cm_iProc

  call cry_moveAM(nam,dz,dei,dly,dlr,xp,yp,pc,cm_anuc,cm_zatom,cm_emr,cm_hcut,cm_bnref,&
                    cm_csref0,cm_csref1,cm_csref4,cm_csref5,cm_collnt,cm_iProc)

end subroutine pyk2_crymoveam


subroutine pyk2_crymovech(nam,dz,ch_x,xp,yp,pc,r,rc,ch_rho,ch_anuc,ch_zatom,ch_emr,ch_hcut,ch_bnref,ch_csect,&
                      ch_csref0,ch_csref1,ch_csref4,ch_csref5,ch_eUm,ch_collnt,ch_iProc)
  use coll_crystal, only: cry_moveCH

  integer,      intent(in)    :: nam
  real(kind=8), intent(in)    :: dz
  real(kind=8), intent(inout) :: ch_x
  real(kind=8), intent(inout) :: xp
  real(kind=8), intent(inout) :: yp
  real(kind=8), intent(inout) :: pc
  real(kind=8), intent(in)    :: r
  real(kind=8), intent(in)    :: rc

  real(kind=8), intent(in)    :: ch_rho
  real(kind=8), intent(in)    :: ch_anuc
  real(kind=8), intent(in)    :: ch_zatom
  real(kind=8), intent(in)    :: ch_emr
  real(kind=8), intent(in)    :: ch_hcut
  real(kind=8), intent(in)    :: ch_bnref

  real(kind=8), intent(in)    :: ch_csref0
  real(kind=8), intent(in)    :: ch_csref1
  real(kind=8), intent(in)    :: ch_csref4
  real(kind=8), intent(in)    :: ch_csref5

  real(kind=8), intent(in)    :: ch_csect
  real(kind=8), intent(in)    :: ch_eUm
  real(kind=8), intent(in)    :: ch_collnt
  integer,      intent(inout) :: ch_iProc

  call cry_moveCH(nam,dz,ch_x,xp,yp,pc,r,rc,ch_rho,ch_anuc,ch_zatom,ch_emr,ch_hcut,ch_bnref,ch_csect,&
                      ch_csref0,ch_csref1,ch_csref4,ch_csref5,ch_eUm,ch_collnt,ch_iProc)

end subroutine pyk2_crymovech
