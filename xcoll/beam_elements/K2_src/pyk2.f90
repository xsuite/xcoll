subroutine pyk2_init(n_alloc, random_generator_seed)
  use floatPrecision
  use numerical_constants
  ! use crcoall    NODIG ??
  use parpro ,           only : npart
  use mod_alloc ,        only : alloc      ! to allocate partID etc
  use mod_common ,       only : iexact, napx, unit208, aa0
  use mod_common_main ,  only : partID, parentID, pairID, naa
  use mod_ranlux ,       only : rluxgo     ! for ranlux init

  use coll_common ,      only : rnd_seed, rcx, rcxp, rcy, rcyp, rcp, rcs, &
                                coll_expandArrays
  !use coll_materials ! for collmat_init
  use coll_k2         ! for scattering

  implicit none

  integer, intent(in)          :: n_alloc
  integer, intent(in)          :: random_generator_seed

  ! Set default values for collimator materials
 ! call collmat_init

  rnd_seed = random_generator_seed

  ! Initialize random number generator
  !if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  call rluxgo(3, rnd_seed, 0, 0)

  call coll_expandArrays(n_alloc)
  call alloc(naa, n_alloc, aa0, "naa")
  call alloc(partID, n_alloc, 0, "partID")
  call alloc(parentID, n_alloc, 0, "parentID")
  call alloc(pairID, 2, n_alloc, 0, "pairID")

end subroutine






subroutine pyk2_start_run(num_particles, cgen, hcut, zatom, emr, &
                  x_part, xp_part, y_part, yp_part, s_part, p_part)

  use floatPrecision
  use numerical_constants

  use parpro ,           only : npart
  use mod_common ,       only : napx, unit208, aa0
  use mod_common_main ,  only : partID, parentID, pairID, naa
  use mod_funlux,        only : funlxp
  use coll_k2,           only : zatom_curr, emr_curr, k2coll_ruth
  use coll_common ,      only : rnd_seed, rcx, rcxp, rcy, rcyp, rcp, rcs
  
  implicit none

  integer, intent(in)             :: num_particles
  real(kind=fPrec), intent(inout) :: cgen(200)
  real(kind=fPrec), intent(in)    :: zatom   !
  real(kind=fPrec), intent(in)    :: emr     !
  real(kind=fPrec), intent(in)    :: hcut    !
  real(kind=8),    intent(inout)  :: x_part(num_particles)
  real(kind=8),    intent(inout)  :: xp_part(num_particles)
  real(kind=8),    intent(inout)  :: y_part(num_particles)
  real(kind=8),    intent(inout)  :: yp_part(num_particles)
  real(kind=8),    intent(inout)  :: s_part(num_particles)
  real(kind=8),    intent(inout)  :: p_part(num_particles)

  integer j
  real(kind=fPrec), parameter :: tlcut = 0.0009982

  npart=num_particles
  unit208=109 

  do j=1,npart
    naa(j) = aa0
    partID(j)   = j
    parentID(j) = j
    pairID(1,j) = (j+1)/2    ! The pairID of particle j
    pairID(2,j) = 2-mod(j,2) ! Either particle 1 or 2 of the pair
  end do
  
  napx=npart  ! this decreases after absorptions!

  ! Prepare for Rutherford differential distribution
  !mcurr = mat ! HACK> mcurr is global, and coll_zatom too which is used inside k2coll_ruth
  zatom_curr = zatom
  emr_curr = emr

  call funlxp(k2coll_ruth, cgen(1), tlcut, hcut)
  
  do j=1,npart
    rcx(j) = x_particles(j)
    rcxp(j) = xp_particles(j)
    rcy(j) = y_particles(j)
    rcyp(j) = yp_particles(j)
    rcs(j) = s_particles(j)
    rcp(j) = p_particles(j)
  end do
  
end subroutine






subroutine pyk2_per_particle(j, cgen, thisp0, nhit, nabs, fracab, mirror, cRot, sRot, cRRot, sRRot, &
            nnuc0, ien0, nnuc1, ien1, &
            coll_exenergy, coll_anuc, coll_zatom, coll_emr, coll_rho, coll_hcut, coll_bnref, & 
            coll_csref0, coll_csref1, coll_csref4, coll_csref5, coll_radl, coll_dlri, coll_dlyi,coll_eUm, coll_ai, &
            coll_collnt, coll_cprob, coll_xintl, coll_bn, coll_ecmsq, coll_xln15s, coll_bpp, is_crystal, c_length, &
            c_aperture, c_offset, c_tilt, x_in, xp_in, y_in, yp_in, p_in, s_in, lhit, part_abs_local, &
            impact, indiv, lint, onesided, nhit_stage, j_slices, nabs_type, linside)


    use, intrinsic :: iso_fortran_env, only : int16
    use parpro
    use crcoall
    !use coll_db
    use coll_common
    use coll_crystal, only : cry_doCrystal
    !use coll_materials
    use mod_common, only : iexact, napx, unit208
    use mod_common_main, only : partID, naa
    use mathlib_bouncer
    use mod_ranlux
    use mod_funlux
    use coll_k2

    implicit none

    integer(kind=4) , intent(in)    :: j
    real(kind=fPrec), intent(inout) :: cgen(200)
    real(kind=fPrec), intent(in)    :: thisp0           ! Reference momentum in GeV
    integer(kind=4),  intent(inout) :: nhit
    integer(kind=4),  intent(inout) :: nabs
    real(kind=fPrec), intent(inout) :: fracab
    real(kind=fPrec), intent(inout) :: mirror
    real(kind=fPrec), intent(in)    :: cRot
    real(kind=fPrec), intent(in)    :: sRot
    real(kind=fPrec), intent(in)    :: cRRot
    real(kind=fPrec), intent(in)    :: sRRot
    integer(kind=4) , intent(inout) :: nnuc0
    real(kind=fPrec), intent(inout) :: ien0
    integer(kind=4) , intent(inout) :: nnuc1
    real(kind=fPrec), intent(inout) :: ien1

    real(kind=fPrec), intent(inout) :: coll_exenergy!
    real(kind=fPrec), intent(in)    :: coll_anuc    ! 
    real(kind=fPrec), intent(in)    :: coll_zatom   !
    real(kind=fPrec), intent(in)    :: coll_emr     !
    real(kind=fPrec), intent(in)    :: coll_rho     ! 
    real(kind=fPrec), intent(in)    :: coll_hcut    ! 
    real(kind=fPrec), intent(in)    :: coll_bnref   !
    real(kind=fPrec), intent(in)    :: coll_csref0  ! 
    real(kind=fPrec), intent(in)    :: coll_csref1  ! 
    real(kind=fPrec), intent(in)    :: coll_csref4  ! 
    real(kind=fPrec), intent(in)    :: coll_csref5  !
    real(kind=fPrec), intent(in)    :: coll_radl    !

    real(kind=fPrec), intent(in)    :: coll_dlri
    real(kind=fPrec), intent(in)    :: coll_dlyi
    real(kind=fPrec), intent(in)    :: coll_eUm
    real(kind=fPrec), intent(in)    :: coll_ai
    real(kind=fPrec), intent(in)    :: coll_collnt
    real(kind=fPrec), intent(in)    :: coll_cprob(0:5)
    real(kind=fPrec), intent(in)    :: coll_xintl
    real(kind=fPrec), intent(inout) :: coll_bn
    real(kind=fPrec), intent(in)    :: coll_ecmsq
    real(kind=fPrec), intent(in)    :: coll_xln15s
    real(kind=fPrec), intent(in)    :: coll_bpp
    logical,          intent(in)    :: is_crystal

    real(kind=fPrec), intent(in)    :: c_length     ! Collimator length in m
    real(kind=fPrec), intent(in)    :: c_aperture   ! Collimator aperture in m
    real(kind=fPrec), intent(in)    :: c_offset     ! Collimator offset in m
    real(kind=fPrec), intent(inout) :: c_tilt(2)    ! Collimator tilt in radians

    real(kind=fPrec), intent(inout) :: x_in(npart)
    real(kind=fPrec), intent(inout) :: xp_in(npart)
    real(kind=fPrec), intent(inout) :: y_in(npart)
    real(kind=fPrec), intent(inout) :: yp_in(npart)
    real(kind=fPrec), intent(inout) :: p_in(npart)
    real(kind=fPrec), intent(inout) :: s_in(npart)

    logical,          intent(in)    :: onesided

    integer,          intent(inout) :: lhit(npart)
    integer,          intent(inout) :: part_abs_local(npart)
    integer,          intent(inout) :: nabs_type(npart)
    integer,          intent(inout) :: nhit_stage(npart)
    real(kind=fPrec), intent(inout) :: indiv(npart)
    real(kind=fPrec), intent(inout) :: lint(npart)
    real(kind=fPrec), intent(inout) :: impact(npart)
    logical,          intent(inout) :: linside(npart)

    logical isImp
    integer j_slices
    real(kind=fPrec) keeps,drift_length,tiltangle
    real(kind=fPrec) x00,z00,p,sp,s,s_impact
    real(kind=fPrec) x_flk,xp_flk,y_flk,yp_flk,s_flk,zpj
    real(kind=fPrec) x_Dump,xpDump,y_Dump,ypDump,s_Dump
    real(kind=fPrec) xIn,xpIn,yIn,ypIn,xOut,xpOut,yOut,ypOut,sImp,sOut
    real(kind=fPrec) x_in0,xp_in0



    length = c_length

    if(part_abs_local == 0) then
      ! Only do scattering process for particles not already absorbed

        impact(j) = -one
        lint(j)   = -one
        indiv(j)  = -one

        x      = x_in(j)
        xp     = xp_in(j)
        xp_in0 = xp_in(j)
        z      = y_in(j)
        zp     = yp_in(j)
        p      = p_in(j)
        sp     = zero
        dpop   = (p - p0)/p0
        x_flk  = zero
        y_flk  = zero
        xp_flk = zero
        yp_flk = zero

        ! Transform particle coordinates to get into collimator coordinate  system
        ! First do rotation into collimator frame
        x  =  x_in(j)*cRot + sRot(j)*y_in
        z  =  y_in(j)*cRot - sRot(j)*x_in
        xp = xp_in(j)*cRot + sRot(j)*yp_in
        zp = yp_in(j)*cRot - sRot(j)*xp_in

        ! For one-sided collimators consider only positive X. For negative X jump to the next particle
        if(.not. onesided .or. x >= zero) then

            ! Log input energy + nucleons as per the FLUKA coupling
            nnuc0   = nnuc0 + naa(j)
            ien0    = ien0 + rcp(j) * c1e3


            ! Now mirror at the horizontal axis for negative X offset
            if(x < zero) then
              mirror    = -one
              tiltangle = -one*c_tilt(2)
            else
              mirror    = one
              tiltangle = c_tilt(1)
            end if
            x  = mirror*x
            xp = mirror*xp

            ! Shift with opening and offset
            x = (x - c_aperture/two) - mirror*c_offset

            ! Include collimator tilt
            if(tiltangle > zero) then
              xp = xp - tiltangle
            end if
            if(tiltangle < zero) then
              x  = x + sin_mb(tiltangle) * c_length
              xp = xp - tiltangle
            end if

            ! CRY Only: x_in0 has to be assigned after the change of reference frame
            x_in0 = x

            ! After finishing the coordinate transformation, or the coordinate manipulations in case of pencil beams,
            ! save the initial coordinates of the impacting particles
            xIn  = x
            xpIn = xp
            yIn  = z
            ypIn = zp

            ! particle passing above the jaw are discarded => take new event
            ! entering by the face, shorten the length (zlm) and keep track of
            ! entrance longitudinal coordinate (keeps) for histograms

            ! The definition is that the collimator jaw is at x>=0.

            ! 1) Check whether particle hits the collimator
            isImp = .false.
            s     = zero
            keeps = zero
            zlm   = -one*c_length

            if(is_crystal) then ! This is a crystal collimator

              call cry_doCrystal(j,x,xp,z,zp,s,p,x_in0,xp_in0,zlm,sImp,isImp,nhit,nabs,lhit,&
                part_abs_local,impact,indiv,c_length,coll_exenergy,coll_anuc,coll_zatom,coll_emr,coll_rho,&
                coll_hcut,coll_bnref,coll_csref0,coll_csref1,coll_csref4,coll_csref5,coll_dlri,coll_dlyi,&
                coll_eUm,coll_ai,coll_collnt,coll_bn)

              if(nabs /= 0) then
                part_abs_local(j)  = 1
                lint(j)                = zlm
              end if

              sImp  = (s - c_length) + sImp
              sOut  = s
              xOut  = x
              xpOut = xp
              yOut  = z
              ypOut = zp

            else ! "Normal" collimator

              if(x >= zero) then
                ! Particle hits collimator and we assume interaction length ZLM equal
                ! to collimator length (what if it would leave collimator after
                ! small length due to angle???)
                zlm       = c_length
                impact(j) = x
                indiv(j)  = xp
              else if(xp <= zero) then
                ! Particle does not hit collimator. Interaction length ZLM is zero.
                zlm = zero
              else
                ! Calculate s-coordinate of interaction point
                s = (-one*x)/xp
                if(s <= zero) then
                  write(lerr,"(a)") "COLLK2> ERROR S <= zero. This should not happen!"
                  !call prror
                end if
                if(s < c_length) then
                  zlm       = c_length - s
                  impact(j) = zero
                  indiv(j)  = xp
                else
                  zlm = zero
                end if
              end if

              ! First do the drift part
              ! DRIFT PART
              drift_length = c_length - zlm
              if(drift_length > zero) then
                if(iexact) then
                  zpj = sqrt(one-xp**2-zp**2)
                  x   = x  + drift_length*(xp/zpj)
                  z   = z  + drift_length*(zp/zpj)
                  sp  = sp + drift_length
                else
                  x  = x  + xp* drift_length
                  z  = z  + zp * drift_length
                  sp = sp + drift_length
                end if
              end if

              ! Now do the scattering part
              if(zlm > zero) then
                if(.not.linside(j)) then
                  ! first time particle hits collimator: entering jaw
                  linside(j) = .true.
                  if(dowrite_impact) then
                    if(tiltangle > zero) then
                      x_Dump = (x + c_aperture/two + tiltangle*sp)*mirror + c_offset
                    else
                      x_Dump = (x + c_aperture/two + tiltangle*(sp - c_length))*mirror + c_offset
                    end if
                    xpDump = (xp + tiltangle)*mirror
                    y_Dump = z
                    ypDump = zp
                    s_Dump = sp+real(j_slices-1,fPrec)*c_length
                  end if
                end if

                s_impact = sp
                nhit = nhit + 1
                call k2coll_jaw(s,nabs,partID(j),coll_exenergy,coll_anuc,coll_zatom,coll_rho,coll_radl,&
                                  coll_cprob,coll_xintl,coll_bn,cgen,coll_ecmsq,coll_xln15s,coll_bpp)
                nabs_type(j) = nabs
                lhit(j)  = 1

                isImp = .true.
                sImp  = s_impact+(real(j_slices,fPrec)-one)*c_length
                sOut  = (s+sp)+(real(j_slices,fPrec)-one)*c_length
                xOut  = x
                xpOut = xp
                yOut  = z
                ypOut = zp

                ! Writeout should be done for both inelastic and single diffractive. doing all transformations
                ! in x_flk and making the set to 99.99 mm conditional for nabs=1
                if(dowrite_impact .or. nabs == 1 .or. nabs == 4) then
                  ! Transform back to lab system for writeout.
                  ! keep x,y,xp,yp unchanged for continued tracking, store lab system variables in x_flk etc

                  x_flk  = xInt
                  xp_flk = xpInt

                  if(tiltangle > zero) then
                    x_flk  = x_flk  + tiltangle*(sInt+sp)
                    xp_flk = xp_flk + tiltangle
                  else if(tiltangle < zero) then
                    xp_flk = xp_flk + tiltangle
                    x_flk  = x_flk  - sin_mb(tiltangle) * (c_length-(sInt+sp))
                  end if

                  x_flk  = (x_flk + c_aperture/two) + mirror*c_offset
                  x_flk  = mirror*x_flk
                  xp_flk = mirror*xp_flk
                  y_flk  = (  yInt*cRRot -  x_flk*sRRot)*c1e3
                  yp_flk = ( ypInt*cRRot - xp_flk*sRRot)*c1e3
                  x_flk  = ( x_flk*cRRot +   yInt*sRRot)*c1e3
                  xp_flk = (xp_flk*cRRot +  ypInt*sRRot)*c1e3
                  s_flk  = (sInt+sp)+(real(j_slices,fPrec)-one)*c_length

                  ! Finally, the actual coordinate change to 99 mm
                  if(nabs == 1) then
                    fracab  = fracab + 1
                    x       = 99.99e-3_fPrec
                    z       = 99.99e-3_fPrec
                    lint(j) = zlm
                    part_abs_local(j)  = 1
                  end if
                end if
              end if ! Collimator jaw interaction

              if(nabs /= 1 .and. zlm > zero) then
                ! Do the rest drift, if particle left collimator early
                drift_length = (c_length-(s+sp))
                if(drift_length > c1m15) then
                  linside(j) = .false.
                  if(dowrite_impact) then
                    if(tiltangle > zero) then
                      x_Dump = (x + c_aperture/two + tiltangle*(s+sp))*mirror + c_offset
                    else
                      x_Dump = (x + c_aperture/two + tiltangle*(s+sp-c_length))*mirror + c_offset
                    end if
                    xpDump = (xp+tiltangle)*mirror
                    y_Dump = z
                    ypDump = zp
                    s_Dump = s+sp+real(j_slices-1,fPrec)*c_length

                  end if
                  if(iexact) then
                    zpj = sqrt(one-xp**2-zp**2)
                    x   = x  + drift_length*(xp/zpj)
                    z   = z  + drift_length*(zp/zpj)
                    sp  = sp + drift_length
                  else
                    x  = x  + xp * drift_length
                    z  = z  + zp * drift_length
                    sp = sp + drift_length
                  end if
                end if
                lint(j) = zlm - drift_length
              end if

            end if ! Collimator isCrystal

            ! Transform back to particle coordinates with opening and offset
            if(x < 99.0e-3_fPrec) then
              ! Include collimator tilt
              if(tiltangle > zero) then
                x  = x  + tiltangle*c_length
                xp = xp + tiltangle
              else if(tiltangle < zero) then
                x  = x  + tiltangle*c_length
                xp = xp + tiltangle
                x  = x  - sin_mb(tiltangle) * c_length
              end if

              ! Transform back to particle coordinates with opening and offset
              z00 = z
              x00 = x + mirror*c_offset
              x   = (x + c_aperture/two) + mirror*c_offset

              ! Now mirror at the horizontal axis for negative X offset
              x  = mirror * x
              xp = mirror * xp

              ! Last do rotation into collimator frame
              x_in(j)  =  x*cRRot +  z*sRRot
              y_in(j)  =  z*cRRot -  x*sRRot
              xp_in(j) = xp*cRRot + zp*sRRot
              yp_in(j) = zp*cRRot - xp*sRRot

            ! Log output energy + nucleons as per the FLUKA coupling
            ! Do not log dead particles
              nnuc1       = nnuc1 + naa(j)                          ! outcoming nucleons
              ien1        = ien1  + rcp(j) * c1e3                   ! outcoming energy

              if(is_crystal) then
                p_in(j) = p
                s_in(j) = s_in(j) + s
              else
                p_in(j) = (one + dpop) * p0
                s_in(j) = sp + (real(j_slices,fPrec)-one) * c_length
              end if
            else
              x_in(j) = x
              y_in(j) = z
            end if

        end if

    end if

end subroutine










subroutine pyk2_finish(num_particles, &
                    x_particles, &
                    xp_particles, &
                    y_particles, &
                    yp_particles, &
                    s_particles, &
                    p_particles)
  use coll_common ,      only : rnd_seed, rcx, rcxp, rcy, rcyp, rcp, rcs, coll_expandArrays

  implicit none


  ! ############################
  ! ## variables declarations ##
  ! ############################

  integer, intent(in)          :: num_particles
  real(kind=8), intent(inout)  :: x_part(num_particles)
  real(kind=8), intent(inout)  :: xp_part(num_particles)
  real(kind=8), intent(inout)  :: y_part(num_particles)
  real(kind=8), intent(inout)  :: yp_part(num_particles)
  real(kind=8), intent(inout)  :: s_part(num_particles)
  real(kind=8), intent(inout)  :: p_part(num_particles)

  do j=1,npart
     x_particles(j) = rcx(j)
     xp_particles(j) = rcxp(j)
     y_particles(j) = rcy(j)
     yp_particles(j) = rcyp(j)
     s_particles(j) = rcs(j)
     p_particles(j) = rcp(j)
  end do
end subroutine 










subroutine pyk2_run(num_particles, &
                    x_particles, &
                    xp_particles, &
                    y_particles, &
                    yp_particles, &
                    s_particles, &
                    p_particles, &
                    part_hit, &
                    part_abs, &
                    part_impact, &
                    part_indiv, &
                    part_linteract, &
                    nhit_stage, &
                    nabs_type, &
                    linside, &
                    run_exenergy, &
                    run_anuc, &
                    run_zatom, &
                    run_emr, &
                    run_rho, &
                    run_hcut, &
                    run_bnref, &
                    run_csref0, &
                    run_csref1, &
                    run_csref4, &
                    run_csref5, &
                    run_radl, &
                    run_dlri, & 
                    run_dlyi, &
                    run_eUm, &
                    run_ai, &
                    run_collnt, &
                    run_cprob, &
                    run_xintl, &
                    run_bn, &
                    run_ecmsq, &
                    run_xln15s, &
                    run_bpp, &
                    is_crystal, &
                    c_length, &
                    c_rotation, &
                    c_aperture, &
                    c_offset, &
                    c_tilt, &
                    c_enom, &
                    onesided, &
                    random_generator_seed)

  use floatPrecision
  use numerical_constants

  use parpro ,           only : npart
  use mod_alloc ,        only : alloc      ! to allocate partID etc
  use mod_common ,       only : iexact, napx, unit208, aa0
  use mod_common_main ,  only : partID, parentID, pairID, naa
  use mod_ranlux ,       only : rluxgo     ! for ranlux init

  use coll_common ,      only : rnd_seed, rcx, rcxp, rcy, rcyp, rcp, rcs, coll_expandArrays
  ! //use coll_materials ! for collmat_init
  use coll_k2        ! for scattering

  implicit none


  ! ############################
  ! ## variables declarations ##
  ! ############################

  integer, intent(in)          :: num_particles
  real(kind=8), intent(inout)  :: x_particles(num_particles)
  real(kind=8), intent(inout)  :: xp_particles(num_particles)
  real(kind=8), intent(inout)  :: y_particles(num_particles)
  real(kind=8), intent(inout)  :: yp_particles(num_particles)
  real(kind=8), intent(inout)  :: s_particles(num_particles)
  real(kind=8), intent(inout)  :: p_particles(num_particles)

  integer(kind=4)  , intent(inout) :: part_hit(num_particles)
  integer(kind=4)  , intent(inout) :: part_abs(num_particles)
  real(kind=8) , intent(inout) :: part_impact(num_particles)
  real(kind=8) , intent(inout) :: part_indiv(num_particles)
  real(kind=8) , intent(inout) :: part_linteract(num_particles)
  integer(kind=4)  , intent(inout) :: nhit_stage(num_particles)
  integer(kind=4)  , intent(inout) :: nabs_type(num_particles)
  logical(kind=4)  , intent(inout) :: linside(num_particles)

  real(kind=8)     , intent(inout) :: run_exenergy
  real(kind=8)     , intent(in) :: run_anuc
  real(kind=8)     , intent(in) :: run_zatom
  real(kind=8)     , intent(in) :: run_emr
  real(kind=8)     , intent(in) :: run_rho
  real(kind=8)     , intent(in) :: run_hcut
  real(kind=8)     , intent(in) :: run_bnref

  real(kind=8)     , intent(in) :: run_csref0
  real(kind=8)     , intent(in) :: run_csref1
  real(kind=8)     , intent(in) :: run_csref4
  real(kind=8)     , intent(in) :: run_csref5

  real(kind=8)     , intent(in) :: run_radl
  real(kind=8)     , intent(in) :: run_dlri
  real(kind=8)     , intent(in) :: run_dlyi
  real(kind=8)     , intent(in) :: run_eUm
  real(kind=8)     , intent(in) :: run_ai
  real(kind=8)     , intent(in) :: run_collnt
  real(kind=8)     , intent(in) :: run_cprob(0:5)
  real(kind=8)     , intent(in) :: run_xintl
  real(kind=8)     , intent(inout) :: run_bn
  real(kind=8)     , intent(in) :: run_ecmsq
  real(kind=8)     , intent(in) :: run_xln15s
  real(kind=8)     , intent(in) :: run_bpp

  logical(kind=4)  , intent(in) :: is_crystal
  real(kind=8) ,    intent(in) :: c_length
  real(kind=8) ,    intent(in) :: c_rotation
  real(kind=8) ,    intent(in) :: c_aperture
  real(kind=8) ,    intent(in) :: c_offset
  real(kind=8) , intent(inout) :: c_tilt(2)
  real(kind=8) ,    intent(in) :: c_enom
  logical(kind=4) ,  intent(in):: onesided
  integer, intent(in)          :: random_generator_seed

  integer j



  ! ####################
  ! ## initialisation ##
  ! ####################
  !character(len=:),    allocatable   :: numpart
  !numpart="20000"
  !read(numpart,*) napx
  npart=num_particles

  if(random_generator_seed .ge. 0) then
        call rluxgo(3, random_generator_seed, 0, 0)
  end if

  do j=1,npart
    naa(j) = aa0
    partID(j)   = j
    parentID(j) = j
    pairID(1,j) = (j+1)/2    ! The pairID of particle j
    pairID(2,j) = 2-mod(j,2) ! Either particle 1 or 2 of the pair
  end do
  
  napx=npart  ! this decreases after absorptions!
  unit208=109 

  do j=1,npart
    rcx(j) = x_particles(j)
    rcxp(j) = xp_particles(j)
    rcy(j) = y_particles(j)
    rcyp(j) = yp_particles(j)
    rcs(j) = s_particles(j)
    rcp(j) = p_particles(j)
  end do

  call k2coll_collimate( &
     run_exenergy, run_anuc, run_zatom, run_emr, run_rho, run_hcut, run_bnref, &
     run_csref0, run_csref1, run_csref4, run_csref5, run_radl, run_dlri, &
     run_dlyi, run_eUm, run_ai, run_collnt, run_cprob, run_xintl, run_bn, &
     run_ecmsq, run_xln15s, run_bpp, is_crystal, &
     c_length, c_rotation, c_aperture, c_offset, c_tilt, &
     rcx, rcxp, rcy, rcyp, rcp, rcs, &
     c_enom, part_hit, part_abs, &
     part_impact, part_indiv, part_linteract, &
     onesided, nhit_stage, 1, nabs_type, linside)

  do j=1,npart
     x_particles(j) = rcx(j)
     xp_particles(j) = rcxp(j)
     y_particles(j) = rcy(j)
     yp_particles(j) = rcyp(j)
     s_particles(j) = rcs(j)
     p_particles(j) = rcp(j)
  end do
end subroutine 
