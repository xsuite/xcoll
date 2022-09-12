subroutine pyk2_init(random_generator_seed)
  use mod_ranlux ,       only : rluxgo     ! for ranlux init

  implicit none

  integer, intent(inout)          :: random_generator_seed

  ! Initialize random number generator
  !if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(random_generator_seed <  0) random_generator_seed = abs(random_generator_seed)
  call rluxgo(3, random_generator_seed, 0, 0)

end subroutine


subroutine make_ruth_dist(cgen, zatom, emr, hcut, cnorm)
  use mod_funlux ,       only : funlxp
  use coll_k2,           only : zatom_curr, emr_curr, cnorm_curr, k2coll_ruth

  implicit none

  real(kind=8), intent(inout) :: cgen(200)
  real(kind=8), intent(in)    :: zatom
  real(kind=8), intent(in)    :: emr
  real(kind=8), intent(in)    :: hcut
  real(kind=8), intent(in)    :: cnorm
  real(kind=8), parameter     :: tlcut = 0.0009982

  zatom_curr = zatom
  emr_curr = emr
  cnorm_curr = cnorm
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

