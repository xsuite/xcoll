
! Manually copied this from other source files from SixTrack to avoid having to compile a bunch of dependencies
subroutine updatePairMap

  use parpro, only : npart
  use mod_common_main, only : pairID, pairMap

  integer j

  pairMap(:,:) = 0
  do j=1,npart
    ! do not update the map in case the particle is not a primary one
    if (pairID(1,j)==0.and.pairID(2,j)==0) cycle
    pairMap(pairID(2,j),pairID(1,j)) = j
  end do

end subroutine updatePairMap
