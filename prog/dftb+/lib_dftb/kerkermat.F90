!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Routines implementing the overlap and gamma matrix for the Kerker convergence accelerator.
module kerkermat
#include "assert.h"
#include "allocate.h"
  use accuracy
  use constants
  use message
  use periodic, only : TNeighborList, getNrOfNeighbors
  use SCC, only : getShellGammaMatrix
  use commontypes
  implicit none
  private

  public :: TKerkerInput, initGtoMatrices, buildGtoMatrices, reduceDq, expandDq, &
             & destructGtoMatrices, infoKerkerLimits, getKerkerMatrices

  !> Input for the kerkerMat module
  type TKerkerInput 
    
    !> Number of atoms
    integer :: nAtom
    !> Hubbard U values. Shape: [nShell, nSpecies]
    real(dp), allocatable :: hubbUs(:,:)
    !> Spin channels
    integer :: nSpin

  end type TKerkerInput

  !! Private module variables (suffixed with "_" for clarity)
  logical :: tInitialised_ = .false.      ! If module is initialised
  integer :: nGauss_  ! Dimension of Kerker matrices (number of unique Hubbards per shell)
  integer :: nAtom_
  integer :: nSpin_
  integer :: nFullOrbResolved_ ! Dimension of mix vector 
  integer, allocatable :: nHubbU_(:)   ! Nr. unique U per species
  integer, allocatable :: iHubbU_(:,:) ! Points to index of unique Hubbard   
  integer, allocatable :: indMat_(:,:) ! Index array for Kerker matrices
  integer, allocatable :: nShell_(:)   ! Number of shells per species
  integer, allocatable :: species_(:)
  real(dp), allocatable :: uniqHubbU_(:,:) ! kerkermat has its own set of unique Hubbards!
  real(dp), allocatable :: sMat_(:,:) ! Overlap of normalized DFTB densities
  real(dp), allocatable :: lMat_(:,:) ! Angular momentum resolved Gamma

contains
  subroutine initGtoMatrices(inp, species, orb, nGauss)
  !!* Computes the dimension of Kerker matrices by identifying the unique 
  !!* Hubbard-Us. Note that this is a non-permanent hack. The rest of the 
  !!* code evaluates the unique U only for spin-unpolarized systems. Similar
  !!* U lead to problems in the LU factorization of the Kerker matrices 
  !!* (linear dependency problem). 
  !!* @param inp Input data (nAtoms, Hubbards, spin).
  !!* @param species Species.
  !!* @param orb Information about the orbitals in the system.
  !!* @param nGauss Dimension of Kerker matrices
    type(TKerkerInput), intent(inout) :: inp 
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    integer, intent(out) :: nGauss
    integer :: mShell, nSpecies, iSp, iL 
    integer :: iU, iAt

    nAtom_ = inp%nAtom
    mShell = orb%mShell
    nSpecies = size(orb%nOrbSpecies)

    ALLOCATE_(nHubbU_, (nSpecies))
    ALLOCATE_(iHubbU_, (mShell, nSpecies))
    ALLOCATE_(uniqHubbU_, (mShell, nSpecies))
    ALLOCATE_(nShell_, (nSpecies))

    nSpin_ = inp%nSpin
    nShell_(:) = orb%nShell
    iHubbU_(:,:) = 0
    iHubbU_(1,:) = 1
    nHubbU_(:) = 1
    uniqHubbU_(:,:) = 0.0_dp
    uniqHubbU_(1,:) = inp%hubbUs(1,:)
    do iSp = 1, nSpecies
      do iL = 2, nShell_(iSp)
        do iU = 1, nHubbU_(iSp)
          if (abs(inp%hubbUs(iL,iSp) - uniqHubbU_(iU,iSp)) < MinHubDiff) then
            iHubbU_(iL,iSp) = iU
            exit
          end if
        end do
        if (iHubbU_(iL,iSp) == 0) then
          nHubbU_(iSp) = nHubbU_(iSp) + 1
          uniqHubbU_(nHubbU_(iSp),iSp) = inp%hubbUs(iL,iSp)
          iHubbU_(iL,iSp) = nHubbU_(iSp)
        end if
      end do
    end do
     
    nFullOrbResolved_ = 0
    do iAt = 1, nAtom_
      iSp = species(iAt)
      do iL = 1, nShell_(iSp)
        nFullOrbResolved_ = nFullOrbResolved_ + 1
      enddo
    enddo  

    nGauss  = 0
    do iAt = 1, nAtom_
      iSp = species(iAt)
      do iU = 1, nHubbU_(iSp)
        nGauss = nGauss + 1
      enddo
    enddo  
    nGauss_ = nGauss 

    ALLOCATE_(sMat_, (nGauss, nGauss))
    ALLOCATE_(lMat_, (nGauss, nGauss))

    tInitialised_ = .true.

  end subroutine initGtoMatrices
 
 !!* Builds the overlap matrix S_ab = \Int n(r-Ra) n(r-Rb) d3r with 
  !!* n(r) = (\tau^3/8Pi) exp(-\tau r), \tau being either atomic specific
  !!* or angular momentum resolved. Also builds the corresponding Gamma matrix.
  !!* For tPeriodic the corresponding lattice sums are performed.   
  !!* @param inp Input data (nAtoms, Hubbards, spin).
  !!* @param neighList Neighbour list.
  !!* @param orb Information about the orbitals in the system.
  !!* @param img2CentCell Array to relate image atoms outside the central cell.
  !!* @param coord Atomic coordinates.
  !!* @param species Species.
 subroutine buildGtoMatrices(inp, neighList, orb, img2CentCell, &
      & coord, species)
    type(TKerkerInput), intent(in) :: inp 
    type(TNeighborList), intent(in) :: neighList
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: species(:)
    real(dp), intent(in) :: coord(:,:)

    !> Local variables
    integer :: mShell, nSpecies, iSp, iU, iAt, iCnt
    integer :: iSp1, iSp2, iU1, iU2, iMat1, iMat2
    integer ::  iAt1, iAt2, iAt2f, iNeigh
    real(dp) :: rA(3), rB(3), zA, zB  

    ASSERT(tInitialised_)
    mShell = orb%mShell
    nSpecies = size(orb%nOrbSpecies)
    if( .not. allocated(indMat_) ) then
      ALLOCATE_(indMat_, (mShell, nAtom_))
      ALLOCATE_(species_, (size(species)))
    endif
    species_(:) = species(:)

    if(inp%nSpin == 4) call error('Kerker does not yet work with Spin-Orbit')

    iCnt = 1
    indMat_(:,:) = 0
    do iAt = 1, nAtom_
      iSp = species_(iAt)
      do iU = 1, nHubbU_(iSp)
        indMat_(iU, iAt) = iCnt
        iCnt = iCnt + 1
      enddo
    enddo
    
    sMat_(:,:) = 0.0_dp
    do iAt1 = 1, nAtom_
      rA(:) =  coord(:,iAt1)
      iSp1 = species_(iAt1)
      do iU1 = 1, nHubbU_(iSp1)
        zA = 3.2_dp*uniqHubbU_(iU1,iSp1) ! Using \tau = 16/5 U 
        iMat1 = indMat_(iU1, iAt1) 
        do iNeigh = 0, neighList%nNeighbor(iAt1)
          iAt2 = neighList%iNeighbor(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          rB(:) =  coord(:,iAt2)
          iSp2 = species_(iAt2f)
          do iU2 = 1, nHubbU_(iSp2)
            zB = 3.2_dp*uniqHubbU_(iU2,iSp2)
            iMat2 = indMat_(iU2, iAt2f) 
            sMat_(iMat1,iMat2) = sMat_(iMat1,iMat2) + sGTO(zA, rA, zB, rB)
          enddo
        enddo
      enddo
    enddo

    call getShellGammaMatrix(lMat_, neighList%iNeighbor, img2CentCell, &
        & coord, species_, nHubbU_, uniqHubbU_, indMat_)
   
    ! Symmetrize matrices
    do iMat1 = 1, maxval(indMat_)
       do iMat2 = iMat1+1, maxval(indMat_)
          sMat_(iMat2, iMat1) = sMat_(iMat1, iMat2)
          lMat_(iMat1, iMat2) = lMat_(iMat2, iMat1)
       enddo
    enddo

  end subroutine buildGtoMatrices

  !!* Computes the overlap matrix S_ab = \Int n(r-Ra) n(r-Rb) d3r by approximating
  !!* the STO by a single GTO. 
  !!* @param zA STO decay constant \tau for orbital/atom A.
  !!* @param rA Position vector for atom A.
  !!* @param zB STO decay constant \tau for orbital/atom B.
  !!* @param rB Position vector for atom B.
  !!* @param iGTO Resulting overlap.
  function sGTO(zA, rA, zB, rB) result (iGTO)
    real(dp),    intent(in)  :: zA
    real(dp),    intent(in)  :: zB
    real(dp),    intent(in)  :: rA(:)
    real(dp),    intent(in)  :: rB(:)
    real(dp), parameter      :: sto2Gto = 0.270104_dp
    real(dp)                 :: alpA, alpB, alpP
    real(dp)                 :: iGTO, rAB2, vect(3), fac

    alpA = sto2Gto*zA*zA
    alpB = sto2Gto*zB*zB
    vect(:) = rA(:)-rB(:)
    rAB2 = sum(vect(:)**2)
    alpP = alpA + alpB
    fac = sqrt(alpA*alpB/(pi*alpP))
    iGTO = exp(-alpA*alpB*rAB2/alpP)*(fac**3)
    return

  end function sGTO

  !!* Reduces dimension of mix vector by summing over entries corresponding
  !!* to equivalent Hubbards
  !!* @param dQin Input vector with dimension nFullOrbResolved_
  !!* @param dQin Output vector with dimension nGauss_
  subroutine reduceDq(dQin, dQout)
    real(dp),    intent(in)  :: dQin(:)
    real(dp),    intent(out)  :: dQout(:)
    integer :: iAt, iL, iSp, iSave, iFull

    ! Quick return if we are not spin-polarized 
    if (size(dQin(:)) == size(dQout(:))) then
      dQout(:) = dQin(:)
      return
    endif

    dQout(:) = 0.0_dp
    iFull = 1
    do iAt = 1, nAtom_
      iSp = species_(iAt)
      iSave = iHubbU_(1, iSp) 
      do iL = 1, nShell_(iSp)
        if(iHubbU_(iL, iSp) == iSave) then
          dQout(indMat_(iHubbU_(iL, iSp), iAt)) = dQout(indMat_(iHubbU_(iL, iSp), iAt)) &
             & +  dQin(iFull)  
        else
          dQout(indMat_(iHubbU_(iL, iSp), iAt)) = dQin(iFull) 
          iSave = iHubbU_(iL, iSp)
        endif
        iFull = iFull + 1
      enddo
    enddo

   end subroutine reduceDq

  !!* Expands dimension of mix vector 
  !!* @param dQin Input vector with dimension nGauss_
  !!* @param dQin Output vector with dimension nFullOrbResolved_
   subroutine expandDq(dQin, dQout)
    real(dp),    intent(in)  :: dQin(:)
    real(dp),    intent(out)  :: dQout(:)
    integer :: iAt, iL, iSp, iSave, iFull

    if (size(dQin(:)) == size(dQout(:))) then
      dQout(:) = dQin(:)
      return
    endif      

    dQout(:) = 0.0_dp
    iFull = 1
    do iAt = 1, nAtom_
      iSp = species_(iAt)
      iSave = 0
      do iL = 1, nShell_(iSp)
        if(iHubbU_(iL, iSp) == iSave) then
          dQout(iFull) = 0.0_dp
        else
          dQout(iFull) = dQin(indMat_(iHubbU_(iL, iSp), iAt)) 
          iSave = iHubbU_(iL, iSp)
        endif
        iFull = iFull + 1
      enddo
    enddo

   end subroutine expandDq

  !!* Provides information on the dimension of the Kerker matrices
  !!* to be used in the mixer routine.
  !!* @param nGauss Dimension of Kerker matrices
  !!* @param nFullOrbResolved Dimension of mix vector
  subroutine infoKerkerLimits(nGauss, nFullOrbResolved)
    integer, intent(out) :: nGauss
    integer, intent(out) :: nFullOrbResolved

    ASSERT(tInitialised_)
    nGauss = nGauss_
    nFullOrbResolved = nFullOrbResolved_

  end subroutine infoKerkerLimits

  !!* Provides the Kerker matrices for the mixer routine
  !!* @param sMat ! Overlap of normalized DFTB densities
  !!* @param lMat ! Angular momentum resolved Gamma
  subroutine getKerkerMatrices(sMat, lMat)
    real(dp)    , intent(out) :: sMat(:,:)
    real(dp)    , intent(out) :: lMat(:,:)

    ASSERT(tInitialised_)
    sMat(:,:) = sMat_(:,:)
    lMat(:,:) = lMat_(:,:)

  end subroutine getKerkerMatrices

  !!* Deallocate local arrays
  subroutine destructGtoMatrices()
    ASSERT(tInitialised_)
    DEALLOCATE_(nHubbU_)
    DEALLOCATE_(iHubbU_)
    DEALLOCATE_(uniqHubbU_)
    DEALLOCATE_(indMat_)
    DEALLOCATE_(nShell_)
    DEALLOCATE_(species_)
    DEALLOCATE_(sMat_)
    DEALLOCATE_(lMat_)
  end subroutine destructGtoMatrices
end module kerkermat
