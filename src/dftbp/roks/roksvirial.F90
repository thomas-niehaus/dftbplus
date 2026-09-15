module dftbp_roks_roksvirial

  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_scc, only : TScc
  use dftbp_io_message, only : error, warning
  use dftbp_math_blasroutines, only : hemv, symm
  use dftbp_roks_roks, only : TRoksCalc
  use dftbp_dftb_onsitecorrection, only : getOnsME
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_constants, only : Hartree__eV

  implicit none
  private

  public :: calculateRoksVirialIntegral

contains
  !> Calculate the occupation-weighted Becke virial integral
  !>
  !> For fractional occupations, the hole and electron characters may be
  !> distributed over several orbitals. The ensemble virial integral is
  !>
  !>   K = sum_he w_h w_e (he|eh),
  !>
  !> where the individual integrals are evaluated from atomic Mulliken
  !> transition charges and optional onsite corrections. The atomic gamma
  !> matrix is constructed with the Hartree-only Hubbard parameters.
  subroutine calculateRoksVirialIntegral(roks, denseDesc, sccCalc, species, iNeighbour,&
      & img2CentCell, orb)

    !> ROKS orbitals and resulting virial integral
    type(TRoksCalc), intent(inout) :: roks

    !> Dense AO-to-atom indexing information
    type(TDenseDescr), intent(in) :: denseDesc

    !> SCC data used to construct the atomic gamma matrix
    type(TScc), intent(in) :: sccCalc

    !> Chemical species index of each central-cell atom
    integer, intent(in) :: species(:)

    !> Neighbour-list indices used to construct gamma
    integer, intent(in) :: iNeighbour(0:,:)

    !> Mapping from periodic images to central-cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital and shell indexing
    type(TOrbitals), intent(in) :: orb

    integer :: iHole, iElectron
    integer :: nOrbTotal, nOrbAtom, nAtom, nClosed
    integer :: iAtom, iFirst, iLast
    integer :: iLocal, jLocal, mu, nu
    integer :: iShell, iShellFirst, iShellLast
    integer :: iSp
    real(dp) :: holeNorm, electronNorm
    real(dp) :: pairWeight
    real(dp) :: transitionChargeSum
    real(dp) :: maxTransitionChargeSum
    real(dp) :: degeneracy, partialTrace
    real(dp) :: pairMullikenIntegral
    real(dp) :: pairOnsiteIntegral
    real(dp) :: virialMullikenIntegral
    real(dp) :: virialOnsiteIntegral
    real(dp), parameter :: weightTolerance = 1.0e-12_dp
    real(dp), allocatable :: spinOccupation(:)
    real(dp), allocatable :: holeWeight(:)
    real(dp), allocatable :: electronWeight(:)
    real(dp), allocatable :: overlapTimesCoefficients(:,:)
    real(dp), allocatable :: qOrb(:)
    real(dp), allocatable :: qAtom(:)
    real(dp), allocatable :: gammaMatrix(:,:)
    real(dp), allocatable :: gammaTimesQ(:)
    real(dp), allocatable :: qBlock(:,:)
    real(dp), allocatable :: onsite(:,:,:)
    real(dp), allocatable :: onsiteTimesQ(:,:)

    if (roks%No /= 2) then
      call error("RoksVirial requires exactly two open-shell orbitals")
    end if

    if (.not. allocated(roks%occupations)) then
      call error("RoksVirial requires ROKS orbital occupations")
    end if

    nOrbTotal = size(roks%coefficients, dim=2)
    nClosed = roks%Nc + roks%No / 2

    if (nClosed < 1 .or. nClosed >= nOrbTotal) then
      call error("Invalid closed-shell boundary for ROKS virial integral")
    end if

    allocate(spinOccupation(nOrbTotal))
    allocate(holeWeight(nOrbTotal), source=0.0_dp)
    allocate(electronWeight(nOrbTotal), source=0.0_dp)

    spinOccupation(:) = roks%occupations(:,1) - roks%occupations(:,2)

    if (any(spinOccupation < -weightTolerance)) then
      call error("ROKS virial requires alpha to be the majority occupation")
    end if

    ! Remove harmless negative roundoff.
    spinOccupation(:) = max(spinOccupation(:), 0.0_dp)

    holeNorm = sum(spinOccupation(1:nClosed))
    electronNorm = sum(spinOccupation(nClosed + 1:nOrbTotal))

    if (holeNorm <= weightTolerance) then
      call error("ROKS virial could not identify the hole subspace")
    end if

    if (electronNorm <= weightTolerance) then
      call error("ROKS virial could not identify the electron subspace")
    end if

    holeWeight(1:nClosed) = spinOccupation(1:nClosed) / holeNorm
    electronWeight(nClosed + 1:nOrbTotal) = &
        & spinOccupation(nClosed + 1:nOrbTotal) / electronNorm

    allocate(overlapTimesCoefficients, mold=roks%coefficients)
    call symm(overlapTimesCoefficients, "L", roks%overlap, roks%coefficients)

    allocate(qOrb(size(roks%coefficients, dim=1)))

    nAtom = size(denseDesc%iAtomStart) - 1
    allocate(qAtom(nAtom))

    allocate(gammaMatrix(nAtom, nAtom))
    call sccCalc%getAtomicGammaMatU(gammaMatrix, roks%hHubbard, species, iNeighbour, img2CentCell)

    allocate(gammaTimesQ(nAtom))

    if (allocated(roks%virialOnSiteElements)) then
      allocate(qBlock(orb%mOrb, orb%mOrb))
      allocate(onsite(orb%mOrb, orb%mOrb, 2))
      allocate(onsiteTimesQ(orb%mOrb, orb%mOrb))
    end if

    virialMullikenIntegral = 0.0_dp
    virialOnsiteIntegral = 0.0_dp
    maxTransitionChargeSum = 0.0_dp

    if (roks%writeDiagnostics) then
      write(stdOut, "(A)") "--> ROKS virial hole orbital weights"
      do iHole = 1, nClosed
        if (holeWeight(iHole) > weightTolerance) then
          write(stdOut, "(A,I6,F14.8)") "--> ROKS virial hole:", iHole, holeWeight(iHole)
        end if
      end do

      write(stdOut, "(A)") "--> ROKS virial electron orbital weights"
      do iElectron = nClosed + 1, nOrbTotal
        if (electronWeight(iElectron) > weightTolerance) then
          write(stdOut, "(A,I6,F14.8)") "--> ROKS virial electron:", &
              & iElectron, electronWeight(iElectron)
        end if
      end do
    end if

    do iHole = 1, nClosed
      if (holeWeight(iHole) <= weightTolerance) then
        cycle
      end if

      do iElectron = nClosed + 1, nOrbTotal
        if (electronWeight(iElectron) <= weightTolerance) then
          cycle
        end if

        pairWeight = holeWeight(iHole) * electronWeight(iElectron)

        qOrb(:) = 0.5_dp * (&
            & roks%coefficients(:,iHole) * overlapTimesCoefficients(:,iElectron)&
            & + roks%coefficients(:,iElectron) * overlapTimesCoefficients(:,iHole))

        do iAtom = 1, nAtom
          iFirst = denseDesc%iAtomStart(iAtom)
          iLast = denseDesc%iAtomStart(iAtom + 1) - 1
          qAtom(iAtom) = sum(qOrb(iFirst:iLast))
        end do

        transitionChargeSum = sum(qAtom)
        maxTransitionChargeSum = max(maxTransitionChargeSum, abs(transitionChargeSum))

        if (abs(transitionChargeSum) > 1.0e-8_dp) then
          call error("ROKS virial transition charges do not sum to zero")
        end if

        call hemv(gammaTimesQ, gammaMatrix, qAtom)
        pairMullikenIntegral = dot_product(qAtom, gammaTimesQ)
        pairOnsiteIntegral = 0.0_dp

        if (allocated(roks%virialOnSiteElements)) then
          do iAtom = 1, nAtom
            iSp = species(iAtom)
            iFirst = denseDesc%iAtomStart(iAtom)
            iLast = denseDesc%iAtomStart(iAtom + 1) - 1
            nOrbAtom = iLast - iFirst + 1

            qBlock(:,:) = 0.0_dp

            do iLocal = 1, nOrbAtom
              mu = iFirst + iLocal - 1

              do jLocal = iLocal, nOrbAtom
                nu = iFirst + jLocal - 1

                qBlock(iLocal,jLocal) = 0.25_dp * (&
                    & roks%coefficients(mu,iHole)&
                    & * overlapTimesCoefficients(nu,iElectron)&
                    & + roks%coefficients(mu,iElectron)&
                    & * overlapTimesCoefficients(nu,iHole)&
                    & + roks%coefficients(nu,iHole)&
                    & * overlapTimesCoefficients(mu,iElectron)&
                    & + roks%coefficients(nu,iElectron)&
                    & * overlapTimesCoefficients(mu,iHole))

                qBlock(jLocal,iLocal) = qBlock(iLocal,jLocal)
              end do
            end do

            call getOnsME(orb, iSp, roks%virialOnSiteElements, nOrbAtom, onsite)
            onsiteTimesQ(:,:) = 0.0_dp
            onsiteTimesQ(:nOrbAtom,:nOrbAtom) = qBlock(:nOrbAtom,:nOrbAtom) *&
                & (onsite(:nOrbAtom,:nOrbAtom,1) + onsite(:nOrbAtom,:nOrbAtom,2))

            do iShell = 1, orb%nShell(iSp)
              iShellFirst = orb%posShell(iShell, iSp)
              iShellLast = orb%posShell(iShell + 1, iSp) - 1

              degeneracy = real(2 * orb%angShell(iShell, iSp) + 1, dp)

              partialTrace = 0.0_dp
              do iLocal = iShellFirst, iShellLast
                partialTrace = partialTrace + onsiteTimesQ(iLocal,iLocal)
              end do
              partialTrace = partialTrace / degeneracy

              do iLocal = iShellFirst, iShellLast
                onsiteTimesQ(iLocal,iLocal) = onsiteTimesQ(iLocal,iLocal) - partialTrace
              end do
            end do

            pairOnsiteIntegral = pairOnsiteIntegral&
                & + sum(qBlock(:nOrbAtom,:nOrbAtom) * onsiteTimesQ(:nOrbAtom,:nOrbAtom))
          end do
        end if

        virialMullikenIntegral = virialMullikenIntegral + pairWeight * pairMullikenIntegral
        virialOnsiteIntegral = virialOnsiteIntegral + pairWeight * pairOnsiteIntegral

        if (roks%writeDiagnostics) then
          write(stdOut, "(A,2(1X,I0),A,F12.8,2(A,ES20.12))") &
              & "--> ROKS virial pair:", iHole, iElectron, ", weight ", pairWeight,&
              & ", Mulliken [eV] ", pairMullikenIntegral * Hartree__eV,&
              & ", onsite [eV] ", pairOnsiteIntegral * Hartree__eV
        end if
      end do
    end do

    roks%virialIntegral = virialMullikenIntegral + virialOnsiteIntegral

    if (roks%writeDiagnostics) then
      write(stdOut, "(A,1X,ES20.12)") "--> ROKS virial maximum transition-charge error:",&
          & maxTransitionChargeSum
      write(stdOut, "(A,1X,ES20.12,A)") "--> ROKS virial Mulliken integral:",&
          & virialMullikenIntegral * Hartree__eV, " eV"

      write(stdOut, "(A,1X,ES20.12,A)") "--> ROKS virial onsite correction:",&
          & virialOnsiteIntegral * Hartree__eV, " eV"

      write(stdOut, "(A,1X,ES20.12,A)") "--> ROKS virial total integral:",&
          & roks%virialIntegral * Hartree__eV, " eV"
    end if

  end subroutine calculateRoksVirialIntegral

end module dftbp_roks_roksvirial
