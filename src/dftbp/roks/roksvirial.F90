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
  !> Calculate the Becke virial integral between the two ROKS open orbitals
  !>
  !> Atomic Mulliken transition charges q_A^{ia} are constructed from the
  !> open-orbital coefficients and overlap matrix. The integral is approximated as
  !>
  !>   K_ia = sum_AB q_A^{ia} gamma_AB q_B^{ia}
  !>
  !> where gamma_AB is built with the Hartree-only Hubbard parameters
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

    integer :: iOrb, fOrb, nOrb, nAtom
    integer :: iAtom, iFirst, iLast
    integer :: iLocal, jLocal, mu, nu
    integer :: iShell, iShellFirst, iShellLast
    integer :: iSp
    real(dp) :: transitionChargeSum
    real(dp) :: degeneracy, partialTrace
    real(dp) :: virialMullikenIntegral
    real(dp) :: virialOnsiteIntegral
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

    iOrb = roks%Nc + 1
    fOrb = roks%Nc + 2

    allocate(overlapTimesCoefficients, mold=roks%coefficients)
    call symm(overlapTimesCoefficients, "L", roks%overlap, roks%coefficients)

    nOrb = size(roks%coefficients, dim=1)
    allocate(qOrb(nOrb))

    qOrb(:) = 0.5_dp * ( roks%coefficients(:,iOrb) * overlapTimesCoefficients(:,fOrb) &
        & + roks%coefficients(:,fOrb) * overlapTimesCoefficients(:,iOrb))

    nAtom = size(denseDesc%iAtomStart) - 1
    allocate(qAtom(nAtom))

    do iAtom = 1, nAtom
      iFirst = denseDesc%iAtomStart(iAtom)
      iLast = denseDesc%iAtomStart(iAtom + 1) - 1
      qAtom(iAtom) = sum(qOrb(iFirst:iLast))
    end do

    transitionChargeSum = sum(qAtom)
    if (abs(transitionChargeSum) > 1.0e-8_dp) then
      call error("ROKS virial transition charges do not sum to zero.")
    end if

    allocate(gammaMatrix(nAtom, nAtom))
    call sccCalc%getAtomicGammaMatU(gammaMatrix, roks%hHubbard, species, iNeighbour, img2CentCell)

    allocate(gammaTimesQ(nAtom))
    call hemv(gammaTimesQ, gammaMatrix, qAtom)

    virialMullikenIntegral = dot_product(qAtom, gammaTimesQ)
    virialOnsiteIntegral = 0.0_dp

    if (allocated(roks%virialOnSiteElements)) then
      allocate(qBlock(orb%mOrb, orb%mOrb))
      allocate(onsite(orb%mOrb, orb%mOrb, 2))
      allocate(onsiteTimesQ(orb%mOrb, orb%mOrb))

      do iAtom = 1, nAtom
        iSp = species(iAtom)
        iFirst = denseDesc%iAtomStart(iAtom)
        iLast = denseDesc%iAtomStart(iAtom + 1) - 1
        nOrb = iLast - iFirst + 1

        qBlock(:,:) = 0.0_dp

        do iLocal = 1, nOrb
          mu = iFirst + iLocal - 1

          do jLocal = iLocal, nOrb
            nu = iFirst + jLocal - 1

            qBlock(iLocal,jLocal) = 0.25_dp * (&
                & roks%coefficients(mu,iOrb)&
                & * overlapTimesCoefficients(nu,fOrb)&
                & + roks%coefficients(mu,fOrb)&
                & * overlapTimesCoefficients(nu,iOrb)&
                & + roks%coefficients(nu,iOrb)&
                & * overlapTimesCoefficients(mu,fOrb)&
                & + roks%coefficients(nu,fOrb)&
                & * overlapTimesCoefficients(mu,iOrb))

            qBlock(jLocal,iLocal) = qBlock(iLocal,jLocal)
          end do
        end do
        call getOnsME(orb, iSp, roks%virialOnSiteElements, nOrb, onsite)
        onsiteTimesQ(:,:) = 0.0_dp
        onsiteTimesQ(:nOrb,:nOrb) = qBlock(:nOrb,:nOrb) *&
            & (onsite(:nOrb,:nOrb,1) + onsite(:nOrb,:nOrb,2))
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
            onsiteTimesQ(iLocal,iLocal) = onsiteTimesQ(iLocal,iLocal)&
                & - partialTrace
          end do
        end do
        virialOnsiteIntegral = virialOnsiteIntegral&
            & + sum(qBlock(:nOrb,:nOrb) * onsiteTimesQ(:nOrb,:nOrb))
      end do
    end if
    roks%virialIntegral = virialMullikenIntegral + virialOnsiteIntegral

    if (roks%writeDiagnostics) then
      write(stdOut, "(A,2(1X,I0))") "--> ROKS virial transition orbitals:", iOrb, fOrb
      write(stdOut, "(A,1X,ES20.12)") "--> ROKS virial transition-charge sum:",&
          & transitionChargeSum
      write(stdOut, "(A,1X,ES20.12,A)") "--> ROKS virial Mulliken integral:",&
          & virialMullikenIntegral * Hartree__eV, " eV"

      write(stdOut, "(A,1X,ES20.12,A)") "--> ROKS virial onsite correction:",&
          & virialOnsiteIntegral * Hartree__eV, " eV"

      write(stdOut, "(A,1X,ES20.12,A)") "--> ROKS virial total integral:",&
          & roks%virialIntegral * Hartree__eV, " eV"
    end if

  end subroutine calculateRoksVirialIntegral

end module dftbp_roks_roksvirial
