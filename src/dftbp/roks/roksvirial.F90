module dftbp_roks_roksvirial

  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_scc, only : TScc
  use dftbp_io_message, only : error, warning
  use dftbp_math_blasroutines, only : hemv, symm
  use dftbp_roks_roks, only : TRoksCalc
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_constants, only : Hartree__eV

  implicit none
  private

  public :: calculateRoksVirialIntegral

contains

  subroutine calculateRoksVirialIntegral(roks, denseDesc, sccCalc, &
      & iNeighbour, img2CentCell)
    
    type(TRoksCalc), intent(inout) :: roks
    type(TDenseDescr), intent(in) :: denseDesc
    type(TScc), intent(in) :: sccCalc
    integer, intent(in) :: iNeighbour(0:,:)
    integer, intent(in) :: img2CentCell(:)

    integer :: iOrb, fOrb, nOrb, nAtom
    integer :: iAtom, iFirst, iLast
    real(dp) :: orbitalOverlap
    real(dp), allocatable :: overlapTimesCoefficients(:,:)
    real(dp), allocatable :: qOrb(:)
    real(dp), allocatable :: gammaMatrix(:,:)
    real(dp), allocatable :: gammaTimesQ(:)

    roks%virialIntegral = 0.0_dp
    roks%virialTransitionChargeSum = 0.0_dp
    roks%virialIntegralAvailable = .false.

    if (allocated(roks%virialTransitionCharges)) then
      deallocate(roks%virialTransitionCharges)
    end if

    if (roks%No /= 2) then
      call warning("Becke virial integral currently requires exactly two open-shell orbitals")
      return
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
    allocate(roks%virialTransitionCharges(nAtom))

    do iAtom = 1, nAtom
      iFirst = denseDesc%iAtomStart(iAtom)
      iLast = denseDesc%iAtomStart(iAtom + 1) - 1
      roks%virialTransitionCharges(iAtom) = sum(qOrb(iFirst:iLast))
    end do

    roks%virialTransitionChargeSum = sum(roks%virialTransitionCharges)
    orbitalOverlap = dot_product(roks%coefficients(:,iOrb), overlapTimesCoefficients(:,fOrb))

    if (abs(roks%virialTransitionChargeSum - orbitalOverlap) > 1.0e-10_dp) then
      call warning("Inconsistent ROKS transition-charge sum")
    end if

    if (abs(roks%virialTransitionChargeSum) > 1.0e-8_dp) then
      call warning("ROKS virial transition charges are not neutral")
    end if

    allocate(gammaMatrix(nAtom,nAtom))
    call sccCalc%getAtomicGammaMatrix(gammaMatrix, iNeighbour, img2CentCell)

    allocate(gammaTimesQ(nAtom))
    call hemv(gammaTimesQ, gammaMatrix, roks%virialTransitionCharges)

    roks%virialIntegral = dot_product(roks%virialTransitionCharges, gammaTimesQ)

    roks%virialIntegralAvailable = .true.

    if (roks%writeDiagnostics) then
      write(stdOut, "(A,2(1X,I0))") "--> ROKS virial transition orbitals:", iOrb, fOrb

      write(stdOut, "(A,1X,ES20.12)") "--> ROKS virial transition-charge sum:", &
          & roks%virialTransitionChargeSum
      
      write(stdOut, "(A,1X,ES20.12,A)") "--> ROKS virial integral:", &
          & roks%virialIntegral * Hartree__eV, " eV"
    end if
    
  end subroutine calculateRoksVirialIntegral

end module dftbp_roks_roksvirial
