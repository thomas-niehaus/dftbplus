!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include 'common.fypp'

!> Restricted open-shell Kohn-Sham data and routines.
!>
!> Constructs a common set of spatial orbitals from conventional alpha-
!> and beta-spin DFTB Hamiltonians. The effective Hamiltonian uses the
!> Roothaan core/open/virtual block structure employed for high-spin
!> restricted open-shell calculations.
!>
!> The present implementation supports real, dense, Gamma-point matrices.
module dftbp_roks_roks

  use dftbp_io_message, only : error
  use dftbp_common_accuracy, only : dp, elecTolMax

  implicit none

  private
  public :: TRoksCalc

  !> Runtime data for a restricted open-shell Kohn-Sham calculation.
  type :: TRoksCalc

    !> Number of doubly occupied core orbitals.
    integer :: Nc = 0

    !> Number of singly occupied open-shell orbitals.
    integer :: No = 0

    !> Number of unoccupied virtual orbitals.
    integer :: Nv = 0

    !> Whether internal ROKS diagnostics should be printed.
    logical :: writeDiagnostics = .false.

    !> Maximum number of fixed-Hamiltonian orbital iterations.
    integer :: maxIterations = 50

    !> Convergence tolerance for the orbital-stationarity residual.
    real(dp) :: tolerance = 1.0e-8_dp

    !> Conventional alpha-spin Hamiltonian in the AO basis.
    real(dp), allocatable :: hamAlpha(:,:)

    !> Conventional beta-spin Hamiltonian in the AO basis.
    real(dp), allocatable :: hamBeta(:,:)

    !> Alpha-spin Hamiltonian in the shared MO basis.
    real(dp), allocatable :: hamAlphaMo(:,:)

    !> Beta-spin Hamiltonian in the shared MO basis.
    real(dp), allocatable :: hamBetaMo(:,:)

    !> Effective ROKS Hamiltonian.
    real(dp), allocatable :: hamEffective(:,:)
    
    !> Original AO overlap matrix.
    real(dp), allocatable :: overlap(:,:)

    !> Effective ROKS Hamiltonian in the MO basis.
    real(dp), allocatable :: hamEffectiveMo(:,:)

    !> Shared ROKS molecular-orbital coefficients.
    real(dp), allocatable :: coefficients(:,:)

    !> Eigenvalues of the effective ROKS Hamiltonian.
    real(dp), allocatable :: eigenvalues(:)
    
    !> Atomic transition charges between the two open-shell orbitals.
    real(dp), allocatable :: virialTransitionCharges(:)

    !> Becke virial two-electron integral in Hartree.
    real(dp) :: virialIntegral = 0.0_dp
    
    !> Sum of the atomic transition charges.
    real(dp) :: virialTransitionChargeSum = 0.0_dp

    !> Whether the integral has been evaluated for the current geometry.
    logical :: virialIntegralAvailable = .false.

  contains

    !> Initialize the core, open-shell and virtual orbital counts.
    procedure :: init => TRoksCalc_init

    !> Allocate dense ROKS Hamiltonian and orbital storage.
    procedure :: allocateMatrices => TRoksCalc_allocateMatrices

    !> Build the spin-averaged Hamiltonian used to obtain trial common orbitals.
    procedure :: buildInitialHamiltonian => TRoksCalc_buildInitialHamiltonian

    !> Transform the alpha and beta Hamiltonians to the current common MO basis.
    procedure :: transformHamiltoniansToMo => TRoksCalc_transformHamiltoniansToMo

    !> Construct the effective ROKS Hamiltonian and transform it to the AO basis.
    procedure :: buildEffectiveHamiltonian => TRoksCalc_buildEffectiveHamiltonian

    !> Insert the spin-dependent core-open and open-virtual MO blocks.
    procedure :: applyMoCouplings => TRoksCalc_applyMoCouplings

    !> Return the largest independent ROKS orbital-gradient element.
    procedure :: getStationarityResidual => TRoksCalc_getStationarityResidual
  end type TRoksCalc


contains

  !> Initialize the high-spin restricted open-shell orbital partition.
  !>
  !> The beta population defines the number of doubly occupied core
  !> orbitals. The excess alpha population defines the number of singly
  !> occupied open-shell orbitals. All remaining orbitals are virtual.
  subroutine TRoksCalc_init(this, nEl, nOrb, maxIterations, tolerance, &
    & writeDiagnostics)
    !> ROKS calculation data.
    class(TRoksCalc), intent(out) :: this

    !> Number of electrons in each spin channel.
    real(dp), intent(in) :: nEl(:)

    !> Number of available spatial orbitals.
    integer, intent(in) :: nOrb

    !> Maximum number of fixed-Hamiltonian ROKS orbital iterations.
    integer, intent(in) :: maxIterations

    !> Convergence tolerance for the orbital-stationarity residual.
    real(dp), intent(in) :: tolerance

    !> Whether internal ROKS diagnostics should be printed.
    logical, intent(in) :: writeDiagnostics

    this%maxIterations = maxIterations
    this%tolerance = tolerance
    this%writeDiagnostics = writeDiagnostics

    if (size(nEl) /= 2) then
      call error("ROKS requires two spin electron populations")
    end if

    if (any(abs(nEl - nint(nEl)) > elecTolMax)) then
      call error("ROKS requires integer spin electron populations")
    end if

    if (nEl(1) < nEl(2)) then
      call error("ROKS currently requires alpha to be the majority spin")
    end if

    ! Orbitals 1:Nc are doubly occupied.
    this%Nc = nint(nEl(2))

    ! Orbitals Nc+1:Nc+No are singly occupied by alpha electrons.
    this%No = nint(nEl(1) - nEl(2))

    ! Remaining orbitals are virtual.
    this%Nv = nOrb - this%Nc - this%No

    if (this%Nc < 0) then
      call error("ROKS has a negative number of core orbitals")
    end if

    if (this%No < 1) then
      call error("ROKS requires at least one open-shell orbital")
    end if

    if (this%Nv < 0) then
      call error("ROKS has more occupied orbitals than available orbitals")
    end if

  end subroutine TRoksCalc_init


  !> Allocate dense ROKS Hamiltonian storage.
  subroutine TRoksCalc_allocateMatrices(this, nRows, nCols)

    !> ROKS calculation data.
    class(TRoksCalc), intent(inout) :: this

    !> Number of locally stored matrix rows.
    integer, intent(in) :: nRows

    !> Number of locally stored matrix columns.
    integer, intent(in) :: nCols

    if (nRows < 1 .or. nCols < 1) then
      call error("Invalid dense matrix dimensions for ROKS")
    end if

    allocate(this%hamAlpha(nRows, nCols), source=0.0_dp)
    allocate(this%hamBeta(nRows, nCols), source=0.0_dp)
    allocate(this%hamAlphaMo(nRows, nCols), source=0.0_dp)
    allocate(this%hamBetaMo(nRows, nCols), source=0.0_dp)
    allocate(this%hamEffective(nRows, nCols), source=0.0_dp)
    allocate(this%overlap(nRows, nCols), source=0.0_dp)
    allocate(this%hamEffectiveMo(nRows, nCols), source=0.0_dp)
    allocate(this%coefficients(nRows, nCols), source=0.0_dp)
    allocate(this%eigenvalues(this%Nc + this%No + this%Nv), source=0.0_dp)

  end subroutine TRoksCalc_allocateMatrices

  !> Build an initial common-orbital Hamiltonian from the alpha and beta Hamiltonians.
  subroutine TRoksCalc_buildInitialHamiltonian(this)

    !> ROKS calculation data.
    class(TRoksCalc), intent(inout) :: this

    @:ASSERT(all(shape(this%hamAlpha) == shape(this%hamEffective)))

    this%hamEffective(:,:) = 0.5_dp * (this%hamAlpha(:,:) + this%hamBeta(:,:))

  end subroutine TRoksCalc_buildInitialHamiltonian
  
  !> Transform the alpha and beta Hamiltonians from the AO basis
  !> to the shared orthonormal ROKS MO basis.
  subroutine TRoksCalc_transformHamiltoniansToMo(this)

    class(TRoksCalc), intent(inout) :: this

    this%hamAlphaMo = matmul(transpose(this%coefficients), matmul(this%hamAlpha, this%coefficients))

    this%hamBetaMo = matmul(transpose(this%coefficients), matmul(this%hamBeta, this%coefficients))

  end subroutine TRoksCalc_transformHamiltoniansToMo

  !> Assemble spin-dependent couplings between the ROKS orbital spaces.
  !>
  !> The alpha/beta mean defines the common part of the effective
  !> Hamiltonian. Couplings between doubly and singly occupied orbitals
  !> are taken from the beta Hamiltonian, while couplings between singly
  !> occupied and empty orbitals are taken from the alpha Hamiltonian.
  !> The reverse blocks are assigned by transposition to retain a real
  !> symmetric matrix.
  !> A related common-orbital construction is used by the PySCF
  !> restricted-open-shell implementation; see pyscf.scf.rohf.
  subroutine TRoksCalc_applyMoCouplings(this)

    class(TRoksCalc), intent(inout) :: this
    integer :: iCoreFirst, iCoreLast
    integer :: iOpenFirst, iOpenLast
    integer :: iVirtualFirst, iVirtualLast
    integer :: nOrb

    nOrb = size(this%hamEffectiveMo, dim=1)

    @:ASSERT(size(this%hamEffectiveMo, dim=2) == nOrb)
    @:ASSERT(this%Nc + this%No + this%Nv == nOrb)

    iCoreFirst = 1
    iCoreLast = this%Nc

    iOpenFirst = this%Nc + 1
    iOpenLast = this%Nc + this%No

    iVirtualFirst = this%Nc + this%No + 1
    iVirtualLast = nOrb

    ! Core-open block: use the beta-spin Hamiltonian.
    if (this%Nc > 0 .and. this%No > 0) then
      this%hamEffectiveMo(iCoreFirst:iCoreLast, iOpenFirst:iOpenLast) = &
          & this%hamBetaMo(iCoreFirst:iCoreLast, iOpenFirst:iOpenLast)

      ! Enforce symmetry explicitly.
      this%hamEffectiveMo(iOpenFirst:iOpenLast, iCoreFirst:iCoreLast) = &
          & transpose(this%hamEffectiveMo(iCoreFirst:iCoreLast, iOpenFirst:iOpenLast))
    end if

    ! Open-virtual block: use the alpha-spin Hamiltonian.
    if (this%No > 0 .and. this%Nv > 0) then
      this%hamEffectiveMo(iOpenFirst:iOpenLast, iVirtualFirst:iVirtualLast) = &
          & this%hamAlphaMo(iOpenFirst:iOpenLast, iVirtualFirst:iVirtualLast)

      ! Enforce symmetry explicitly.
      this%hamEffectiveMo(iVirtualFirst:iVirtualLast, iOpenFirst:iOpenLast) = &
          & transpose(this%hamEffectiveMo(iOpenFirst:iOpenLast, iVirtualFirst:iVirtualLast))
    end if

  end subroutine TRoksCalc_applyMoCouplings

  !> Return the maximum residual coupling between different orbital spaces.
  !>
  !> At a stationary high-spin restricted-open-shell solution, the
  !> beta core-open, spin-averaged core-virtual, and alpha open-virtual
  !> Hamiltonian blocks vanish in the common MO basis.
  function TRoksCalc_getStationarityResidual(this) result(residual)

    !> ROKS calculation data.
    class(TRoksCalc), intent(in) :: this

    !> Largest absolute independent orbital-gradient element.
    real(dp) :: residual

    integer :: iCoreFirst, iCoreLast
    integer :: iOpenFirst, iOpenLast
    integer :: iVirtualFirst, iVirtualLast
    integer :: nOrb

    nOrb = size(this%hamAlphaMo, dim=1)

    iCoreFirst = 1
    iCoreLast = this%Nc

    iOpenFirst = this%Nc + 1
    iOpenLast = this%Nc + this%No

    iVirtualFirst = this%Nc + this%No + 1
    iVirtualLast = nOrb

    residual = 0.0_dp

    ! Core-open rotations are governed by the beta Hamiltonian.
    if (this%Nc > 0 .and. this%No > 0) then
      residual = max(residual, maxval(abs(this%hamBetaMo(iCoreFirst:iCoreLast,&
         & iOpenFirst:iOpenLast))))
    end if

    ! Core-virtual rotations are governed by the spin average.
    if (this%Nc > 0 .and. this%Nv > 0) then
      residual = max(residual, maxval(abs(0.5_dp * (this%hamAlphaMo(iCoreFirst:iCoreLast, &
         & iVirtualFirst:iVirtualLast) + this%hamBetaMo(iCoreFirst:iCoreLast, &
         & iVirtualFirst:iVirtualLast)))))
    end if

    ! Open-virtual rotations are governed by the alpha Hamiltonian.
    if (this%No > 0 .and. this%Nv > 0) then
      residual = max(residual, maxval(abs(this%hamAlphaMo(iOpenFirst:iOpenLast,&
         & iVirtualFirst:iVirtualLast))))
    end if

  end function TRoksCalc_getStationarityResidual


  !> Form the common ROKS Hamiltonian and return it to the AO basis.
  subroutine TRoksCalc_buildEffectiveHamiltonian(this)

    class(TRoksCalc), intent(inout) :: this

    real(dp), allocatable :: tmp(:,:)

    ! Initialize all MO blocks with the alpha/beta average.
    this%hamEffectiveMo(:,:) = 0.5_dp * (this%hamAlphaMo(:,:) + this%hamBetaMo(:,:))

    call this%applyMoCouplings()

    ! H_AO = S C H_MO C^T S
    tmp = matmul(this%overlap, this%coefficients)

    this%hamEffective(:,:) = matmul(tmp, matmul(this%hamEffectiveMo, transpose(tmp)))

  end subroutine TRoksCalc_buildEffectiveHamiltonian
end module dftbp_roks_roks
