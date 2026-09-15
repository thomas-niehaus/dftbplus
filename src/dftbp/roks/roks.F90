!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include 'common.fypp'

!> Restricted open-shell Kohn-Sham data and routines
!>
!> Constructs a common set of spatial orbitals from conventional alpha-
!> and beta-spin DFTB Hamiltonians. The effective Hamiltonian uses the
!> Roothaan core/open/virtual block structure employed for high-spin
!> restricted open-shell calculations.
!>
!> The current implementation uses real-valued dense matrices and is
!> restricted to calculations at the Gamma point for periodic systems.
module dftbp_roks_roks

  use dftbp_io_message, only : error
  use dftbp_common_accuracy, only : dp, elecTolMax

  implicit none

  private
  public :: TRoksCalc

  !> Runtime data for a restricted open-shell Kohn-Sham calculation
  type :: TRoksCalc

    !> Number of doubly occupied core orbitals
    integer :: Nc = 0

    !> Number of singly occupied open-shell orbitals
    integer :: No = 0

    !> Number of unoccupied virtual orbitals
    integer :: Nv = 0

    !> Whether internal ROKS diagnostics should be printed
    logical :: writeDiagnostics = .false.

    !> Whether initial common orbitals are read from eigenvec.bin
    logical :: readEigenvectors = .false.

    !> Species-resolved Hubbard-like parameters containing only the Hartree kernel
    real(dp), allocatable :: hHubbard(:)

    !> Onsite kernel parameters used only for the virial integral
    real(dp), allocatable :: virialOnSiteElements(:,:,:,:)

    !> Maximum number of fixed-Hamiltonian orbital iterations
    integer :: maxIterations = 50

    !> Convergence tolerance for the orbital-stationarity residual
    real(dp) :: tolerance = 1.0e-8_dp

    !> Weight of the newly constructed effective Hamiltonian.
    real(dp) :: damping = 1.0_dp

    !> Alpha and beta occupations used to construct the effective Hamiltonian
    real(dp), allocatable :: occupations(:,:)

    !> Conventional alpha-spin Hamiltonian in the AO basis
    real(dp), allocatable :: hamAlpha(:,:)

    !> Conventional beta-spin Hamiltonian in the AO basis
    real(dp), allocatable :: hamBeta(:,:)

    !> Alpha-spin Hamiltonian in the shared MO basis
    real(dp), allocatable :: hamAlphaMo(:,:)

    !> Beta-spin Hamiltonian in the shared MO basis
    real(dp), allocatable :: hamBetaMo(:,:)

    !> Effective ROKS Hamiltonian
    real(dp), allocatable :: hamEffective(:,:)

    !> Original AO overlap matrix
    real(dp), allocatable :: overlap(:,:)

    !> Effective ROKS Hamiltonian in the MO basis
    real(dp), allocatable :: hamEffectiveMo(:,:)

    !> Shared ROKS molecular-orbital coefficients
    real(dp), allocatable :: coefficients(:,:)

    !> Eigenvalues of the effective ROKS Hamiltonian
    real(dp), allocatable :: eigenvalues(:)

    !> Becke virial two-electron integral in a.u.
    real(dp) :: virialIntegral = 0.0_dp

  contains

    !> Initialize the core, open-shell and virtual orbital counts
    procedure :: init => TRoksCalc_init

    !> Allocate dense ROKS Hamiltonian and orbital storage
    procedure :: allocateMatrices => TRoksCalc_allocateMatrices

    !> Build the spin-averaged Hamiltonian used to obtain trial common orbitals
    procedure :: buildInitialHamiltonian => TRoksCalc_buildInitialHamiltonian

    !> Transform the alpha and beta Hamiltonians to the current common MO basis
    procedure :: transformHamiltoniansToMo => TRoksCalc_transformHamiltoniansToMo

    !> Construct the effective ROKS Hamiltonian and transform it to the AO basis
    procedure :: buildEffectiveHamiltonian => TRoksCalc_buildEffectiveHamiltonian

    !> Insert occupation-weighted off-diagonal MO couplings
    procedure :: applyMoCouplings => TRoksCalc_applyMoCouplings

    !> Store occupations for the following ROKS inner iteration
    procedure :: setOccupations => TRoksCalc_setOccupations

    !> Return the largest independent ROKS orbital-gradient element
    procedure :: getStationarityResidual => TRoksCalc_getStationarityResidual
  end type TRoksCalc


contains

  !> Initialize the high-spin restricted open-shell orbital partition
  !>
  !> The beta population defines the number of doubly occupied core
  !> orbitals. The excess alpha population defines the number of singly
  !> occupied open-shell orbitals. All remaining orbitals are virtual.
  subroutine TRoksCalc_init(this, nEl, nOrb, maxIterations, tolerance, damping, writeDiagnostics,&
      & readEigenvectors)

    !> ROKS calculation data
    class(TRoksCalc), intent(out) :: this

    !> Number of electrons in each spin channel
    real(dp), intent(in) :: nEl(:)

    !> Number of available spatial orbitals
    integer, intent(in) :: nOrb

    !> Maximum number of fixed-Hamiltonian ROKS orbital iterations
    integer, intent(in) :: maxIterations

    !> Convergence tolerance for the orbital-stationarity residual
    real(dp), intent(in) :: tolerance

    !> Damping of the inner ROKS loop.
    real(dp), intent(in) :: damping

    !> Whether internal ROKS diagnostics should be printed
    logical, intent(in) :: writeDiagnostics

    !> Whether initial common orbitals are read from eigenvec.bin
    logical, intent(in) :: readEigenvectors

    this%maxIterations = maxIterations
    this%tolerance = tolerance
    this%writeDiagnostics = writeDiagnostics
    this%readEigenvectors = readEigenvectors
    this%damping = damping

    if (size(nEl) /= 2) then
      call error("ROKS requires two spin electron populations")
    end if

    if (any(abs(nEl - nint(nEl)) > elecTolMax)) then
      call error("ROKS requires integer spin electron populations")
    end if

    if (nEl(1) < nEl(2)) then
      call error("ROKS currently requires alpha to be the majority spin")
    end if

    ! Orbitals 1:Nc are doubly occupied
    this%Nc = nint(nEl(2))

    ! Orbitals Nc+1:Nc+No are singly occupied by alpha electrons
    this%No = nint(nEl(1) - nEl(2))

    ! Remaining orbitals are virtual
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

    allocate(this%occupations(nOrb, 2), source=0.0_dp)

    this%occupations(1:this%Nc + this%No, 1) = 1.0_dp
    this%occupations(1:this%Nc, 2) = 1.0_dp

  end subroutine TRoksCalc_init


  !> Allocate dense ROKS Hamiltonian storage
  subroutine TRoksCalc_allocateMatrices(this, nRows, nCols)

    !> ROKS calculation data.
    class(TRoksCalc), intent(inout) :: this

    !> Number of locally stored matrix rows
    integer, intent(in) :: nRows

    !> Number of locally stored matrix columns
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

  !> Store alpha and beta occupations for the next ROKS orbital optimization
  subroutine TRoksCalc_setOccupations(this, occupations)

    class(TRoksCalc), intent(inout) :: this
    real(dp), intent(in) :: occupations(:,:)

    if (size(occupations, dim=1) /= size(this%occupations, dim=1)) then
      call error("Invalid number of ROKS orbital occupations")
    end if

    if (size(occupations, dim=2) /= 2) then
      call error("ROKS requires alpha and beta orbital occupations")
    end if

    if (any(occupations < -elecTolMax) .or. &
        any(occupations > 1.0_dp + elecTolMax)) then
      call error("ROKS occupations must lie between zero and one")
    end if

    this%occupations(:,:) = occupations(:,:)

  end subroutine TRoksCalc_setOccupations

  !> Build an initial common-orbital Hamiltonian from the alpha and beta Hamiltonians
  subroutine TRoksCalc_buildInitialHamiltonian(this)

    !> ROKS calculation data
    class(TRoksCalc), intent(inout) :: this

    @:ASSERT(all(shape(this%hamAlpha) == shape(this%hamEffective)))

    this%hamEffective(:,:) = 0.5_dp * (this%hamAlpha(:,:) + this%hamBeta(:,:))

  end subroutine TRoksCalc_buildInitialHamiltonian

  !> Transform alpha and beta Hamiltonians from the AO basis to the shared orthonormal ROKS MO basis
  subroutine TRoksCalc_transformHamiltoniansToMo(this)

    class(TRoksCalc), intent(inout) :: this

    this%hamAlphaMo = matmul(transpose(this%coefficients), matmul(this%hamAlpha, this%coefficients))

    this%hamBetaMo = matmul(transpose(this%coefficients), matmul(this%hamBeta, this%coefficients))

  end subroutine TRoksCalc_transformHamiltoniansToMo

  !> Assemble occupation-weighted couplings between common spatial orbitals
  !>
  !> For orbitals p and q, the coupling is constructed from
  !>
  !>   (f_p_alpha - f_q_alpha) F_alpha_pq
  !> + (f_p_beta  - f_q_beta)  F_beta_pq.
  !>
  !> The scaling recovers the conventional core/open/virtual ROKS
  !> blocks for integer occupations without amplifying small fractional
  !> occupation differences.
  !
  !> A related common-orbital construction is used by the PySCF
  !> restricted-open-shell implementation; see pyscf.scf.rohf.
  subroutine TRoksCalc_applyMoCouplings(this)

    class(TRoksCalc), intent(inout) :: this

    integer :: p, q, nOrb
    real(dp) :: deltaAlpha, deltaBeta
    real(dp) :: occupationDifference
    real(dp) :: coupling
    real(dp), parameter :: occupationTolerance = 1.0e-10_dp

    nOrb = size(this%hamEffectiveMo, dim=1)

    @:ASSERT(size(this%hamEffectiveMo, dim=2) == nOrb)
    @:ASSERT(size(this%occupations, dim=1) == nOrb)
    @:ASSERT(size(this%occupations, dim=2) == 2)

    do q = 2, nOrb
      do p = 1, q - 1

        deltaAlpha = this%occupations(p,1) - this%occupations(q,1)
        deltaBeta = this%occupations(p,2) - this%occupations(q,2)

        occupationDifference = abs(deltaAlpha) + abs(deltaBeta)

        if (occupationDifference > occupationTolerance) then
          coupling = (deltaAlpha * this%hamAlphaMo(p,q) &
              & + deltaBeta * this%hamBetaMo(p,q)) &
              & / max(1.0_dp, occupationDifference)

          this%hamEffectiveMo(p,q) = coupling
          this%hamEffectiveMo(q,p) = coupling
        end if
      end do
    end do

  end subroutine TRoksCalc_applyMoCouplings

  !> Return the largest independent occupation-weighted coupling
  !>
  !> Pairs with identical alpha and beta occupations do not contribute,
  !> since rotations within an equally occupied subspace leave the
  !> density unchanged.
  function TRoksCalc_getStationarityResidual(this) result(residual)

    class(TRoksCalc), intent(in) :: this

    real(dp) :: residual
    real(dp) :: deltaAlpha, deltaBeta
    real(dp) :: occupationDifference
    real(dp) :: coupling
    real(dp), parameter :: occupationTolerance = 1.0e-10_dp
    integer :: p, q, nOrb

    nOrb = size(this%hamAlphaMo, dim=1)
    residual = 0.0_dp

    do q = 2, nOrb
      do p = 1, q - 1

        deltaAlpha = this%occupations(p,1) - this%occupations(q,1)
        deltaBeta = this%occupations(p,2) - this%occupations(q,2)

        occupationDifference = abs(deltaAlpha) + abs(deltaBeta)

        if (occupationDifference > occupationTolerance) then
          coupling = (deltaAlpha * this%hamAlphaMo(p,q) &
              & + deltaBeta * this%hamBetaMo(p,q)) &
              & / max(1.0_dp, occupationDifference)

          residual = max(residual, abs(coupling))
        end if
      end do
    end do

  end function TRoksCalc_getStationarityResidual

  !> Form the common ROKS Hamiltonian and return it to the AO basis
  subroutine TRoksCalc_buildEffectiveHamiltonian(this)

    class(TRoksCalc), intent(inout) :: this

    real(dp), allocatable :: tmp(:,:)

    ! Initialize all MO blocks with the alpha/beta average
    this%hamEffectiveMo(:,:) = 0.5_dp * (this%hamAlphaMo(:,:) + this%hamBetaMo(:,:))

    call this%applyMoCouplings()

    ! H_AO = S C H_MO C^T S
    tmp = matmul(this%overlap, this%coefficients)

    this%hamEffective(:,:) = matmul(tmp, matmul(this%hamEffectiveMo, transpose(tmp)))

  end subroutine TRoksCalc_buildEffectiveHamiltonian
end module dftbp_roks_roks
