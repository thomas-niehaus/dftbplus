!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!--------------------------------------------------------------------------------------------------!

!> Restricted open-shell Kohn-Sham data and routines.
module dftbp_roks_roks

  use dftbp_common_globalenv, only : stdOut
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

  contains

    !> Initialize the core, open-shell and virtual orbital counts.
    procedure :: init => TRoksCalc_init

    !> Allocate dense Hamiltonian storage.
    procedure :: allocateMatrices => TRoksCalc_allocateMatrices

    !> Build the initial common-orbital Hamiltonian.
    procedure :: buildInitialHamiltonian => TRoksCalc_buildInitialHamiltonian

    procedure :: transformHamiltoniansToMo => TRoksCalc_transformHamiltoniansToMo

    procedure :: buildEffectiveHamiltonian => TRoksCalc_buildEffectiveHamiltonian

    procedure :: applyMoCouplings => TRoksCalc_applyMoCouplings

  end type TRoksCalc


contains


  !> Initialize the ROKS orbital spaces from the spin populations.
  subroutine TRoksCalc_init(this, nEl, nOrb)

    !> ROKS calculation data.
    class(TRoksCalc), intent(out) :: this

    !> Number of electrons in each spin channel.
    real(dp), intent(in) :: nEl(:)

    !> Number of available spatial orbitals.
    integer, intent(in) :: nOrb

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

    if (allocated(this%hamAlpha) .or. &
        & allocated(this%hamBeta) .or. &
        & allocated(this%hamAlphaMo) .or. &
        & allocated(this%hamBetaMo) .or. &
        & allocated(this%hamEffective) .or. &
        & allocated(this%coefficients) .or. &
        & allocated(this%overlap) .or. &
        & allocated(this%hamEffectiveMo) .or. &
        & allocated(this%eigenvalues)) then
      call error("ROKS matrices have already been allocated")
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

    if (.not. allocated(this%hamAlpha)) then
      call error("ROKS alpha Hamiltonian has not been allocated")
    end if

    if (.not. allocated(this%hamBeta)) then
      call error("ROKS beta Hamiltonian has not been allocated")
    end if

    if (.not. allocated(this%hamEffective)) then
      call error("ROKS effective Hamiltonian has not been allocated")
    end if

    if (any(shape(this%hamAlpha) /= shape(this%hamBeta)) .or. &
        & any(shape(this%hamAlpha) /= shape(this%hamEffective))) then
      call error("Inconsistent ROKS Hamiltonian dimensions")
    end if

    this%hamEffective(:,:) = 0.5_dp * &
        & (this%hamAlpha(:,:) + this%hamBeta(:,:))

  end subroutine TRoksCalc_buildInitialHamiltonian
  
  !> Transform the alpha and beta Hamiltonians from the AO basis
  !> to the shared orthonormal ROKS MO basis.
  subroutine TRoksCalc_transformHamiltoniansToMo(this)

    class(TRoksCalc), intent(inout) :: this

    if (.not. allocated(this%coefficients)) then
      call error("ROKS coefficients have not been allocated")
    end if

    if (.not. allocated(this%hamAlpha)) then
      call error("ROKS alpha Hamiltonian has not been allocated")
    end if

    if (.not. allocated(this%hamBeta)) then
      call error("ROKS beta Hamiltonian has not been allocated")
    end if

    if (.not. allocated(this%hamAlphaMo) .or. &
        & .not. allocated(this%hamBetaMo)) then
      call error("ROKS MO-basis Hamiltonians have not been allocated")
    end if

    this%hamAlphaMo = matmul(transpose(this%coefficients),&
        & matmul(this%hamAlpha, this%coefficients))

    this%hamBetaMo = matmul(transpose(this%coefficients),&
        & matmul(this%hamBeta, this%coefficients))

  end subroutine TRoksCalc_transformHamiltoniansToMo

  !> Apply ROKS coupling rules in the core/open/virtual MO blocks.
  subroutine TRoksCalc_applyMoCouplings(this)

    class(TRoksCalc), intent(inout) :: this

    real(dp) :: coreOpenCorrection
    real(dp) :: openVirtualCorrection
    integer :: iCoreFirst, iCoreLast
    integer :: iOpenFirst, iOpenLast
    integer :: iVirtualFirst, iVirtualLast
    integer :: nOrb

    coreOpenCorrection = 0.0_dp
    openVirtualCorrection = 0.0_dp
    
    nOrb = size(this%hamEffectiveMo, dim=1)

    if (size(this%hamEffectiveMo, dim=2) /= nOrb) then
      call error("ROKS effective MO Hamiltonian must be square")
    end if

    if (this%Nc + this%No + this%Nv /= nOrb) then
      call error("ROKS orbital partition does not match Hamiltonian size")
    end if

    iCoreFirst = 1
    iCoreLast = this%Nc

    iOpenFirst = this%Nc + 1
    iOpenLast = this%Nc + this%No

    iVirtualFirst = this%Nc + this%No + 1
    iVirtualLast = nOrb

    ! Core-open block: use the beta-spin Hamiltonian.
    if (this%Nc > 0 .and. this%No > 0) then
      coreOpenCorrection = maxval(abs( &
        & this%hamBetaMo(iCoreFirst:iCoreLast, iOpenFirst:iOpenLast) - &
        & this%hamEffectiveMo(iCoreFirst:iCoreLast, iOpenFirst:iOpenLast)))
      
      this%hamEffectiveMo(iCoreFirst:iCoreLast, iOpenFirst:iOpenLast) = &
          & this%hamBetaMo(iCoreFirst:iCoreLast, iOpenFirst:iOpenLast)

      ! Enforce symmetry explicitly.
      this%hamEffectiveMo(iOpenFirst:iOpenLast, iCoreFirst:iCoreLast) = &
          & transpose(this%hamEffectiveMo( &
          & iCoreFirst:iCoreLast, iOpenFirst:iOpenLast))
    end if

    ! Open-virtual block: use the alpha-spin Hamiltonian.
    if (this%No > 0 .and. this%Nv > 0) then
      openVirtualCorrection = maxval(abs( &
          & this%hamAlphaMo(iOpenFirst:iOpenLast, iVirtualFirst:iVirtualLast) - &
          & this%hamEffectiveMo(iOpenFirst:iOpenLast, iVirtualFirst:iVirtualLast)))
      
      this%hamEffectiveMo(iOpenFirst:iOpenLast, &
          & iVirtualFirst:iVirtualLast) = &
          & this%hamAlphaMo(iOpenFirst:iOpenLast, &
          & iVirtualFirst:iVirtualLast)

      ! Enforce symmetry explicitly.
      this%hamEffectiveMo(iVirtualFirst:iVirtualLast, &
          & iOpenFirst:iOpenLast) = &
          & transpose(this%hamEffectiveMo( &
          & iOpenFirst:iOpenLast, iVirtualFirst:iVirtualLast))
    end if

    write(stdOut, "(A,ES12.4)") "ROKS core-open correction: ", coreOpenCorrection
    write(stdOut, "(A,ES12.4)") "ROKS open-virtual correction: ", openVirtualCorrection

  end subroutine TRoksCalc_applyMoCouplings
  
  subroutine TRoksCalc_buildEffectiveHamiltonian(this)

    class(TRoksCalc), intent(inout) :: this

    real(dp), allocatable :: tmp(:,:)

    ! Initialize all MO blocks with the alpha/beta average.
    this%hamEffectiveMo(:,:) = 0.5_dp * &
        & (this%hamAlphaMo(:,:) + this%hamBetaMo(:,:))

    call this%applyMoCouplings()

    write(stdOut, "(A,ES12.4)") "ROKS effective MO symmetry error: ", &
        & maxval(abs(this%hamEffectiveMo - transpose(this%hamEffectiveMo)))

    ! H_AO = S C H_MO C^T S
    tmp = matmul(this%overlap, this%coefficients)

    this%hamEffective(:,:) = matmul(tmp,&
        & matmul(this%hamEffectiveMo, transpose(tmp)))

  end subroutine TRoksCalc_buildEffectiveHamiltonian
end module dftbp_roks_roks
