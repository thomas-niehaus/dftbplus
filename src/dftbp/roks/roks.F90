!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!--------------------------------------------------------------------------------------------------!

!> Restricted open-shell Kohn-Sham data and routines.
module dftbp_roks_roks

  use dftbp_common_accuracy, only : dp, elecTolMax
  use dftbp_io_message, only : error

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

    !> Effective ROKS Hamiltonian.
    real(dp), allocatable :: hamEffective(:,:)

  contains

    !> Initialize the core, open-shell and virtual orbital counts.
    procedure :: init => TRoksCalc_init

    !> Allocate dense Hamiltonian storage.
    procedure :: allocateMatrices => TRoksCalc_allocateMatrices

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
        & allocated(this%hamEffective)) then
      call error("ROKS Hamiltonian matrices have already been allocated")
    end if

    allocate(this%hamAlpha(nRows, nCols), source=0.0_dp)
    allocate(this%hamBeta(nRows, nCols), source=0.0_dp)
    allocate(this%hamEffective(nRows, nCols), source=0.0_dp)

  end subroutine TRoksCalc_allocateMatrices


end module dftbp_roks_roks
