!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations
!  Copyright (C) 2006 - 2025 DFTB+ developers group
!
!  See the LICENSE file for terms of usage and distribution.
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Dense diagonalization driver for restricted-open-shell calculations
module dftbp_dftbplus_roksdiag

  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_dftbplus_eigenvects, only : diagDenseMtx
  use dftbp_elecsolvers_elecsolvers, only : TElectronicSolver
  use dftbp_io_message, only : warning
  use dftbp_roks_roks, only : TRoksCalc
  use dftbp_type_densedescr, only : TDenseDescr

#:if WITH_SCALAPACK
  use dftbp_dftbplus_eigenvects, only : diagDenseMtxBlacs
#:endif

  implicit none

  private
  public :: diagDenseRoksHamiltonian

contains

!> Converge and diagonalize the dense effective ROKS Hamiltonian
!>
!> The spin-averaged Hamiltonian provides preliminary common orbitals.
!> At fixed alpha and beta Hamiltonians, the effective ROKS Hamiltonian
!> is then rebuilt and diagonalized until the independent orbital-space
!> couplings satisfy the stationarity tolerance.
subroutine diagDenseRoksHamiltonian(env, denseDesc, electronicSolver, iScc, roks,&
    & hamiltonian, overlap, errStatus)

  !> Execution environment
  type(TEnvironment), intent(inout) :: env

  !> Dense matrix distribution descriptor
  type(TDenseDescr), intent(in) :: denseDesc

  !> Electronic eigensolver.
  type(TElectronicSolver), intent(inout) :: electronicSolver

  !> Outer SCC iteration counter.
  integer, intent(in) :: iScc

  !> ROKS matrices, orbitals and diagnostics settings
  type(TRoksCalc), intent(inout) :: roks

  !> Dense Hamiltonian work array
  real(dp), intent(inout) :: hamiltonian(:,:)

  !> Dense overlap work array
  real(dp), intent(inout) :: overlap(:,:)

  !> Error status returned by the eigensolver
  type(TStatus), intent(out) :: errStatus

  real(dp), allocatable :: previousHamiltonian(:,:)
  real(dp) :: roksStationarityResidual
  integer :: iRoks
  logical :: roksConverged

  @:ASSERT(iScc > 0)

  if (iScc == 1) then
    if (roks%writeDiagnostics) then
      write(stdOut, "(A)") "--> ROKS: initializing orbitals from spin-averaged Hamiltonian"
    end if
    ! Obtain preliminary common orbitals. These orbitals define the
    ! core, open-shell and virtual subspaces in which the spin-dependent
    ! effective Hamiltonian is assembled.
    hamiltonian(:,:) = roks%hamEffective(:,:)
    overlap(:,:) = roks%overlap(:,:)

#:if WITH_SCALAPACK
    call diagDenseMtxBlacs(electronicSolver, 1, 'V', denseDesc%blacsOrbSqr, hamiltonian,&
        & overlap, roks%eigenvalues, roks%coefficients, errStatus)
    @:PROPAGATE_ERROR(errStatus)
#:else
    call diagDenseMtx(env, electronicSolver, 'V', hamiltonian, overlap, roks%eigenvalues, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    roks%coefficients(:,:) = hamiltonian(:,:)
#:endif
  else
    if (roks%writeDiagnostics) then
      write(stdOut, "(A)") "--> ROKS: reusing orbitals from previous SCC iteration"
    end if
  end if

  roksConverged = .false.
  roksStationarityResidual = huge(1.0_dp)
  allocate(previousHamiltonian, mold=roks%hamEffective)

  do iRoks = 1, roks%maxIterations

    ! Construct the effective Hamiltonian using the current common
    ! orbitals and the fixed alpha and beta Hamiltonians.
    call roks%transformHamiltoniansToMo()
    call roks%buildEffectiveHamiltonian()

    if (iRoks > 1) then
      roks%hamEffective(:,:) = roks%damping * roks%hamEffective(:,:) &
          & + (1.0_dp - roks%damping) * previousHamiltonian(:,:)
    end if
    previousHamiltonian(:,:) = roks%hamEffective(:,:)

    ! Diagonalize the current effective ROKS Hamiltonian
    hamiltonian(:,:) = roks%hamEffective(:,:)
    overlap(:,:) = roks%overlap(:,:)

#:if WITH_SCALAPACK
    call diagDenseMtxBlacs(electronicSolver, 1, 'V', denseDesc%blacsOrbSqr, &
        & hamiltonian, overlap, roks%eigenvalues, roks%coefficients, errStatus)
    @:PROPAGATE_ERROR(errStatus)
#:else
    call diagDenseMtx(env, electronicSolver, 'V', hamiltonian, overlap, roks%eigenvalues, errStatus)
    @:PROPAGATE_ERROR(errStatus)
      
    ! The serial solver returns eigenvectors in the Hamiltonian array
    roks%coefficients(:,:) = hamiltonian(:,:)
#:endif

    ! Evaluate the orbital gradient using the newly updated common
    ! orbitals, while keeping H_alpha and H_beta fixed.
    call roks%transformHamiltoniansToMo()
    roksStationarityResidual = roks%getStationarityResidual()

    if (roks%writeDiagnostics) then
      write(stdOut, "(A,I3,A,ES12.4)") "--> ROKS: inner iteration ", iRoks, &
          & ", residual ", roksStationarityResidual
    end if

    if (roksStationarityResidual < roks%tolerance) then
      roksConverged = .true.
      exit
    end if
      
  end do

  if (roks%writeDiagnostics .and. roksConverged) then
    write(stdOut, "(A,I0,A,ES12.4)") "--> ROKS: converged in ", min(iRoks, roks%maxIterations), &
       & " inner iterations; final residual ", roksStationarityResidual
  end if

  if (.not. roksConverged) then
    call warning("ROKS inner orbital iteration did not converge")
    if (roks%writeDiagnostics) then
      write(stdOut, "(A,I0,A,ES12.4)") "--> ROKS: inner orbital iteration did not converge in ", &
         & roks%maxIterations, " iterations; residual = ", roksStationarityResidual
    end if
  end if
  
  end subroutine diagDenseRoksHamiltonian

end module dftbp_dftbplus_roksdiag
