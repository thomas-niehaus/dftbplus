!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements the Langevin thermostat
!!
module dftbp_md_langevintherm
  use dftbp_common_accuracy, only : dp, minTemp
  use dftbp_io_message, only : error
  use dftbp_math_ranlux, only : getRandom, TRanlux
  use dftbp_md_mdcommon, only : MaxwellBoltzmann, rescaleTokT, restFrame, TMDCommon, boxMueller
  use dftbp_md_tempprofile, only : TTempProfile
  use dftbp_md_thermostat, only : TThermostat
  implicit none

  private
  public :: TLangevinThermInput
  public :: TLangevinTherm, TLangevinTherm_init
  
  !> Thermostat specific input data for the Langevin thermostat
  type :: TLangevinThermInput

    !> Coupling strength
    real(dp) :: gamma

  end type TLangevinThermInput


  !> Langevin thermostat
  type, extends(TThermostat) :: TLangevinTherm
    private

    !> Nr. of atoms
    integer :: nAtom

    !> Random number generator
    type(TRanlux), allocatable :: pRanlux

    !> Mass of the atoms
    real(dp), allocatable :: mass(:)

    !> Temperature generator
    type(TTempProfile), pointer :: pTempProfile

    !> Coupling strength
    real(dp) :: gamma

    !> MD framework
    type(TMDCommon) :: pMDFramework

    !> MD timestep
    real(dp) :: deltaT

    ! Regional data storage
    logical :: tRegioTherm
    integer :: startReg(2), endReg(2)
    real(dp) :: regKT(2), energyExchange(2)

  contains

    procedure :: getInitVelocities => TLangevinTherm_getInitVelocities
    procedure :: updateVelocities => TLangevinTherm_updateVelocities
    procedure :: writeState => TLangevinTherm_writeState

  end type TLangevinTherm

contains


  !> Creates an Langevin thermostat instance.
  subroutine TLangevinTherm_init(this, input, pRanlux, masses, tempProfile, pMDFramework, deltaT)

    !> Initialised instance on exit.
    type(TLangevinTherm), intent(out) :: this

    !> Thermostat specific input data
    type(TLangevinThermInput), intent(in) :: input

    !> Random generator
    type(TRanlux), allocatable, intent(inout) :: pRanlux

    !> Masses of the atoms.
    real(dp), intent(in) :: masses(:)

    !> Pointer to a temperature profile object.
    type(TTempProfile), pointer, intent(in) :: tempProfile

    !> Molecular dynamics general specifications
    type(TMDCommon), intent(in) :: pMDFramework

    !> MD time step
    real(dp), intent(in) :: deltaT

    integer :: unit, iStat, iReg, iAt
    real(dp) :: tempK
    ! Boltzmann constant in atomic units: approx 3.166811e-6 Hartree/K
    real(dp), parameter :: KB_AU = 3.16681153d-6

    call move_alloc(pRanlux, this%pRanlux)
    this%nAtom = size(masses)
    this%mass = masses
    this%pTempProfile => tempProfile
    this%gamma = input%gamma
    this%pMDFramework = pMDFramework
    this%deltaT = deltaT
    
    ! Read regional information once during initialization
    open(newunit=unit, file='thermoRange.dat', status='old', iostat=iStat) 
    if (iStat /= 0) then
      this%tRegioTherm = .false.
    else
      this%tRegioTherm = .true.
      this%energyExchange = 0.0_dp
      
      do iReg = 1, 2
        read(unit, *, iostat=iStat) this%regKT(iReg), this%startReg(iReg), this%endReg(iReg)
        if (iStat /= 0) call error("Error: Could not read line for region in thermoRange.dat")
        if (this%startReg(iReg) < 1 .or. this%endReg(iReg) > this%nAtom) then
          call error("Thermostat region indices out of physical bounds (1 to nAtoms).")
        end if
        if (this%startReg(iReg) > this%endReg(iReg)) then
          call error("startReg > endReg in thermoRange.dat.")
        end if
        ! Convert K to atomic units (Hartree)
        this%regKT(iReg) = this%regKT(iReg) * KB_AU
      end do
      if (.not. (this%endReg(1) < this%startReg(2) .or. this%endReg(2) < this%startReg(1))) then
        call error("Thermostat regions 1 and 2 overlap in thermoRange.dat.")
      end if
    end if
    close(unit)
    
  end subroutine TLangevinTherm_init


  !> Returns the initial velocities.
  subroutine TLangevinTherm_getInitVelocities(this, velocities)

    !> Instance
    class(TLangevinTherm), intent(inout) :: this

    !> Velocities on return.
    real(dp), intent(out) :: velocities(:,:)

    real(dp) :: kT
    integer :: ii, iReg

    @:ASSERT(all(shape(velocities) <= [3, this%nAtom]))

    call this%pTempProfile%getTemperature(kT)
    if (kT < minTemp) then
      call error("Langevin thermostat not supported at zero temperature")
    end if

    if(this%tRegioTherm) then
      kT = sum(this%regKT)/2.0_dp
    end if 
    
    do ii = 1, this%nAtom
      call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), kT, this%pRanlux)
    end do
    call restFrame(this%pMDFramework, velocities, this%mass)
    call rescaleTokT(this%pMDFramework, velocities, this%mass, kT)

  end subroutine TLangevinTherm_getInitVelocities


  !> Updates the provided velocities according the current temperature.
  subroutine TLangevinTherm_updateVelocities(this, velocities)

    !> Instance
    class(TLangevinTherm), intent(inout) :: this

    !> Updated velocities on exit.
    real(dp), intent(inout) :: velocities(:,:)

    real(dp), allocatable :: uniform(:), gaussian(:)
    real(dp) :: kT, xi
    integer :: ii, nDeg, iReg
    real(dp) :: v_sq_old, v_sq_new

    @:ASSERT(all(shape(velocities) <= [3, this%nAtom]))
    nDeg = 3 * this%nAtom
    nDeg = nDeg + mod(nDeg, 2)
    allocate(uniform(nDeg))
    allocate(gaussian(nDeg))
    call this%pTempProfile%getTemperature(kT)
    call getRandom(this%pRanlux, uniform)
    do ii = 1, nDeg, 2
      call boxMueller(gaussian(ii), gaussian(ii+1), uniform(ii), uniform(ii+1))
    end do
    
    if(this%tRegioTherm) then
      do iReg = 1, 2
        do ii = this%startReg(iReg), this%endReg(iReg)
          v_sq_old = sum(velocities(:, ii)**2)
          
          velocities(:, ii) = velocities(:, ii) * (1.0_dp - this%gamma * this%deltaT)
          xi = sqrt(2.0_dp * this%mass(ii) * this%gamma * this%regKT(iReg) * this%deltaT) / this%mass(ii)
          velocities(:, ii) = velocities(:, ii) + xi * gaussian((ii-1)*3+1:ii*3)
          
          v_sq_new = sum(velocities(:, ii)**2)
          this%energyExchange(iReg) = this%energyExchange(iReg) + &
              0.5_dp * this%mass(ii) * (v_sq_new - v_sq_old)
        end do
      end do
    else
      velocities(:, :) = velocities * (1.0_dp - this%gamma * this%deltaT)
      do ii = 1, this%nAtom
        xi = sqrt(2.0_dp * this%mass(ii) * this%gamma * kT * this%deltaT) / this%mass(ii)
        velocities(:, ii) = velocities(:, ii) + xi * gaussian((ii-1)*3+1:ii*3)
      end do
    end if 

    call restFrame(this%pMDFramework, velocities, this%mass)
 
  end subroutine TLangevinTherm_updateVelocities


  !> Writes internals of thermostat
  subroutine TLangevinTherm_writeState(this, fd)

    !> Instance
    class(TLangevinTherm), intent(in) :: this

    !> File unit to write thermostat state out to
    integer, intent(in) :: fd

    write(fd,"(a,1x,3e20.10,2x,3e20.10)") 'Acc. energy loss in reservoirs:', this%energyExchange

  end subroutine TLangevinTherm_writeState

end module dftbp_md_langevintherm
