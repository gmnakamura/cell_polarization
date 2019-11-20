program main
  use omp_lib
  use subroutines
  implicit real*8(a-h,o-z)
  real*8       ::params(7)
  character*128::ifilename
  real*8 ,allocatable::data(:,:),datum(:,:)
  
  call read_args(params,ifilename,idata_skip)
  Lx = int(params(1))
  Ly = int(params(2))
  isteps   = int(params(3))
  isamples = int(params(4))
  L = Lx*Ly
  !
  ! Data storage can grow really fast with
  ! lattice size as the time interval shrinks
  ! with L. Therefore, it is convenient to
  ! only store data every other idata_skip steps
  !
  is_zero = min(1,idata_skip)
  idata_skip = (1-is_zero)+is_zero*(int(L/max(1,idata_skip))-1)

  
  !idata_skip = int(L/4) ! set it to 1 to storage all time steps
  !
  !
  !data initialization
  m = int(isteps/idata_skip) +1 !- 1

  idata_size_base = 3
  idata_size = idata_size_base + L
  allocate(data(0:m,idata_size),datum(0:m,idata_size))
  data = 0d0
  datum = 0d0
  !!end of initialization

  !$OMP PARALLEL DO private(datum,lattice) 
  do isample=1,isamples
     call sample_fixedtime(Lx,Ly,params,idata_skip,datum)
     !$OMP CRITICAL
     data = data + datum*1d0/isamples
     !$OMP END CRITICAL
  end do
  !$OMP END PARALLEL DO
  
  open(9, file=trim(ifilename)//"_density.dat")
  islices = int(isteps/idata_skip)
  do islice=0,islices - 1
     write(9 ,*) islice*idata_skip,data(islice,1:idata_size_base)    
  end do
  close(9)

  open(10, file=trim(ifilename)//"_lattice.dat")
  islices = int(isteps/idata_skip)
  do islice=0,islices - 1
     write(10 ,*) islice*idata_skip,data(islice,idata_size_base+1:idata_size)    
  end do
  close(10)

end program main

!
!CHECKS:
!
! ==================================================
! TEST #1::
!
! - parameters
! a) params(idepolarization) = params(ipolarization)*icoordination
! b) params(igapjunction)    = 0.5
!
! - expected results
! a) polarized particles should be equal to the number of
!    nonpolarized cells
! b) count(lattice == ipolarized) = count(lattice > ipolarized)
!

! ==================================================
! TEST #2::
!
! - parameters
! a) params(idepolarization) = 0
! b) params(ipolarization)   = 0
! c) params(igapjunction)    = 0.5
! d) Ly = 4 (for the boundary condition with 2 initial rows) 
!
! - expected results
! a) equilibrium density of particles in direction e2
! b) <lattice>/ipolarized_e2 = (1/2) 
! 

! ==================================================
! TEST #3::
!
! - parameters
! a) params(idepolarization) = 0
! b) params(ipolarization)   = 0
! c) params(igapjunction)    = 0.5
!
! - expected results
! a) equilibrium density of particles in direction e2
! b) <lattice>/ipolarized_e2 = (2/Ly) 
! 
