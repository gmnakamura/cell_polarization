program main
  use omp_lib
  use subroutines
  implicit real*8(a-h,o-z)
  !-------------------
  !PARAMETERS
  !-------------------
  character*256::ifilename
  real*8       ::params(iparams_size)
  real*8 ,allocatable::data(:,:),datum(:,:)
  integer,allocatable::icells(:,:)

  idata_skip = 0
  call read_args(params,ifilename,idata_skip)
  Lx       = int(params(Lx_))
  Ly       = int(params(Ly_))
  N0       = int(params(n_))
  isteps   = int(params(isteps_))
  isamples = int(params(isamples_))
  is_notzero = min(1,idata_skip)
  idata_skip = is_notzero*idata_skip+(1-is_notzero)
  ! data initialization
  m = int(isteps/idata_skip) +1 !- 1
  idata_size = idata_size_base      
  allocate(data(0:m,idata_size),datum(0:m,idata_size))
  data = 0d0
  datum = 0d0
  ! end of initialization

  n = N0

  
  !$OMP PARALLEL DO private(datum,icells,lattice)
  do isample = 1,isamples
     datum = 0d0
     call sample_fixedtime(n,Lx,Ly,params,idata_skip,datum)
     !$OMP CRITICAL
     data = data + datum*1d0/isamples
     !$OMP END CRITICAL
  end do
  !$OMP END PARALLEL DO
  
  
  open(9, file=trim(ifilename)//".dat")
  islices = int(isteps/idata_skip)
  do islice=0,islices - 1
     write(9 ,*) islice*idata_skip,real(data(islice,1:idata_size),4)
  end do
  close(9)
  
  
end program main

!---------------------------------------
! TEST 1:: particle conservation
! result:: ...OK!
!
!---------------------------------------
! TEST 2:: diffusion single particle
!
! expl:: one expects an effective diffusion coeff
!       for a single particle,
!
!            D(eff)~ D*(Gamma/(Gamma+Lambda))
!
!       Gamma = pol. rate, Lambda = depol rate
!       a) compute the squared displacement of a
!          single particle with Gamma1 and Gamma2.
!       b) compute the linear fit of r^2 for each Gamma
!       c) let a_k be the angular coeff for each Gamma_k
!       d) check if a_2/a_1 ~ D_2/D_1
!       e) EX: Gamma_1/Gamma_2 = 0.1/1.0, Lambda = 0.1
!              a2/a1 ~ 1.83
!              D2/D1 ~ 1.8181818
! result:: ...OK
