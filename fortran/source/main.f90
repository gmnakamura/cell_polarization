program main
!  use omp_lib
  use subroutines
  implicit real*8(a-h,o-z)
  real*8       ::params(7)
  character*128::ifilename
  real*8 ,allocatable::data(:),datum(:)
  
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
  !data initializaion
  m = int(isteps/idata_skip) - 1 
  allocate(data(0:m),datum(0:m))
  data = 0d0
  datum = 0d0
  !!end of initialization
  do isample=1,isamples     
     call sample_fixedtime(Lx,Ly,params,idata_skip,datum)
     data = data + datum/isamples
  end do

  
  open(9, file=trim(ifilename)//"_density.dat")
  islices = int(isteps/idata_skip)
  do islice=0,islices - 1
     ! istep = islice*idata_skip
     write(9 ,*) islice*idata_skip ,data(islice)    
  end do
  close(9)

  ! open(10, file=trim(ifilename)//"_average.dat")
  ! open(11, file=trim(ifilename)//"_density.dat")
  ! do i=0,Ly-1     
  !    write(10,*) i,sum([( buffer(j+ i*Lx ),j=0,Lx-1  )])/real(Lx,8)    
  !    write(11,*) i,sum([( density(j+ i*Lx ),j=0,Lx-1 )])/real(Lx,8)
  ! end do
  ! close(10)
  ! close(11)

  ! open(12, file=trim(ifilename)//"_psi.dat")  
  ! do i=0,Ly-1
  !    ctmp = sum([( cpsi(j+ i*Lx ),j=0,Lx-1  )])/real(Lx,8)
  !    write(12,*) i,real(ctmp),imag(ctmp)
  ! end do
  ! close(12)  
end program main
