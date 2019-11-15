program main
!  use omp_lib
  use subroutines
  implicit real*8(a-h,o-z)
  real*8       ::params(7)
  character*128::ifilename
  integer,allocatable::lattice(:)
  real*8 ,allocatable::data(:),datum(:)
  integer::ineigh(0:6)
  
  call read_args(params,ifilename)
  Lx = int(params(1))
  Ly = int(params(2))
  isteps   = int(params(3))
  isamples = int(params(4))
  L = Lx*Ly
  !initializaion
  allocate(lattice(0:L-1))
  lattice = 0
  allocate(data(0:isteps-1),datum(0:isteps-1))
  data = 0d0
  datum = 0d0
  !!end of initialization
  do isample=1,isamples     
     call sample_fixedtime(Lx,Ly,lattice,params,datum)
     data = data + datum/isamples
  end do
  
  call pretty_printing(Lx,Ly,lattice)

  
  open(9, file=trim(ifilename)//"_density.dat") 
  do istep=0,isteps-1
     write(9 ,*) istep,data(istep)    
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
