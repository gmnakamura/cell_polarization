program main
  use omp_lib
  use subroutines
  implicit real*8(a-h,o-z)
  integer,parameter::istatistics = (2+icoordination)
  real*8       ::params(3)
  character*128::ifilename
  character*128 ::fmt
  integer   ,allocatable::lattice(:),isort(:)
  real*8    ,allocatable::buffer(:),stats(:,:),density(:)
  complex*16,allocatable::cPsi(:)
  complex*16::ctmp
  
  call read_args(lx,ly,isteps,isamples,params,ifilename)
  L = Lx*Ly
  !
  ! allocate lattice and buffer 
  !
  allocate(lattice(0:(L - 1) ))
  allocate(buffer(0:L-1))
  allocate(stats(istatistics,0:isteps-1))
  allocate(cPsi(0:L-1))
  allocate(density(0:L-1))
  buffer =0d0
  density = 0d0
  rinv   = 1d0/isamples
  cpsi = (0d0,0d0)
  rfac = acos(-1d0)/3d0
  
  call OMP_SET_DYNAMIC(.true.)
  !$OMP PARALLEL DO PRIVATE(lattice,istep,iflag,prob,rng,isort,k,itrial,total) shared(buffer,stats,cPsi)
  do isample=1,isamples          
     lattice = 0
     call boundary_condition(Lx,Ly,lattice)  
     do istep = 1,isteps
        total = sum(min(lattice,1))
        !
        ! bottleneck
        !
        !$OMP CRITICAL
        stats(1,istep-1) = stats(1,istep-1) + rinv*(total**1d0)
        stats(2,istep-1) = stats(2,istep-1) + rinv*(total**2d0)  
        !$OMP END CRITICAL           
        iflag = 0
        prob  = 0d0
        call random_number(rng)
        itrial= 0
        isort = melanger_sansreplacement(L)
        do while (iflag.eq.0)
           !
           ! choose one site at random
           !
           itrial = itrial + 1
           k      = isort(itrial)
           call update(k,Lx,Ly,lattice,rng,prob,params,iflag)        
           !
           ! exit condition:
           !                either iflag or itrial/L>1
           !
           iflag = iflag  -  400*(itrial/L)
        end do
        call boundary_condition(Lx,Ly,lattice)        
     end do
     !$OMP CRITICAL
     buffer = buffer + rinv*lattice
!     buffer2= buffer2+ rinv*lattice*lattice
     density= density+ rinv*min(1,lattice) !not necessary to compute the square
     cpsi   = cpsi   + rinv*min(1,lattice)*exp( (0d0,1d0)*rfac*(lattice -1) )
!     cpsi2  = cpsi2  + rinv*min(1,lattice)*exp( (0d0,2d0)*rfac*(lattice -1) )
     !$OMP END CRITICAL
  end do
  !$OMP END PARALLEL DO  



  open(9, file=trim(ifilename)//"_statistics.dat") 
  do istep=0,isteps-1
     write(9 ,*) istep,stats(:,istep)    
  end do
  close(9)

  open(10, file=trim(ifilename)//"_average.dat")
  open(11, file=trim(ifilename)//"_density.dat")
  do i=0,Ly-1     
     write(10,*) i,sum([( buffer(j+ i*Lx ),j=0,Lx-1  )])/real(Lx,8)    
     write(11,*) i,sum([( density(j+ i*Lx ),j=0,Lx-1 )])/real(Lx,8)
  end do
  close(10)
  close(11)

  open(12, file=trim(ifilename)//"_psi.dat")  
  do i=0,Ly-1
     ctmp = sum([( cpsi(j+ i*Lx ),j=0,Lx-1  )])/real(Lx,8)
     write(12,*) i,real(ctmp),imag(ctmp)
  end do
  close(12)

  
  
  
end program main
