program test
  use subroutines
  implicit real*8(a-h,o-z)
  integer,parameter::Lx=5,Ly=5
  integer::lattice(0:Lx*Ly-1),ihist(10)
  integer,allocatable::iorder(:)
  real*8::params(7)
  lattice = 0
  call boundary_conditions(Lx,Ly,lattice)
  call pretty_printing(Lx,Ly,lattice)
  params = (/ Lx*1d0, Ly*1d0, 1d0,1d0,1d0, 0d0, 1d0*Lx*Ly/)
  print *,"---------------------"
  do isample=1,5
     lattice = 0
     lattice(0) = 1       
     prob = 0d0
     rng = 0d0
     call update_fixedtime(0,Lx,Ly,lattice,params,rng,prob,iflag)
     print *,">>"
     print *,iflag,prob,rng
     call pretty_printing(Lx,Ly,lattice)
  end do

  
  print *,"---------------------"
  do isample=1,5
     iflag = 0
     params(7)=0d0
     lattice  = 0
     call boundary_conditions(Lx,Ly,lattice)
     prob = 0d0
     rng = 1d-8
     call update_fixedtime(Lx,Lx,Ly,lattice,params,rng,prob,iflag)
     print *,">>"
     print *,iflag,prob,rng
     call pretty_printing(Lx,Ly,lattice)
  end do

  print *,"---------------------"
  do isample=1,5
     iflag = 0
     params(7)=0d0
     params(6)=1d0*Lx*Ly
     params(5)=0d0
     lattice  = 0
     call boundary_conditions(Lx,Ly,lattice)
     prob = 0d0
     rng = 1d-8
     call update_fixedtime(Lx,Lx,Ly,lattice,params,rng,prob,iflag)
     print *,">>"
     print *,iflag,prob,rng
     call pretty_printing(Lx,Ly,lattice)
  end do
  
  
  print *,"---------------------"
  print *,"ichoice",[1,2,3,4,5,6]
  do i=1,int(1d4)
     k = ichoice([(i,i=1,10)])
     ihist(k) = ihist(k)+1
  end do
  print *,ihist/1d4

  print *, [ipolarized_e1:ipolarized_e6]

  print *,'--------------------------'
  print *,"melanger"
  lattice = 0
  call boundary_conditions(Lx,Ly,lattice)
  call pretty_printing(Lx,Ly,lattice)
  iorder =melanger_sansreplacement(Lx,Ly,lattice)
  !if (allocated(iorder)) deallocate(iorder)
  !iorder =melanger_sansreplacement(Lx,Ly,lattice)
  print *,int(iorder,1)
  print *,'size::',size(iorder)
  print *,'--------------------------'
  print *,"initial conditions"
  lattice = 0
  call boundary_conditions(Lx,Ly,lattice)
  call pretty_printing(Lx,Ly,lattice)
  params = (/ Lx*1d0, Ly*1d0, 1d0,1d0,1d0, 0d0, 0d0/)
  call sample_fixedtime(Lx,Ly,10,params)
end program test
