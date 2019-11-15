module subroutines
  implicit real*8(a-h,o-z)
  !
  ! coordination number of the hexanonal lattice
  !
  integer,parameter::icoordination = 6

  ! sketch of the lattice enumeration
  !
  !  00  01  02  03  04  
  !05  06  07  08  09
  !  10  11  12  13  14
  !15  16  17  18  19
  ! 

  integer,parameter::iempty        = 0
  integer,parameter::inonpolarized = 1
  integer,parameter::ipolarized_e1 = 2
  integer,parameter::ipolarized_e2 = 3
  integer,parameter::ipolarized_e3 = 4
  integer,parameter::ipolarized_e4 = 5
  integer,parameter::ipolarized_e5 = 6
  integer,parameter::ipolarized_e6 = 7

  integer,parameter::igapjunction    = 5
  integer,parameter::idepolarization = 6
  integer,parameter::ipolarization   = 7
  
contains
  function ichoice(ivector)
    !================================================
    ! returns a single sample chosen from ivector
    !================================================
    integer::ichoice
    integer,intent(in)::ivector(:)
    real*8 ::rng
    call random_number(rng)
    ichoice = ivector(int(rng*size(ivector))+1)
  end function ichoice
  function icheck(x,y)    
    !================================================
    ! returns 1 if y < x, 0 otherwise 
    !================================================
    integer::icheck
    icheck = (int(sign(1d0,x-y))+1)/2
  end function icheck
  subroutine read_args(fargs,ifilename)
    !================================================
    ! collect arguments from command line
    !================================================
    character(len=32)::arg
    real*8          ,intent(out)::fargs(7)    
    character(len=*),intent(out)::ifilename

    ! default values
    ! N = 10 , steps = 100 , samples =100 ,
    ! gapjunction = 0.5, depolarization = 0.05 , polarization = 0.5
    fargs = (/ 1d1, 1d1, 1d2 , 1d2 , 5d-1 , 5d-2, 5d-2/)

    ifilename =  "data"    
    
    ntotal = command_argument_count() 
    do i = 1,ntotal,2
       call get_command_argument(i,arg)
       select case(trim(arg))
       case ('-x','--Lx')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(1)
          ifilename = trim(ifilename)//'_x'//trim(arg)
       case ('-y','--Ly')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(2)
          ifilename = trim(ifilename)//'_y'//trim(arg)          
       case ('--steps','-s')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(3)
          ifilename = trim(ifilename)//'_steps'//trim(arg)
       case ('--samples','-mc')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(4)
          ifilename = trim(ifilename)//'_samples'//trim(arg)
       case ('--gap-junction','-g')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(5)
          ifilename = trim(ifilename)//'_gapjunction'//trim(arg)
       case ('--depolarization','-d')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(6)
          ifilename = trim(ifilename)//'_depolarization'//trim(arg)
       case ('--polarization','-p')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(7)
          ifilename = trim(ifilename)//'_polarization'//trim(arg)
       case('-h','--help')
          print *,'example usage'
          print *,'      main -x 10 -y 10 --samples 100 --steps 100 --gap-junction 0.5 --polarization 0.05 --depolarization 0.1'
          call exit(0)
       end select
    end do
    ifilename = trim(ifilename) !//'.dat'

    i = time()
    call ctime(i,arg)
    print *,"************************"
    print *,trim(arg)
    print *,"************************"
    print *," "
    print *,"Starting parameters for hexagonal lattice :: "
    print *,"Lx        = ",int(fargs(1))
    print *,"Ly        = ",int(fargs(2))
    print *,"isteps    = ",int(fargs(3))
    print *,"isamples  = ",int(fargs(4))
    print *,""
    print *,"gap-junction   rate = ",fargs(5)
    print *,"depolarization rate = ",fargs(6)
    print *,"polarization   rate = ",fargs(7)
    print *," "
    print *,"filename = ",ifilename
    print *,"************************"
  end subroutine read_args
  !================================================
  subroutine get_neighbours(k,Lx,Ly,ineigh)
    !================================================
    ! returns a vector with index of neighbours of 
    ! site k, in the cylindrical geometry with hexagonal
    ! lattice. 
    !================================================
    ! NOTE: periodic boundary conditions
    !
    ! NOTE: sketch of the lattice enumeration
    !
    !  00  01  02  03  04  
    !05  06  07  08  09
    !  10  11  12  13  14
    !15  16  17  18  19
    !
    integer,intent(in) ::k,Lx,Ly
    integer,intent(out)::ineigh(0:icoordination)
    ineigh(0) = k ! index-0 points to the site itself
    
    kx = mod(k,Lx)
    ky = k/Lx
    ! OBS : direction chart 
    !
    !     e3  e2
    !  e4        e1    --> x-direction
    !     e5  e6
    !      
    ! neighs in the same x-layer
    ineigh(1) = mod(kx+1 + Lx, Lx) + ky*Lx
    ineigh(4) = mod(kx-1 + Lx, Lx) + ky*Lx    
    
    ieven = mod(ky,2)
    
    kx1 = mod(kx+ieven  +Lx,Lx)
    ineigh(2) = kx1 + mod(ky+1+Ly,Ly)*Lx
    ineigh(6) = kx1 + mod(ky-1+Ly,Ly)*Lx

    kx1 = mod(kx+ieven-1+Lx,Lx)
    ineigh(3) = kx1 + mod(ky+1+Ly,Ly)*Lx
    ineigh(5) = kx1 + mod(ky-1+Ly,Ly)*Lx
    
  end subroutine get_neighbours
  !================================================
  subroutine select(k,idirection,Lx,Ly,inext,ishared)
    !================================================
    ! return the index inext, which is the neighbour
    ! of k in the idirection-th direction.
    !
    ! It also returns the ishared neighbours between
    ! k and k+e(idirection)
    !================================================
    integer,intent(in)   ::k,idirection,Lx,Ly
    integer,intent(inout)::inext,ishared(2)
    integer              ::ineigh(0:icoordination)
    call get_neighbours(k,Lx,Ly,ineigh)    
    inext = ineigh(idirection)    
    !
    ! idirection can be 0,1,...,icoordination    
    !
    ! shared table
    ! e1 - > e2, e6
    ! e2 - > e3, e1
    ! e3 - > e4, e2
    ! e4 - > e5, e3
    ! e5 - > e6, e4
    ! e6 - > e1, e5

    ! shift by 1 
    m = idirection - 1
    ! max added just to avoid negative index if
    ! idirection = 0
    m1= max(mod(m+1+icoordination,icoordination) + 1,0)
    m2= max(mod(m-1+icoordination,icoordination) + 1,0)    
    ! ishared = k if idirection = 0
    ishared(1) = ineigh(m1)*min(idirection,1)
    ishared(2) = ineigh(m2)*min(idirection,1)    
  end subroutine select
  !================================================
  subroutine boundary_conditions(Lx,Ly,lattice)
    !================================================
    ! set all lattice sites with ky = 0 to the
    ! direction e2;
    ! set empty sites with ky = 1 to direction e2 
    ! and removes cells with ky=Ly-1
    !================================================
    integer,intent(in)   ::Lx,Ly
    integer,intent(inout)::lattice(0:(Lx*Ly-1))
    integer              ::itmp(0:Lx-1) 
    lattice(0:(Lx-1)) = ipolarized_e2
    idx = 2*Lx - 1
    itmp = min(1,lattice(Lx:idx))
    ! if empty, add a cell in e2 direction
    lattice(Lx:idx) = lattice(Lx:idx)*itmp+(1-itmp)*ipolarized_e2
    ! last Lx vanishes
    lattice(Lx*(Ly-1):Lx*Ly-1)= iempty
  end subroutine boundary_conditions
  !================================================
  subroutine update_fixedtime(k,Lx,Ly,lattice,params,rng,prob,iflag)
    !================================================
    ! updates the status of lattice(k).
    !
    ! lattice(k) = 0, empty site
    !            = (1,2,3,4,5,6,7) otherwise
    !
    ! rng is a constant input, common to all tests
    ! until a configurational change is detected or
    ! exhausted (without changes)
    !
    ! prob is the likelihood of a configurational change
    ! it increases with each call of the subroutine
    !================================================
    integer,intent(in)   ::Lx,Ly,k
    integer,intent(inout)::lattice(0:(Lx*Ly-1)),iflag
    real*8 ,intent(in)   ::params(7),rng
    real*8 ,intent(inout)::prob
    integer              ::ishared(2)    
    iflag=0
    ! get the target cell current status
    icurrent = lattice(k)    
    ! determine whether the cell exists or not (binary)
    istatus  = min(1,icurrent)

    
    if (istatus.eq.iempty) return ! nothing to do here
    L = Lx*Ly    
    ! params hold single cell transition rates
    ! rates  hold the minimum time interval
    !        assuming L distinct cells
    dt = 1d0/L
    !! a better time interval is given by
    !
    ! dt = (1d0/L)*(sum(params(igapjunction:ipolarization)))
    !
    
    !
    ! try to polarize the cell if lattice(k) = inonpolarized
    !
    if (lattice(k).eq.inonpolarized) then
       prob = prob + params(ipolarization)*dt
       ! check if the test succeeded
       itmp = icheck(prob,rng) ! itmp = 1 if rng < prob
       !
       ! pick a new direction
       !
       inew = ichoice([(iz,iz=ipolarized_e1,ipolarized_e6)])
       lattice(k) = inew*itmp + (1-itmp)
       iflag = inonpolarized
       return
    end if
    !
    ! cell is polarized this point forward
    !
    
    !
    ! try to depolarize it
    !    
    prob   = prob+params(idepolarization)*dt
    idepol = icheck(prob,rng)
    if (idepol.eq.1) then
       lattice(k) = inonpolarized
       iflag      = idepolarization
       return
    end if
    !
    ! cell is polarized and failed to depolarized
    !
    ! NOTE: idirection = icurrent - 1
    call select(k,icurrent-1,Lx,Ly,inext,ishared)
    !
    ! check whether the target site is empty
    !
    istatus = istatus * (1-min(1,lattice(inext)))
    istatus_contact = min(sum(lattice(ishared)),1)    
    !if         istatus_contact = 1, use p
    !otherwise, istatus_contact = 0, use q
    p =      params(igapjunction) *dt
    q = (1d0-params(igapjunction))*dt
    prob = prob +istatus*(  istatus_contact)*p
    prob = prob +istatus*(1-istatus_contact)*q

    icell_move = icheck(prob,rng)    
    lattice(inext) = lattice(k)*icell_move+lattice(inext)*(1-icell_move)
    lattice(k)     = lattice(k)*(1-icell_move)
    iflag = igapjunction*icell_move
    return    
  end subroutine update_fixedtime
  !================================================
  subroutine sample_fixedtime(Lx,Ly,lattice,params,data)
    real*8 ,intent(in)   ::params(7)
    integer,intent(in)   ::Lx,Ly
    integer,intent(inout)::lattice(0:(Lx*Ly-1))
    real*8 ,intent(out)  ::data(:)
    integer,allocatable  ::iorder(:)        
    isteps = int(params(3))
    !initialization
    lattice = 0
    call boundary_conditions(Lx,Ly,lattice)
    do istep = 0,isteps-1
       !
       ! measurements NOTE: deferred datatype starts as 1-index
       !
       call measurements(Lx,Ly,lattice,data(istep+1:istep+1))
       !
       ! update one time interval
       !
       call random_number(rng)
       prob  = 0d0
       iflag = 0
       isite = 0       
       allocate(iorder,source=melanger_sansreplacement(Lx,Ly,lattice))
       isize=size(iorder)
       
       do while ((isite.LT.isize).and.(iflag.lt.1))
          isite = isite+1
          k     = iorder(isite)          
          call update_fixedtime(k,Lx,Ly,lattice,params,rng,prob,iflag)
       end do
       call boundary_conditions(Lx,Ly,lattice)
       deallocate(iorder)
    end do
    
  end subroutine sample_fixedtime
  !================================================
  function melanger_sansreplacement(Lx,Ly,lattice) result(itmp)
    integer,intent(in)   ::lattice(0:(Lx*Ly-1))
    real*8 ,allocatable::rng(:)
    integer,allocatable::itmp(:)

    itmp = [(k, k=Lx,Lx*(Ly-1)-1)]
    !
    ! NOTE:: skipping the source and the sink in
    !        the line above
    !
    itmp = pack(itmp,lattice(itmp)>0)
    L = size(itmp)
    allocate(rng(L))
    call random_number(rng)
    rng = rng*L
    do i=1,L
       iq = int(rng(i)) + 1 ! index shift
       ia = itmp(iq)
       itmp(iq) = itmp(i)
       itmp(i)  = ia
    end do
  end function melanger_sansreplacement
!================================================
  subroutine measurements(Lx,Ly,lattice,data)
    !================================================
    ! returns vector data extracted from lattice
    !================================================
    integer,intent(in) ::Lx,Ly,lattice(0:Lx*Ly-1)
    real*8 ,intent(out)::data(:)

    ! total number of particles in the lattice
    data(1) = get_number_particles(lattice)
    
  end subroutine measurements
  !================================================
  function get_number_particles(lattice) result(x)
    integer,intent(in)::lattice(:)
    integer::x
    x = sum(min(lattice,1))
  end function get_number_particles
  !================================================  
  subroutine pretty_printing(Lx,Ly,lattice)
    integer,intent(in)::Lx,Ly,lattice(0:(Lx*Ly-1))
    do ky =0,Ly-1
       if (mod(ky,2).eq.0) then
          print *,(int( lattice(k+Lx*ky) ,2),k=0,Lx-1)
       else
          print *,"   ",(int( lattice(k+Lx*ky) ,2),k=0,Lx-1)
       end if
    end do
  end subroutine pretty_printing
end module subroutines