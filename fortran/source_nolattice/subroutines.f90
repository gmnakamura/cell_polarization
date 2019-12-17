module subroutines
  implicit real*8(a-h,o-z)
  !
  ! coordination number of the square lattice
  !
  integer,parameter::icoords       = 2
  integer,parameter::icoordination = 2*icoords
 
  integer,parameter::iempty        = 0
  integer,parameter::inonpolarized = 1
  integer,parameter::ipolarized_e1 = 2
  integer,parameter::ipolarized_e2 = 3
  integer,parameter::ipolarized_e3 = 4
  integer,parameter::ipolarized_e4 = 5

  
  integer,parameter::iparams_size = 9
  integer,parameter::Lx_=1
  integer,parameter::Ly_=2
  integer,parameter::n_ =3
  integer,parameter::isteps_  =4
  integer,parameter::isamples_=5
  integer,parameter::ihopping_        = 6
  integer,parameter::igapjunction_    = 7
  integer,parameter::idepolarization_ = 8
  integer,parameter::ipolarization_   = 9

  
  ! parameters related to measurements
  integer,parameter::idata_size_base = 4

  integer::iversor(0:icoords,0:icoordination)

  ! maximum number of particles
  integer,parameter::nmax=256
contains
  !================================================
  function ichoice(n)
    !-------------------------------------------
    !
    !-------------------------------------------
    integer,intent(in)::n
    integer           ::ichoice
    real*8::rng
    call random_number(rng)
    ichoice = int(rng*n)
  end function ichoice
  !================================================  
  function icheck(x,y)    
    !-------------------------------------------
    ! returns 1 if y < x, 0 otherwise
    !-------------------------------------------
    integer::icheck
    icheck = int( (int(sign(1d0,x-y))+1)/2 )
  end function icheck
  !================================================
    subroutine read_args(fargs,ifilename,idata_skip_factor)
    !================================================
    ! collect arguments from command line
    !================================================
    character(len=32)::arg
    real*8          ,intent(out)::fargs(iparams_size)
    integer         ,intent(out)::idata_skip_factor
    character(len=*),intent(out)::ifilename

    ! default values
    ! Lx,Ly = 11 , steps = 100 , samples =100 ,
    ! hopping        = 1.0, gapjunction  = 0.5,
    ! depolarization = 0.0, polarization = 0.0
    fargs = (/ 11d0, 11d0, 1d0, 1d2 , 1d2 , 1d0, 5d-1 , 1d-1, 1d-1/)
    ! default data skip
    idata_skip_factor = 0
    
    ifilename =  "data"    
    
    ntotal = command_argument_count() 
    do i = 1,ntotal,2
       call get_command_argument(i,arg)
       select case(trim(arg))
       case ('-x','--Lx')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(Lx_)
          ifilename = trim(ifilename)//'_x'//trim(arg)
       case ('-y','--Ly')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(Ly_)
          ifilename = trim(ifilename)//'_y'//trim(arg)
       case ('-n','--N0')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(n_)
          ifilename = trim(ifilename)//'_N'//trim(arg)
       case ('--steps','-s')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(isteps_)
          ifilename = trim(ifilename)//'_steps'//trim(arg)
       case ('--samples','-mc')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(isamples_)
          ifilename = trim(ifilename)//'_samples'//trim(arg)
       case ('--hopping','-o')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(ihopping_)
          ifilename = trim(ifilename)//'_hop'//trim(arg)
       case ('--gap-junction','-g')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(igapjunction_)
          ifilename = trim(ifilename)//'_gap'//trim(arg)
       case ('--depolarization','-d')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(idepolarization_)
          ifilename = trim(ifilename)//'_depol'//trim(arg)
       case ('--polarization','-p')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(ipolarization_)
          ifilename = trim(ifilename)//'_pol'//trim(arg)
       case ('--skip')          
          call get_command_argument(i+1,arg)
          read(arg,*) idata_skip_factor          
          ! if 0, then use all data
          idata_skip_factor = int(max(0,idata_skip_factor))
          if (idata_skip_factor.EQ.0) arg='0'
!          write(arg,'(I1)') int(idata_skip_factor,1)
          ifilename = trim(ifilename)//'_dataskip'//trim(arg)
       case('-h','--help')
          print *,'example usage'
          print *,'main -x 10 -y 10 --n 256 --samples 100 --steps 100 --gap-junction 0.5 --polarization 0.05 --depolarization 0.1 --skip 0'
          print *,' '
          print *,'OBS: --skip 0 or any negative number forces collection of data at each time interval'
          call exit(0)
       end select
    end do
    ifilename = trim(ifilename) !//'.dat'

    !NOTE:
    ! polarization is divided by the coordination number
    fargs(ipolarization_)=fargs(ipolarization_)/icoordination
          
    print *,"************************"
    print *," "
    print *,"Starting parameters for hexagonal lattice :: "
    print *,"Lx        = ",int(fargs(Lx_))
    print *,"Ly        = ",int(fargs(Ly_))
    print *,"N0        = ",int(fargs(n_))
    print *,"isteps    = ",int(fargs(isteps_))
    print *,"isamples  = ",int(fargs(isamples_))
    print *,""
    print *,"hopping        rate = ",fargs(ihopping_)
    print *,"gap-junction   rate = ",fargs(igapjunction_)
    print *,"depolarization rate = ",fargs(idepolarization_)
    print *,"polarization   rate = ",fargs(ipolarization_)
    print *," "
    print *,"filename = ",trim(ifilename)
    print *,"************************"
  end subroutine read_args
  !================================================
  subroutine init_shared()
    iversor(:,0) = (/ 0, 0, 0/) 
    iversor(:,ipolarized_e1-inonpolarized) = (/ 0,  1, 0 /)
    iversor(:,ipolarized_e2-inonpolarized) = (/ 0,  0, 1 /)
    iversor(:,ipolarized_e3-inonpolarized) = (/ 0, -1, 0 /)
    iversor(:,ipolarized_e4-inonpolarized) = (/ 0,  0,-1 /)
  end subroutine init_shared
  !================================================
  function get_neighbour(idirection,icell)
    !-------------------------------------------
    ! returns an icell vector with the coordinates
    ! of the neighbour of icell in the idirection
    !-------------------------------------------
    integer,intent(in)::icell(0:icoords)
    integer::get_neighbour(0:icoords)
    get_neighbour = iversor(:,idirection)
  end function get_neighbour
  !================================================
  subroutine update_fixedtime(k,n,icells,params,rng,prob,iflag)    
    !-------------------------------------------
    ! updates the k-th cell in icells during the
    ! time interval delta t, with probability prob
    !
    ! upon success, returns non-null iflag
    !------------------------------------------
    integer,intent(in)   ::k,n
    integer,intent(inout)::iflag,icells(0:icoords,nmax)
    real*8 ,intent(in)   ::params(iparams_size),rng
    real*8 ,intent(inout)::prob
    integer::iaux(0:icoords),itmp
    iflag = 0
    istatus = min(1,icells(0,k)) ! safety-check
    if (istatus.eq.iempty) return
    ! define the time interval HERE
    ! NOTE:: the time interval affects the probability
    !        that a transition occurs
    dt = 1d0 / (params(Lx_)*params(Ly_))

    dt = 1d0/n ! oversimplification 

    !------------------------------------------------
    ! randow_walker behavior (uncomment for testing)
    !
    ! rw = 1d0 ! or 0d0
    ! prob = prob + rw
    ! if (icheck(prob,rng).eq.1) then
    !    idir = ichoice(icoordination)+1
    !    icells(:,k) = icells(:,k)+iversor(:,idir)
    !    iflag = 100
    !    return
    ! end if
    !------------------------------------------------

    
    ! check whether the cell is polarized or not
    ! if non-polarized, try to polarize it in a
    ! random direction
    if (icells(0,k).eq.inonpolarized) then
       prob = prob +icoordination*dt*params(ipolarization_)
       itmp = icheck(prob,rng) ! itmp =1 if rng < prob
       ! new direction chosen from uniform distribution
       ! [0,icoordination-1] + shift 
       icells(0,k) = (1-itmp)+ itmp*(ichoice(icoordination)+ipolarized_e1)
       iflag = ipolarization_*itmp
       return
    end if
    !
    ! cell is polarized
    !

    ! try to depolarize it
    prob = prob + dt*params(idepolarization_)
    
    if (icheck(prob,rng).eq.1) then
       icells(0,k) = inonpolarized
       iflag = idepolarization_
       return
    end if

    ! NOTE:: a random-walker would try to move in all
    !        directions with equal chance (one test
    !        followed by a random choice in direction)
    !        Polarization prevents multiple directions

    ! polarized_cell behavior
    icurrent = icells(0,k) - inonpolarized
    iaux = icells(:,k) + iversor(:,icurrent)
    istatus = is_free(iaux,icells,n) ! TODO:: implement hash or a lattice
    !
    ! ignore the neighbouring problem for now (so it's not gapjunction)
    !
    p = params(igapjunction_)*params(ihopping_)
    q = (1d0 - params(igapjunction_))*params(ihopping_)
    !
    ! TODO::
    !
    is_shared = 0 ! ignore the neighbouring problem for now

    
    prob = prob + dt*p*istatus
    prob = prob + dt*q*istatus*is_shared
    
    itmp = icheck(prob,rng)
    icells(:,k) = icells(:,k)+itmp*iversor(:,icurrent)
    iflag = itmp*igapjunction_    
  end subroutine update_fixedtime
  !================================================
  subroutine sample_fixedtime(params,idata_skip,data)    
    !-------------------------------------------
    !-------------------------------------------
    real*8 ,intent(in)   ::params(iparams_size)
    integer,intent(in)   ::idata_skip
    real*8 ,intent(inout)::data(:,:)
    integer::icells(0:icoords,nmax)
    icells = 0
    n = int(params(n_))
    call init_cells(n,icells)
    idx = 1
    
    call measurements(n,icells,data(idx,:))
    isteps = int(params(isteps_))
    do istep = 1,isteps
       call random_number(rng)
       prob  = 0d0
       iflag = 0
       k  = 0
       n0 = n ! dirty fix: make sure particles created within
              ! updates are not updated 
       do while ((iflag.eq.0).and.(k < n0)) 
          k = k + 1 ! TODO:: shuffle indices to remove bias
          call update_fixedtime(k,n,icells,params,rng,prob,iflag)
       end do
       
       if (mod(istep,idata_skip).eq.0) then
          idx = idx+1
          call measurements(n,icells,data(idx,:)) ! entry data
       end if
    end do

  end subroutine sample_fixedtime
  !================================================
  subroutine  sample_fixedtime_config(params,idata_skip,filename)
    !-------------------------------------------
    !-------------------------------------------
    real*8 ,intent(in)   ::params(iparams_size)
    integer,intent(in)   ::idata_skip
    character(len=*),intent(in)::filename
    integer::icells(0:icoords,nmax)
    icells = 0
    n = int(params(n_))
    call init_cells(n,icells)
    idx = 1

    open(unit=10,file=trim(filename)//"_config.dat")
    
    isteps = int(params(isteps_))
    do istep = 1,isteps
       call random_number(rng)
       prob  = 0d0
       iflag = 0
       k  = 0
       n0 = n ! dirty fix: make sure particles created within
              ! updates are not updated 
       do while ((iflag.eq.0).and.(k < n0)) 
          k = k + 1 ! TODO:: shuffle indices to remove bias
          call update_fixedtime(k,n,icells,params,rng,prob,iflag)
       end do
       
       if (mod(istep,idata_skip).eq.0) then
          idx = idx+1
          write(10,*) istep,(icells(:,k),k=1,n)
       end if
    end do
    close(10)
  end subroutine sample_fixedtime_config
  !================================================
  function is_free(iaux,icells,n)
    !-------------------------------------------
    ! returns 1 if iaux is in icells
    !
    ! NOTE:: overall elapsed time per call sits around 700ns
    !        = chokepoint
    !-------------------------------------------
    integer,intent(in)::n,iaux(0:icoords),icells(0:icoords,n)
    integer::is_free
    real*8::terrible(n),beta
    beta = 1d10
    terrible = [( sum((icells(1:icoords,j)-iaux(1:icoords))**2),j=1,n)]
    terrible = (exp(-beta*terrible))
    ! terrible should be equal to 0 for most cases
    ! except when iaux is occupied in icells
    !
    ! sum(terrible) = 0 if free,
    ! sum(terrible) > 0 if not-free
    is_free = 1 - min(1,int(sum(terrible)))
  end function is_free
  !================================================
  subroutine init_cells(n,icells)
    !-------------------------------------------
    ! initialize icells
    !-------------------------------------------
    integer,intent(in)   ::n
    integer,intent(inout)::icells(0:icoords,n)    
    icells = 0
    k0 = 0
    L = int(sqrt(n*1d0))
    ispacing =  1 !2 !4
    ihalfspacing = max(ispacing/2,1)
    ishift = int(L/2)*ispacing 
    
    do k = 0,n-1
       kx = mod(k,L)
       ky = k/L
       k1 = k+1
       icells(0,k1) = ichoice(icoordination+1)+inonpolarized
       icells(1,k1) = kx*ispacing - ishift !+ mod(ky,2)*ihalfspacing
       icells(2,k1) = ky*ispacing - ishift
    end do
    
    ! do k =1,n
    !    icells(0,k) = ichoice(icoordination+1)+inonpolarized
    !    icells(1,k) = k0
    !    k0=k0+2
    ! end do
    
  end subroutine init_cells
  !================================================
  subroutine measurements(n,icells,data)
    !-------------------------------------------
    !-------------------------------------------
    integer,intent(in) ::icells(0:icoords,n)
    real*8 ,intent(out)::data(:)

    data(1) = sum(min(1,icells(0,:))) !icells(0,1) -1  
    data(2) = sum(icells(1,:))
    data(3) = sum(icells(2,:))
    data(4) = sum(icells(1,:)**2+ icells(2,:)**2)

  end subroutine measurements
end module subroutines
