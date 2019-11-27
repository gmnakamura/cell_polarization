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

  integer,parameter::Lx_=1
  integer,parameter::Ly_=2
  integer,parameter::ihopping        = 5
  integer,parameter::igapjunction    = 6
  integer,parameter::idepolarization = 7
  integer,parameter::ipolarization   = 8
  
  integer,parameter::iparams_size = 8
  
  ! parameters related to measurements
  integer,parameter::idata_size_base = 4

  integer::iversor(0:icoords,0:icoordination)
  
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
    fargs = (/ 11d0, 11d0, 1d2 , 1d2 , 1d0, 5d-1 , 1d-1, 1d-1/)
    ! default data skip
    idata_skip_factor = 4
    
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
       case ('--steps','-s')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(3)
          ifilename = trim(ifilename)//'_steps'//trim(arg)
       case ('--samples','-mc')
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(4)
          ifilename = trim(ifilename)//'_samples'//trim(arg)
       case ('--hopping','-o')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(ihopping)
          ifilename = trim(ifilename)//'_hop'//trim(arg)
       case ('--gap-junction','-g')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(igapjunction)
          ifilename = trim(ifilename)//'_gap'//trim(arg)
       case ('--depolarization','-d')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(idepolarization)
          ifilename = trim(ifilename)//'_depol'//trim(arg)
       case ('--polarization','-p')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(ipolarization)
          ifilename = trim(ifilename)//'_pol'//trim(arg)
       case ('--skip')          
          call get_command_argument(i+1,arg)
          read(arg,*) idata_skip_factor          
          ! if 0, then use all data
          idata_skip_factor = max(0,idata_skip_factor)
          if (idata_skip_factor.EQ.0) arg='0'
!          write(arg,'(I1)') int(idata_skip_factor,1)
          ifilename = trim(ifilename)//'_dataskip'//trim(arg)
       case('-h','--help')
          print *,'example usage'
          print *,'main -x 10 -y 10 --samples 100 --steps 100 --gap-junction 0.5 --polarization 0.05 --depolarization 0.1 --skip 4'
          print *,' '
          print *,'OBS: --skip 0 or any negative number forces collection of data at each time interval'
          call exit(0)
       end select
    end do
    ifilename = trim(ifilename) !//'.dat'

    !NOTE:
    ! polarization is divided by the coordination number
    fargs(ipolarization)=fargs(ipolarization)/icoordination
          
    print *,"************************"
    print *," "
    print *,"Starting parameters for hexagonal lattice :: "
    print *,"Lx        = ",int(fargs(1))
    print *,"Ly        = ",int(fargs(2))
    print *,"isteps    = ",int(fargs(3))
    print *,"isamples  = ",int(fargs(4))
    print *,""
    print *,"hopping        rate = ",fargs(ihopping)
    print *,"gap-junction   rate = ",fargs(igapjunction)
    print *,"depolarization rate = ",fargs(idepolarization)
    print *,"polarization   rate = ",fargs(ipolarization)
    print *," "
    print *,"filename = ",ifilename
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
  subroutine update_fixedtime(k,icells,n,params,rng,prob,iflag)    
    !-------------------------------------------
    ! updates the k-th cell in icells during the
    ! time interval delta t, with probability prob
    !
    ! upon success, returns non-null iflag
    !-------------------------------------------
    integer,intent(in)   ::k,n
    integer,intent(inout)::iflag,icells(0:icoords,n)
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
    
    ! check whether the cell is polarized or not
    ! if non-polarized, try to polarize it in a
    ! random direction
    if (icells(0,k).eq.inonpolarized) then
       prob = prob +icoordination*dt*params(ipolarization)
       itmp = icheck(prob,rng) ! itmp =1 if rng < prob
       ! new direction chosen from uniform distribution
       ! [0,icoordination-1] + shift 
       icells(0,k) = (1-itmp)+ itmp*(ichoice(icoordination)+ipolarized_e1)
       iflag = ipolarization
       return
    end if
    !
    ! cell is polarized
    !

    ! try to depolarize it
    prob = prob + dt*params(idepolarization)
    if (icheck(prob,rng).eq.1) then
       icells(0,k) = inonpolarized
       iflag = idepolarization
       return
    end if

    ! NOTE:: a random-walker would try to move in all
    !        directions with equal chance (one test
    !        followed by a random choice in direction)
    !        Polarization prevents multiple directions

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
    ! polarized_cell behavior
    icurrent = icells(0,k) - inonpolarized
    iaux = icells(:,k) + iversor(:,icurrent)
    istatus = is_free(iaux,icells,n) ! TODO:: implement hash or a lattice
    !
    ! ignore the neighbouring problem for now (so it's not gapjunction)
    !
    p = params(igapjunction)*params(ihopping)
    q = (1d0 - params(igapjunction))*params(ihopping)
    !
    ! TODO::
    !
    is_shared = 0 ! ignore the neighbouring problem for now

    
    prob = prob + dt*p*istatus
    prob = prob + dt*q*istatus*is_shared
    
    itmp = icheck(prob,rng)
    icells(:,k) = icells(:,k)+itmp*iversor(:,icurrent)
    iflag = itmp*igapjunction    
  end subroutine update_fixedtime
  !================================================
  subroutine sample_fixedtime(params,idata_skip,data)    
    !-------------------------------------------
    !-------------------------------------------
    integer,parameter::nmax=10
    real*8 ,intent(in)   ::params(iparams_size)
    integer,intent(in)   ::idata_skip
    real*8 ,intent(inout)::data(:,:)
    integer::icells(0:icoords,nmax)
    icells = 0
    ! initial condition
    n = 1
    icells(0,1:n) = 2
!    call init_cells(icells,n)
    idx = 1
    call measurements(n,icells,data(idx,:))
    isteps = int(params(3))
    do istep = 1,isteps
       call random_number(rng)
       prob = 0d0
       iflag= 0
       k = 0
       n0 = n
       do while ((iflag.eq.0).and.(k < n0)) 
          k = k + 1
          !print *,'entrou',k,istep
          call update_fixedtime(k,icells,n,params,rng,prob,iflag)
!          print *,'saiu',k,istep,iflag,int(icells(1:2,1),1)
       end do
       if (mod(istep,idata_skip).eq.0) then
          idx = idx+1
          !print *,'entrou 2',istep
          call measurements(n,icells,data(idx,:)) ! entry data
          !print *,'saiu 2',istep
       end if
    end do 
  end subroutine sample_fixedtime
  !================================================
  function is_free(iaux,icells,n)
    !-------------------------------------------
    ! returns 1 if iaux is in icells
    !-------------------------------------------
    integer,intent(in)::n,iaux(0:icoords),icells(0:icoords,n)
    integer::is_free
    real*8::terrible(n),beta
    beta = 1d6
    terrible = [( sum((icells(1:icoords,j)-iaux(1:icoords))**2),j=1,n)]
    terrible = exp(-beta*terrible)
    ! terrible should be equal to 0 for most cases
    ! except when iaux is occupied in icells
    !
    ! sum(terrible) = 0 if free,
    ! sum(terrible) > 0 if not-free
    is_free = 1 - min(1,int(sum(terrible)))
  end function is_free
  !================================================
  subroutine measurements(n,icells,data)
    !-------------------------------------------
    !-------------------------------------------
    integer,intent(in) ::icells(0:icoords,n)
    real*8 ,intent(out)::data(:)

    data(1) = sum(min(1,icells(0,:)))
    data(2) = icells(1,1)
    data(3) = icells(2,1)
    data(4) = icells(1,1)**2+ icells(2,1)**2
  end subroutine measurements
end module subroutines