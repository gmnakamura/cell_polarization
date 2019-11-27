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
  integer,parameter::n_ =3
  integer,parameter::isteps_=4
  integer,parameter::isamples_=5
  integer,parameter::ihopping        = 6
  integer,parameter::igapjunction    = 7
  integer,parameter::idepolarization = 8
  integer,parameter::ipolarization   = 9
  
  integer,parameter::iparams_size = 9
  
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
    ! Lx,Ly,n = 11 , steps = 100 , samples =100 ,
    ! hopping        = 1.0, gapjunction  = 0.5,
    ! depolarization = 0.0, polarization = 0.0
    fargs = (/ 11d0, 11d0, 11d0, 1d2 , 1d2 , 1d0, 5d-1 , 1d-1, 1d-1/)
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
          idata_skip_factor = int(max(0,idata_skip_factor))
          if (idata_skip_factor.EQ.0) arg='0'
!          write(arg,'(I1)') int(idata_skip_factor,1)
          ifilename = trim(ifilename)//'_dataskip'//trim(arg)
       case('-h','--help')
          print *,'example usage'
          print *,'main -x 10 -y 10 --samples 100 --steps 100 --gap-junction 0.5 --polarization 0.05 --depolarization 0.1 --skip 0'
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
    print *,"Lx        = ",int(fargs(Lx_))
    print *,"Ly        = ",int(fargs(Ly_))
    print *,"N0        = ",int(fargs(n_))
    print *,"isteps    = ",int(fargs(isteps_))
    print *,"isamples  = ",int(fargs(isamples_))
    print *,""
    print *,"hopping        rate = ",fargs(ihopping)
    print *,"gap-junction   rate = ",fargs(igapjunction)
    print *,"depolarization rate = ",fargs(idepolarization)
    print *,"polarization   rate = ",fargs(ipolarization)
    print *," "
    print *,"filename = ",trim(ifilename)
    print *,"************************"
  end subroutine read_args
  !================================================
  function get_neighbour(idirection,k,Lx,Ly)
    !-------------------------------------------
    ! returns an icell vector with the coordinates
    ! of the neighbour of icell in the idirection
    !-------------------------------------------
    integer,intent(in)::idirection,k,Lx,Ly
    integer::itmp(0:icoordination)
    integer::get_neighbour

    kx = mod(k,Lx)
    ky = k/Lx
    itmp(0) = k
    itmp(1) = mod(kx+1+Lx,Lx) + ky*Lx
    itmp(2) = kx + mod(ky+1+Ly,Ly)*Lx
    itmp(3) = mod(kx-1+Lx,Lx) + ky*Lx
    itmp(4) = kx + mod(ky-1+Ly,Ly)*Lx 
    
    get_neighbour = itmp(idirection)
  end function get_neighbour
  !================================================
  subroutine update_fixedtime(k,n,Lx,Ly,icells,lattice,params,rng,prob,iflag)    
    !-------------------------------------------
    ! updates the k-th cell in icells during the
    ! time interval delta t, with probability prob
    !
    ! upon success, returns non-null iflag
    !------------------------------------------
    integer  ,intent(in)   ::k,n,Lx,Ly
    integer  ,intent(inout)::iflag,icells(0:1,Lx*Ly)
    integer*1,intent(inout)::lattice(0:Lx*Ly-1)
    real*8   ,intent(in)   ::params(iparams_size),rng
    real*8   ,intent(inout)::prob
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
       prob = prob +icoordination*dt*params(ipolarization)
       itmp = icheck(prob,rng) ! itmp =1 if rng < prob
       ! new direction chosen from uniform distribution
       ! [0,icoordination-1] + shift 
       icells(0,k) = (1-itmp)+ itmp*(ichoice(icoordination)+ipolarized_e1)
       iflag = ipolarization*itmp
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

    ! polarized_cell behavior
    icurrent= icells(0,k) - inonpolarized ! current direction
    iold    = icells(1,k)
    inew    = get_neighbour(icurrent,iold,Lx,Ly) 
    istatus = 1-lattice(inew) ! 0  if occupied
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
    icells(1,k)   = iold*(1-itmp)+itmp*inew
    lattice(iold) = 1-itmp
    lattice(inew) = lattice(inew)*(1-itmp) + itmp
    iflag = itmp*igapjunction    
  end subroutine update_fixedtime
  !================================================
  subroutine sample_fixedtime(n,Lx,Ly,params,idata_skip,data)    
    !-------------------------------------------
    !-------------------------------------------
    real*8 ,intent(in)   ::params(iparams_size)
    integer,intent(in)   ::idata_skip,n
    real*8 ,intent(inout)::data(:,:)
    integer  ::icells(0:1,Lx*Ly)
    integer*1::lattice(0:Lx*Ly-1)
    
    call init_cells(n,Lx,Ly,icells,lattice)
    idx = 1
    
    call measurements(n,Lx,Ly,icells,data(idx,:))
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
          call update_fixedtime(k,n,Lx,Ly,icells,lattice,&
                                params,rng,prob,iflag)
       end do
       
       if (mod(istep,idata_skip).eq.0) then
          idx = idx+1
          call measurements(n,Lx,Ly,icells,data(idx,:)) ! entry data
       end if
    end do
    
  end subroutine sample_fixedtime
  !================================================
  subroutine init_cells(n,Lx,Ly,icells,lattice)
    !-------------------------------------------
    ! initialize icells as follows
    ! x_1 0   x_2  0   x_3 ...
    !  0  x_5  0  x_6   0  ...
    ! ...          0   x_n ...
    !-------------------------------------------
    integer  ,intent(in)   ::n,Lx,Ly
    integer  ,intent(inout)::icells(0:1,Lx*Ly)
    integer*1,intent(inout)::lattice(0:Lx*Ly-1)
    integer::itmp(0:Lx*Ly-1)
    real*8 ::rng(0:Lx*Ly-1)

    lattice =0    
    icells = 0
    !--------------------------------------
    ! test to match the initial conditions of
    ! the version without lattice/hash support
    !--------------------------------------
    ! k0=Lx/2 + (Ly/2)*Lx
    ! do k=1,n      
    !    iposition = k0 + mod((k0/Ly),2)
    !    k0=k0+2
    !    icells(0,k) = ichoice(icoordination+1)+inonpolarized
    !    icells(1,k) = iposition
    !    lattice(iposition) = 1
    ! end do
    !--------------------------------------

    iremaining = 1
    L = Lx*Ly

    itmp = [(k,k=0,L-1)]
    call random_number(rng)


    ! first particle will be put in the center
    kx = Lx / 2
    ky = Ly / 2
    kcenter = kx + ky*Lx
    !iaux = itmp(0)
    itmp(0) = kcenter !itmp(kcenter)
    itmp(kcenter) = 0 !iaux

    ! only shuffle the others
    do i=1,L-1
       idx  = int(rng(i)*(L-1)) + 1 ! always skip index 0 
       iaux = itmp(i)
       itmp(i)  = itmp(idx)
       itmp(idx)= iaux
    end do


    
    do i=0,n-1
       icells(0,i+1) = ichoice(icoordination+1)+inonpolarized
       icells(1,i+1) = itmp(i)
       lattice(itmp(i)) = 1       
    end do

       
  end subroutine init_cells
  !================================================
  subroutine measurements(n,Lx,Ly,icells,data)
    !-------------------------------------------
    ! extract data from icells
    !
    ! OBS: the subroutine assumes the center of
    !      the lattice sits over (Lx/2,Ly/2)
    !-------------------------------------------
    integer,intent(in) ::icells(0:1,n),Lx,Ly
    real*8 ,intent(out)::data(:)

    Lx2= Lx/2
    Ly2= Ly/2

    ! measuremts to compare against results using
    ! the version without lattice
    !
    data(1) = sum(min(1,icells(0,:)))
    data(2) = sum(mod(icells(1,:),Lx) - Lx2)
    data(3) = sum((icells(1,:)/Lx)    - Ly2)
    data(4) =sum((mod(icells(1,:),Lx)-Lx2)**2+((icells(1,:)/Lx)-Ly2)**2)

    
  end subroutine measurements
  !================================================
  subroutine entropy(n,Lx,Ly,icells,lattice,data)
    !-------------------------------------------
    ! extracts the entropy of the configuration
    !-------------------------------------------
    integer  ,intent(in) ::icells(0:1,n),Lx,Ly
    integer*1,intent(in) ::lattice(0:Lx*Ly-1)
    real*8 ,intent(out)::data(:)
    integer,parameter::ibsize=8
    integer::ibyte(0:ibsize-1) ! consider adding to the module and ini
    
    
    ibyte = [(2**k,k=0,ibsize-1)]

    L = Lx*Ly
    islices = L - ibsize

    do k=0,islices
    end do
    
  end subroutine entropy
end module subroutines
