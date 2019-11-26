module subroutines
  implicit real*8(a-h,o-z)
  !
  ! coordination number of the square lattice
  !
  integer,parameter::icoordination = 4
 
  integer,parameter::iempty        = 0
  integer,parameter::inonpolarized = 1
  integer,parameter::ipolarized_e1 = 2
  integer,parameter::ipolarized_e2 = 3
  integer,parameter::ipolarized_e3 = 4
  integer,parameter::ipolarized_e4 = 5

  integer,parameter::igapjunction    = 5
  integer,parameter::idepolarization = 6
  integer,parameter::ipolarization   = 7


  ! parameters related to measurements
  integer,parameter::idata_size_base = 4
contains
  !================================================
  function ichoice(isize,ishift)
    !================================================
    ! returns a single sample chosen from ivector
    !================================================
    integer::ichoice
    integer,intent(in)::isize,ishift
    real*8 ::rng
    call random_number(rng)
    ichoice = int(rng*isize)+ishift
  end function ichoice
  !================================================
  function icheck(x,y)    
    !================================================
    ! returns 1 if y < x, 0 otherwise 
    !================================================
    integer::icheck
    icheck = int( (int(sign(1d0,x-y))+1)/2 )
  end function icheck
  !================================================
  subroutine read_args(fargs,ifilename,idata_skip_factor)
    !================================================
    ! collect arguments from command line
    !================================================
    character(len=32)::arg
    real*8          ,intent(out)::fargs(7)
    integer         ,intent(out)::idata_skip_factor
    character(len=*),intent(out)::ifilename

    ! default values
    ! Lx,Ly = 11 , steps = 100 , samples =100 ,
    ! gapjunction = 0.5, depolarization = 0.1 , polarization = 0.1
    fargs = (/ 11d0, 11d0, 1d2 , 1d2 , 5d-1 , 1d-1, 1d-1/)
    ! default data skip
    idata_skip_factor = 4
    
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
          read(arg,*) fargs(igapjunction)
          ifilename = trim(ifilename)//'_gapjunction'//trim(arg)
       case ('--depolarization','-d')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(idepolarization)
          ifilename = trim(ifilename)//'_depolarization'//trim(arg)
       case ('--polarization','-p')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(ipolarization)
          ifilename = trim(ifilename)//'_polarization'//trim(arg)
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
    print *,"gap-junction   rate = ",fargs(igapjunction)
    print *,"depolarization rate = ",fargs(idepolarization)
    print *,"polarization   rate = ",fargs(ipolarization)
    print *," "
    print *,"filename = ",ifilename
    print *,"************************"
  end subroutine read_args
  !================================================
  subroutine get_neighbours(k,Lx,Ly,ineigh)
    !================================================
    ! returns a vector with index of neighbours of 
    ! site k in the square lattice. 
    !================================================
    ! NOTE: periodic boundary conditions
    !
    ! NOTE: sketch of the lattice enumeration
    !
    !  00  01  02  03  04  
    !  05  06  07  08  09
    !  10  11  12  13  14
    !  15  16  17  18  19
    !
    integer,intent(in) ::k,Lx,Ly
    integer,intent(out)::ineigh(0:icoordination)
    ineigh(0) = k ! index-0 points to the site itself
    
    kx = mod(k,Lx)
    ky = k/Lx
    ! OBS : direction chart 
    !
    !       e2
    !  e3        e1    --> x-direction
    !       e4
    !      
    ! neighs in the same x-layer
    ineigh(1) = mod(kx+1 + Lx, Lx) + ky*Lx
    ineigh(3) = mod(kx-1 + Lx, Lx) + ky*Lx    
    ineigh(2) = kx + mod(ky+1 +Ly,Ly)*Lx
    ineigh(4) = kx + mod(ky-1 +Ly,Ly)*Lx    
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
    ! irrelevant for this test !
    ishared(1) = 0
    ishared(2) = 0
  end subroutine select
  !================================================
  subroutine boundary_conditions(Lx,Ly,icells)
    !================================================
    ! set all lattice sites with ky = 0 to the
    ! direction e2;
    ! set empty sites with ky = 1 to direction e2 
    ! and removes cells with ky=Ly-1
    !================================================
    integer,intent(in)   ::Lx,Ly
    integer,intent(inout)::icells(-1:(Lx*Ly-1),0:1)
    kx = Lx/2 
    ky = Ly/2
    ! icells(-1,0)= 1
    ! icells(0,0) = kx+ky*Lx
    ! icells(0,1) = ipolarized_e2

    icells(-1,0) = Lx
    do k=0, Lx-1
       icells(k,0) = k
       icells(k,1) = ichoice(icoordination+1,1)
    end do        
  end subroutine boundary_conditions
  !================================================
  subroutine update_fixedtime(kcell,Lx,Ly,icells,params,rng,prob,iflag)
    !================================================
    ! updates the status of lattice(k).
    !
    ! lattice(k) = 0, empty site
    !            = (1,2,3,...) otherwise
    !
    ! rng is a constant input, common to all tests
    ! until a configurational change is detected or
    ! exhausted (without changes)
    !
    ! prob is the likelihood of a configurational change
    ! it increases with each call of the subroutine
    !================================================
    integer,intent(in)   ::Lx,Ly,kcell
    integer,intent(inout)::iflag,icells(-1:Lx*Ly-1,0:1)
    real*8 ,intent(in)   ::params(7),rng
    real*8 ,intent(inout)::prob
    integer              ::ishared(2),icandidate(icoordination)    
    iflag=0
    k = icells(kcell,0)      
    ! get the target cell current status
    icurrent = icells(kcell,1)    
    ! determine whether the cell exists or not (binary)
    istatus  = min(1,icurrent)    
    if (istatus.eq.iempty) return ! nothing to do here
    L = Lx*Ly    
    ! params hold single cell transition rates
    ! rates  hold the minimum time interval
    !        assuming L distinct cells
    dt = 1d0/L

    !NOTE
    !dt = 1d0/((1d0+icoordination*params(ipolarization)+params(idepolarization)))
    
    !! a better time interval is given by
    !
    ! dt = (1d0/L)*(sum(params(igapjunction:ipolarization)))
    !
    
    !
    ! try to polarize the cell if lattice(k) = inonpolarized
    !
    if (icells(kcell,1).eq.inonpolarized) then
       prob = prob + icoordination*params(ipolarization)*dt
       ! check if the test succeeded
       itmp = icheck(prob,rng) ! itmp = 1 if rng < prob       
       icells(kcell,1) = (1-itmp) + itmp*ichoice(icoordination,inonpolarized+1)
       iflag = ipolarization
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
       icells(kcell,1) = inonpolarized
       iflag           = idepolarization
       return
    end if
    !
    ! cell is polarized and failed to depolarized
    !
    ! NOTE: idirection = icurrent - 1
    !    call select(k,icurrent-1,Lx,Ly,inext,ishared)       
    inext = ichoice(icoordination,1)    
    kx = mod(k,Lx)
    ky = k/Lx
    icandidate(1) = mod(kx + 1+Lx,Lx)+ky*Lx
    icandidate(2) = kx+Lx*mod(ky+1+Ly,Ly)
    icandidate(3) = mod(kx - 1+Lx,Lx)+ky*Lx
    icandidate(4) = kx+Lx*mod(ky-1+Ly,Ly)
    inext = icandidate(inext)    
    !
    ! check whether the target site is empty
    !
    istatus = istatus * is_free(inext,Lx,Ly,icells)
    ! istatus_contact = min(sum(lattice(ishared)),1)    
    ! !if         istatus_contact = 1, use p
    ! !otherwise, istatus_contact = 0, use q
    ! p =      params(igapjunction) *dt   
    ! q = (1d0-params(igapjunction))*dt 
    
    ! prob = prob +istatus*(  istatus_contact)*p
    ! prob = prob +istatus*(1-istatus_contact)*q

    prob = prob +istatus*dt
 
    icell_move = icheck(prob,rng)    
    icells(kcell,0) = k*(1-icell_move) + inext*icell_move
    iflag = igapjunction*icell_move
    return    
  end subroutine update_fixedtime
  !================================================
  subroutine sample_fixedtime(Lx,Ly,params,idata_skip,data)
    real*8 ,intent(in)   ::params(7)
    integer,intent(in)   ::Lx,Ly,idata_skip
    real*8 ,intent(out)  ::data(:,:)
    integer              ::icells(-1:Lx*Ly-1,0:1)
    integer,allocatable  ::iorder(:)   
    isteps = int(params(3))

    icells = 0
    call boundary_conditions(Lx,Ly,icells)

    ! icells(-1,0)= 1
    ! icells(0,0) = (Lx/2)+(Ly/2)*Lx
    ! icells(0,1) = ichoice(icoordination,2)
    !
    ! measurements NOTE: deferred datatype starts as 1-index
    !

    idx = 1
    call measurements(Lx,Ly,icells,data(idx,:)) ! entry data

    do istep = 1,isteps              
       call random_number(rng)
       prob  = 0d0
       iflag = 0
       inext_cell = 0
       !
       ! update sides according to order iorder
       ! iorder : shuffled index of non-empty sites
       !
       !allocate(iorder,source=melanger_sansreplacement(Lx,Ly,lattice))
       isize=icells(-1,0)       
       do while ((inext_cell.LT.isize).and.(iflag.lt.1))      
          call update_fixedtime(inext_cell,Lx,Ly,icells,params,rng,prob,iflag)
          inext_cell = inext_cell+1
       end do
       if (mod(istep,idata_skip).eq.0) then
          idx = idx+1
          call measurements(Lx,Ly,icells,data(idx,:)) ! entry data
       end if
    end do

  end subroutine sample_fixedtime
  !================================================
  function is_free(k,Lx,Ly,icells)
    integer,intent(in)::k,icells(-1:Lx*Ly,0:1)
    integer::is_free
    n = icells(-1,0)
    is_free=1
    do m =0,n-1
       if (k.eq.icells(m,0)) then
          is_free =1
          return
       end if
    end do
  end function is_free
  !================================================
  subroutine measurements(Lx,Ly,icells,data)
    !================================================
    ! returns vector data extracted from lattice
    !================================================
    integer,intent(in) ::Lx,Ly,icells(-1:Lx*Ly-1,0:1)
    real*8 ,intent(out)::data(:)       
    data(1) = min(1,icells(0,1)) !count(lattice.gt.1)
    !data(3) = data(1)+data(2)
    k = icells(0,0)
    kx = mod(k,Lx) - (Lx/2)
    ky = k/Lx  - (Ly/2)    
    data(2) = kx 
    data(3) = ky
    data(4) = (kx**2+ky**2)    
  end subroutine measurements
end module subroutines
