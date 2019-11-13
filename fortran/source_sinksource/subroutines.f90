module subroutines
  implicit real*8(a-h,o-z)

  ! coordination number of the hexanonal lattice
  integer,parameter::icoordination = 6

  ! sketch of the lattice enumeration
  !
  !  00  01  02  03  04  
  !05  06  07  08  09
  !  10  11  12  13  14
  !15  16  17  18  19
  !
  ! 
  
contains
  subroutine read_args(lx,ly,isteps,isamples,params,ifilename)
    !================================================
    ! collect arguments from command line
    !================================================
    character(len=32)::arg
    real*8 ::fargs(6)
    real*8 ,intent(out)::params(3)
    integer,intent(out)::isteps,isamples,lx,ly
    character(len=*),intent(out)::ifilename

    ! default values
    ! N = 10 , steps = 100 , samples =100 , p = 0.5, direction = 0.05
    fargs = (/ 1d1,1d1, 1d5 , 1d0 , 0.5d0 , 5d-2/)

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
       case ('--contact','-p')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(5)
          ifilename = trim(ifilename)//'_contact'//trim(arg)
       case ('--direction','-d')          
          call get_command_argument(i+1,arg)
          read(arg,*) fargs(6)
          ifilename = trim(ifilename)//'_direction'//trim(arg)
       case('-h','--help')
          print *,'example usage'
          print *,'      main -x 10 -y 10 --samples 100 --steps 100 --contact 0.5 --direction 0.05'
          call exit(0)
       end select
    end do
    ifilename = trim(ifilename) !//'.dat'


    Lx = int(fargs(1))
    Ly = int(fargs(2))
    isteps   = int(fargs(3))
    isamples = int(fargs(4))
    params(1)=      fargs(5) /(1d0*icoordination)
    params(2)= (1d0-fargs(5))/(1d0*icoordination)
    params(3)= fargs(6)      !/(1d0*icoordination)

    i = time()
    call ctime(i,arg)
    print *,"************************"
    print *,trim(arg)
    print *,"************************"
    print *," "
    print *,"Starting parameters for hexagonal lattice :: "
    print *,"Lx        = ",Lx
    print *,"Ly        = ",Ly
    print *,"isteps    = ",isteps
    print *,"isamples  = ",isamples
    print *,"contact   = ",params(1)*icoordination
    print *,"direction = ",params(3) !*icoordination
    print *," "
    print *,"filename  = ",ifilename
    print *,"************************"
  end subroutine read_args
  !================================================
  subroutine get_neighbours(k,Lx,Ly,ineigh)
    integer,intent(out)::ineigh(icoordination)
    !================================================
    ! returns a vector with index of neighbours of 
    ! site k, in the cylindrical geometry with hexagonal
    ! lattice. 
    !================================================
    ! OBS: indices can be larger than or equal to Lx*Ly
    !      at the open border. This is intended.
    !
    ! OBS: sketch of the lattice enumeration
    !
    !  00  01  02  03  04  
    !05  06  07  08  09
    !  10  11  12  13  14
    !15  16  17  18  19
    !
    !
    ! Intent :: lattice(ineigh) = status of neighbours
    kx = mod(k,Lx)
    ky = k/Lx
    ! OBS : direction chart 
    !
    !     e3  e2
    !  e4        e1
    !     e5  e6
    !      
    ! neighs in the same x-layer
    ineigh(1) = mod(kx+1 + Lx, Lx) + ky*Lx
    ineigh(4) = mod(kx-1 + Lx, Lx) + ky*Lx    
    ! effectively increase the lattice size
    ! on the y-direction by 1, so that open
    ! boundary conditions produce indices
    ! larger/equal Lx*Ly
    Ly1 = Ly+1

    ieven = mod(ky,2)
    
    kx1 = mod(kx+ieven  +Lx,Lx)
    ineigh(2) = kx1 + mod(ky+1+Ly1,Ly1)*Lx
    ineigh(6) = kx1 + mod(ky-1+Ly1,Ly1)*Lx

    kx1 = mod(kx+ieven-1+Lx,Lx)
    ineigh(3) = kx1 + mod(ky+1+Ly1,Ly1)*Lx
    ineigh(5) = kx1 + mod(ky-1+Ly1,Ly1)*Lx

    
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
    ! OBS : inext > Lx*Ly means it is not available
    !
    ! OBS : direction chart 
    !
    !     e3  e2
    !  e4        e1
    !     e5  e6
    !
    ! e1 = ( 1   ,     0     )
    ! e2 = ( 1/2 , sqrt(3)/2 )
    integer,intent(out)::inext,ishared(2)
    integer            ::ineigh(icoordination)
    call get_neighbours(k,Lx,Ly,ineigh)
    
    inext = ineigh(idirection)

    
    !
    ! idirection can be 1,2,3,4,5,6
    ! idirection = 0 implies the site is not occupied
    !
    !        shared
    ! e1 - > e2, e6
    ! e2 - > e3, e1
    ! e3 - > e4, e2
    ! e4 - > e5, e3
    ! e5 - > e6, e4
    ! e6 - > e1, e5
    m = idirection - 1
    m1= mod(m+1+icoordination,icoordination) + 1
    m2= mod(m-1+icoordination,icoordination) + 1
    ishared(1) = ineigh(m1)
    ishared(2) = ineigh(m2)       
  end subroutine select
  !================================================
  subroutine boundary_condition(Lx,Ly,lattice)
    !================================================
    ! set all lattice sites with ky = 0 to either
    ! direction 2 or 3 ; and removes cells in ky=Ly-1
    !================================================
    integer,intent(inout)::lattice(0: (Lx*Ly-1) )
    integer              ::itmp(0:Lx-1)
    real*8               ::rng(Lx)
    call random_number(rng)
    ! either 2 or 3
    ! lattice(0:(Lx-1)) = 2 + ( int(2*rng) )

    itmp = min(1,lattice(0:Lx-1))
    lattice(0:(Lx-1)) = lattice(0:Lx-1)*itmp + (2+(int(2*rng)))*(1 - itmp)
    
    ! last Lx vanishes
    lattice(Lx*(Ly-1):)= 0

    
  end subroutine boundary_condition
  !================================================
  subroutine update(k,Lx,Ly,lattice,rng,prob,params,iflag)
    !================================================
    ! updates the status of site k using stochastic
    ! rules.
    !
    ! k       :: site
    ! lattice :: lattice with all cell configurations
    ! params  :: params(1) = p , params(2) = change dir
    !
    ! lattice(i,j) = 0, no direction
    !              = (1,2,3,4,5,6) otherwise
    !================================================
    integer,intent(inout)::lattice(0: (Lx*Ly-1) )
    real*8 ,intent(in)   ::params(3)
    integer              ::ishared(2)

    L = Lx*Ly
    ! get the target cell current status
    icurrent = lattice(k)
    ! get if the cell exists or not (binary)
    istatus  = min(1,icurrent)
    !first, try to change the direction of the cell
    prob = prob + istatus*params(3) ! simplification


    if (rng.LT.prob) then
       call random_number(rdirection)
       ! succeed! exit and inform the new direction
       lattice(k) = 1 + int(rdirection*icoordination)       
       iflag = lattice(k)       
       return
       ! OBS : set -p 0.0 and -d 1.0
       !       check if directions are distributed
       !       according to the uniform distribution
    end if    
    ! failed! then tries to move in the direction
    ! IF istatus is not 0
    if (istatus.eq.1) then
       call select(k,icurrent,Lx,Ly,inext,ishared)
       if (inext<L) then
          where (ishared < L)
             ! we only want to determine whether
             ! cells exists neighbouring both
             ! k and k+e(icurrent)
             ishared = min(lattice(ishared),1)
          elsewhere
             ishared = 0
          end where
          ! sum ishared and returns 1 if at least one
          ! neighbour exists 
          istatus_contact = min(sum(ishared),1)       
          !
          ! check whether the target site is empty
          !
          istatus = istatus * (1-min(1,lattice(inext)))
          !if         istatus_contact = 0, use 1-p
          !otherwise, istatus_contact = 1, use p    
          prob = prob +istatus*params(1)*(  istatus_contact)
          prob = prob +istatus*params(2)*(1-istatus_contact)

          if (rng.LT.prob) then          
             lattice(inext) = lattice(k)       
             lattice(k)     = 0          
             iflag          = lattice(inext) + istatus_contact*100


             
             return 
          end if          
       end if
    end if

    !else, everything has failed
    iflag = 0
    return
    
  end subroutine update
  !================================================
  function melanger_sansreplacement(L) result(itmp)      
    real*8 ::rng(0:(L-1))
    integer::itmp(0:(L-1))
    call random_number(rng)
    
    itmp = [ (i, i=0,(L-1))]

    do i=0,(L-1)
       iq = int(rng(i)*L)
       ia = itmp(iq)
       itmp(iq) = itmp(i)
       itmp(i)  = ia
    end do
  end function melanger_sansreplacement 
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
