import numpy as np
#import numba

def get_particles(lattice):
    """
    return the number of particles, regarless of
    its state.

    occupied state == lattice[k] > 0
    empty    state == lattice[k] = 0
    """
    return np.sum(lattice > 0)

def get_state(iQ,Lx,Ly,lattice):
    """
    return a vector with the number of particles 
    in the state iQ, along the x-direction
    """
    return np.array([np.sum([lattice[i+j*Lx]==iQ for j in range(Ly)]) for i in range(Lx)],dtype=int)

def get_states(iQ,Lx,Ly,lattice):
    """
    return a vector with the number of particles
    along the x-direction, with total iQ number
    of states.
    """
    x = np.zeros(Lx,dtype=int)
    for i in range(1,iQ+1): # skip 0-state 
        x += get_state(i,Lx,Ly,lattice)
    return x

def get_neighbours_square(k,Lx,Ly):
    """
    return the neighbours of site k in the
    square lattice
    """
    ky,kx = np.divmod(k,Lx)
    return [ (kx+1+Lx)%Lx+ky*Lx,
             (kx-1+Lx)%Lx+ky*Lx,
             kx+((ky+1+Ly)%Ly)*Lx ,
             kx+((ky-1+Ly)%Ly)*Lx]

def get_neighbours(k,Lx,Ly,iopen = 1):
    """
    return the neighbours of site k in the
    hexagonal lattice, with open boundary
    conditions on the y-direction for odd iopen

    00  01  02  03  04
      05  06  07  08  09
    10  11  12  13  14
      15  16  17  18  19

    directions (e1 = x-versor)

        y
        |
      e3 e2
    e4     e1 --x
      e5 e6

    """
    ky,kx = np.divmod(k,Lx)
    ieven = ky%2

    Ly1 = Ly + (iopen%2) 

    return np.array([(kx+1      +Lx)%Lx+ky*Lx,
            (kx+ieven  +Lx)%Lx+((ky+1+Ly1)%Ly1)*Lx,
            (kx+ieven-1+Lx)%Lx+((ky+1+Ly1)%Ly1)*Lx,
            (kx-1+Lx)%Lx +ky*Lx,
            (kx+ieven-1+Lx)%Lx+((ky-1+Ly1)%Ly1)*Lx,
            (kx+ieven  +Lx)%Lx+((ky-1+Ly1)%Ly1)*Lx ],dtype=int)

def get_neighbours_shared(k1,k2,Lx,Ly,iopen=1):
    """
    return the shared neighbours between sites
    k1 and k2.
    """
    tmp = np.intersect1d(get_neighbours(k1,Lx,Ly,iopen),
                         get_neighbours(k2,Lx,Ly,iopen))
    return tmp[ tmp < Lx*Ly ] # only valid neighbours

def select_next(k,idirection,Lx,Ly):
    """
    return the site the particle is supposed to move
    next

    idirection = 0,1,2,3,4,5 
    corresponds to e1,e2,...,e6 in the hexagonal lattice
    """
    return get_neighbours(k,Lx,Ly)[idirection]


def boundary_conditions(Lx,Ly,lattice):
    """
    set all sites with ky = 0 to direction e2
    and remove all cells with ky = Ly-1

    note:: e2 = status 3 
    """
    #identify all occupied sites. Add new cells 
    #at the non-occupied ones (source)
    itmp = (lattice[:Lx] > 0).astype(int)
    lattice[:Lx] = lattice[:Lx] *itmp + 3*(1 - itmp)
    #remove cells (sink)
    lattice[-Lx:] = 0
    return None

def update_fixedtime(k,Lx,Ly,lattice,rng,prob,param):
    """
    update the k-th site

    0 : not occupied
    1 : de-polarized
    2 : e0
    3 : e1
    .
    .
    .
    8 : e6
    """
    L = Lx*Ly
    istatus=lattice[k]
    if istatus == 0 : return 0,prob # nothing to do here
    # 
    # try to polarize
    #
    icell = min(istatus -1,1)
    prob +=  param[3]*(1-icell) # 1 - min(1,1-1) = 1
                                # 1 - min(1,2-1) = 0 ....
    if rng < prob:
        lattice[k] = np.random.randint(2,8)
        return 200+lattice[k],prob # polarized in this direction
    #
    # failed to polarize. if icell = 0, then whatever follows
    # should not add to the probability. 
    #
    # consider adding a if + return statement here
    #

    #
    # try to depolarize
    #
    prob += param[2]*icell
    if rng < prob:
        lattice[k] = 1 # depolarized state
        return 100,prob

    # idirection = 0, 1, ... ,5
    idirection = max(istatus - 2,0) # NOTE:: in case icell = 0
                                    # NOTE:: consider adding an if above
    target = select_next(k,idirection,Lx,Ly)
    # check if target is not in the border, otherwise
    # just return
    if not (target < L): return -target,prob # NOTE:: it should not be necessary
                                             # NOTE:: because border sites should not be selected
    shared=[lattice[j] for j in get_neighbours_shared(k,target,Lx,Ly,iopen=1)]
    #
    # try to move to target
    #
    istatus_target = min(1,lattice[target]) # is it occupied 
    istatus_shared = min(1,sum(shared))     # is the neighbourhood occupied
    #
    # OBS: move with neighbours around or not, does not change the 
    #      outcome. The probability, however, it does change.
    #
    prob += param[0]*(1-istatus_target)*icell*istatus_shared
    prob += param[1]*(1-istatus_target)*icell*(1-istatus_shared)

    if rng < prob:
        lattice[target] = lattice[k] # caution here (python-numpy)
        lattice[k]      = 0
        return 1+idirection,prob
    return 0,prob

def update_variabletime(k,Lx,Ly,lattice,rates):
    """
    updates the k-th site using the Gillespie algorithm
    and returns the resulting time increase

    OBS: the routine assumes the k-th site is occupied 
         by a cell, regardless of polarization
    OBS: rates refers to the single cell transition
         rates
    """
    iflag       = 0
    delta_time  = 0
    total_rate  = 0
    icoord = 6 # hexagon
    if (lattice[k] < 2) & (rates['polarization'] > 0e0):
        #
        # cell has no polarization
        #
        total_rate = rates['polarization']
        delta_time = -(1e0/total_rate)*np.log(np.random.rand())
        iflag = 100+np.random.choice(range(2,icoord),1) 
        return iflag,delta_time        
    #
    # cell has polarization
    #
    idirection = lattice[k] - 2 # NOTE: assuming lattice[k] > 1
    next_site  = select_next(k,idirection,Lx,Ly) 
    next_avail = 1 - min(lattice[next_site],1)        
    has_shared = min(1, sum(
        [lattice[j] for j in get_neighbours_shared(k,next_site,Lx,Ly,iopen=1)]))
    #
    # calculate the total transition rate
    #
    prob_dist = np.array([rates['gapjunction_positive']*next_avail*has_shared,
                          rates['gapjunction_negative']*next_avail*(1-has_shared),
                          rates['depolarization']],dtype=np.float64)
    total_rate  = np.sum(prob_dist) 
    #
    # generate a random number from the exponential distribution
    #
    if total_rate > 0:
        prob_dist   = prob_dist / total_rate
        #
        # select a transition according to the partial weights
        # of the rates
        #
        transition = int(np.random.choice(range(len(prob_dist)),
                                      p=prob_dist))
        #
        # use a dictionary to replace the switch case, given
        # lattice[k] and lattice[next_site]
        #
        # lattice[k],lattice[next_site] = \
        #     {0: lambda x,y: (0,x),
        #      1: lambda x,y: (0,x),
        #      2: lambda x,y: (1,y)}[transition](lattice[k],lattice[next_site])

        depol = int(transition > 1)
        lattice[next_site]= lattice[k]*(1-depol) + depol
        lattice[k]        = depol
        
        delta_time = -(1e0/total_rate)*np.log(np.random.rand())
        iflag = transition+200
        #print(iflag,type(iflag))
    return iflag,delta_time

#==================================================
def sample_fixedtime(Lx=5,Ly=10,steps=100,rates={}):
    """
    returns a sample lattice after steps timesteps have
    elapsed, given the transition rates rates, via 
    gillespie algorithm. Note that the time interval between
    timesteps are *not* equal.
    """
    lattice = np.zeros(Lx*Ly,dtype=int)
    boundary_conditions(Lx,Ly,lattice)
        
    data_density = np.zeros((steps,2),dtype=np.float64)
    dt     = 1e0/(Lx*Ly)
    params = [rates['gapjunction_positive']*dt,
              rates['gapjunction_negative']*dt,
              rates['depolarization'      ]*dt,
              rates['polarization'        ]*dt]
    
    for step in range(steps):
        prob = 0e0
        rng  = np.random.rand()

        data_density[step,:]=(step*dt,get_particles(lattice)) 
        #
        # list instead of a generator! slow down...
        #
        order = np.arange((Lx-1)*Ly)
        np.random.shuffle(order)
        for k in order:
            iflag,prob = update_fixedtime(k,Lx,Ly,lattice,rng,prob,params)
            if iflag>0: break
        boundary_conditions(Lx,Ly,lattice)
    return lattice,data_density

def sample_variabletime(Lx=5,Ly=10,steps=100,rates={}):
    """
    returns a sample lattice after steps timesteps have
    elapsed, given the transition rates rates, via 
    gillespie algorithm. Note that the time interval between
    timesteps are *not* equal.
    """
    lattice = np.zeros(Lx*Ly,dtype=int)
    boundary_conditions(Lx,Ly,lattice)        
    data_density = np.zeros((steps,2),dtype=np.float64)
    current = 0
    for step in range(steps):
        data_density[step,:]=(current,get_particles(lattice)) 
        k = np.random.choice( [j for j in range((Lx-1)*Ly) if lattice[j] > 0]  )
        iflag,dt = update_variabletime(k,Lx,Ly,lattice,rates)
        boundary_conditions(Lx,Ly,lattice)
        current += dt
    return lattice,data_density








#==================================================
if __name__ == "__main__":
    
    print("::Starting unitary tests::")
    def test_neighbours():
        """
        test all neighbours from a hexagonal lattice
        with Lx = 5 and Ly = 4 (open boundary)
    
        00  01  02  03  04
          05  06  07  08  09
        10  11  12  13  14
          15  16  17  18  19
        --------------------
        20  21  22  23  24
          20  21  22  23  24
        """
        Lx=5
        Ly=4
        neighs = {0: [1,5,9,4,24,20],
                  1: [2,6,5,0,20,21],
                  2: [3,7,6,1,21,22],
                  3: [4,8,7,2,22,23],
                  4: [0,9,8,3,23,24],
                  5: [6,11,10,9,0,1],
                  6: [7,12,11,5,1,2],
                  12:[13,17,16,11,6,7],
                  14:[10,19,18,13,8,9],
                  15:[16,21,20,19,10,11],
                  17:[18,23,22,16,12,13] }
        v = [ (x,all(neighs[x]==get_neighbours(x,Lx,Ly))) for x in neighs]
        #for entry in v:
        #    print ("... site %s :: %s" % (entry[0],entry[1]))
        return all([x[1] for x in v])
#--------------------------------------
    def test_update_fixedtime_1():
        istatus = 0

        Lx=5
        Ly=10
        L =Lx*Ly
        lattice = np.zeros(L,dtype=int)
        boundary_conditions(Lx,Ly,lattice)
        params = [0e0,0e0,0e0,1e0] #single cell will always polarize
        lattice[0] = istatus
        prob = 0e0
        rng  = 1e-12
        iflag,prob = update_fixedtime(0,Lx,Ly,lattice,rng,prob,params)
        return all( [iflag == 0, prob < 1e-10, lattice[0]==0 ])
#--------------------------------------
    def test_update_fixedtime_2():
        istatus = 1

        Lx=5
        Ly=10
        L =Lx*Ly
        lattice = np.zeros(L,dtype=int)
        boundary_conditions(Lx,Ly,lattice)
        params = [0e0,0e0,0e0,1e0] #single cell will always polarize
        params = np.array(params,dtype=np.float64)/L
        lattice[0] = istatus
        prob = 0e0
        rng  = 1e-12
        iflag,prob = update_fixedtime(0,Lx,Ly,lattice,rng,prob,params)
        return all( [iflag != 0, prob <= (1e0/L)+1e-10,(lattice[0] !=(istatus|0))])
#--------------------------------------
    def test_update_fixedtime_3():        
        iflag = False
        Lx=5
        Ly=10
        L = Lx*Ly

        params = [1e0,0e0,0e0,0e0] #single cell only moves if neighs around
        params = np.array(params,dtype=np.float64)/L
        
        lattice = np.zeros(L,dtype=int)
        boundary_conditions(Lx,Ly,lattice)
        prob = 0e0
        rng  = 0e0 # np.random.rand()
        iflag,prob = update_fixedtime(0,Lx,Ly,lattice,rng,prob,params)

        return all((iflag==2,lattice[5]==3,lattice[0]==0,prob<1e-6+params[0]))
#--------------------------------------
    def test_update_fixedtime_4():        
        iflag = False
        Lx=5
        Ly=10
        L = Lx*Ly

        params = [1e0,0e0,0e0,0e0] #single cell only moves if neighs around
        params = np.array(params,dtype=np.float64)/L
        
        lattice = np.zeros(L,dtype=int)
        boundary_conditions(Lx,Ly,lattice)
        prob = 0e0
        rng  = 0e0 # np.random.rand()
        lattice[5]=1
        iflag,prob = update_fixedtime(0,Lx,Ly,lattice,rng,prob,params)
        return all((iflag==0,lattice[5]==1,lattice[0]==3,prob<1e-6+params[0]))
#--------------------------------------
    def test_update_fixedtime_5():        
        iflag = False
        Lx=5
        Ly=10
        L = Lx*Ly

        params = [0.5e0,0.5e0,0e0,0e0] #does not matter if a neigh is present
        params = np.array(params,dtype=np.float64)/L
        
        lattice = np.zeros(L,dtype=int)
        boundary_conditions(Lx,Ly,lattice)
        prob = 0e0
        rng  = 0e0 # np.random.rand()
        lattice[11]=2
        iflag,prob = update_fixedtime(11,Lx,Ly,lattice,rng,prob,params)
        return all((iflag==1,lattice[11]==0,lattice[12]==2,prob<1e-6+params[0]))
#--------------------------------------
    def test_update_fixedtime_6():        
        iflag = False
        Lx=5
        Ly=10
        L = Lx*Ly

        params = [0e0,0e0,1e0,0e0] # cell always depolarize
        params = np.array(params,dtype=np.float64)/L
        
        lattice = np.zeros(L,dtype=int)
        boundary_conditions(Lx,Ly,lattice)
        prob = 0e0
        rng  = 0e0 # np.random.rand()
        iflag,prob = update_fixedtime(0,Lx,Ly,lattice,rng,prob,params)
        return all((iflag==100,lattice[0]==1,prob<1e-6+params[2]))
#--------------------------------------
    print(">>> Test neighbours :: %s" % (test_neighbours()))
    print(">>> Test update 1   :: %s" % (test_update_fixedtime_1() ))
    print(">>> Test update 2   :: %s" % (test_update_fixedtime_2() ))
    print(">>> Test update 3   :: %s" % (test_update_fixedtime_3() ))
    print(">>> Test update 4   :: %s" % (test_update_fixedtime_4() ))
    print(">>> Test update 5   :: %s" % (test_update_fixedtime_5() ))
    print(">>> Test update 6   :: %s" % (test_update_fixedtime_6() ))
