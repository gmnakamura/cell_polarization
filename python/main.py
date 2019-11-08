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


def update(k,Lx,Ly,lattice,rng,prob,param):
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
    if not (target < L): return -target,prob
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
#==================================================
def run_sample(Lx=5,Ly=10,steps=100,params=[1e0,0e0,0e0,0e0]):
    """
    returns the lattice and density data after 
    steps steps have elapsed
    """
    lattice = np.zeros(Lx*Ly,dtype=int)
    boundary_conditions(Lx,Ly,lattice)
        
    data_density = np.zeros(steps,dtype=np.float64)

    for step in range(steps):
        prob = 0e0
        rng  = np.random.rand()

        data_density[step]=get_particles(lattice) 
        #
        # list instead of a generator! slow down...
        #
        order = np.arange(Lx*Ly)
        np.random.shuffle(order)
        for k in order:
            iflag,prob = update(k,Lx,Ly,lattice,rng,prob,params)
            if iflag>0: break
        boundary_conditions(Lx,Ly,lattice)
    return lattice,data_density
#==================================================
#==================================================
def pretty_printing(Lx,Ly,lattice):
    """
    print the lattice for display in regular terminal
    """
    for ky in range(Ly):
        print(" "*(ky%2),*(lattice[i +ky*Lx] for i in range(Lx)))
    return None
#==================================================
def get_args():
    import argparse
    parser = argparse.ArgumentParser(description='''
    Stochastic simulation of cell dynamics

    a) gap junction model
    b) cell polarization

    ''')
    parser.add_argument('-x',metavar='Lx',nargs=1,type=int,
                default=5,help='lattice size in x-direction (default 5)')
    parser.add_argument('-y',metavar='Ly',nargs=1,type=int,
                default=10,help='lattice size in x-direction (default 10)')
    parser.add_argument('--samples',metavar='num samples',nargs=1,type=int,
                default=1,help='number of sample runs (default 1)')
    parser.add_argument('--steps',metavar='num steps',nargs=1,type=int,
                default=10,help='number of time steps (default 10)')
    parser.add_argument('--contact',metavar='rate',nargs=1,type=float,
                default=0.9,
                help='single cell | transition rate | gap junction (default 0.9)')
    parser.add_argument('--depolarization',metavar='rate',
                nargs=1,type=float, default=0.01,
                help='single cell | transition rate | depolarization (default 0.01)')
    parser.add_argument('--polarization',metavar='rate',
                nargs=1,type=float, default=0.01,
                help='single cell | transition rate | polarization (default 0.9)')
    args = parser.parse_args()
    return args
#==================================================
#==================================================
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
    def test_update_1():
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
        iflag,prob = update(0,Lx,Ly,lattice,rng,prob,params)
        return all( [iflag == 0, prob < 1e-10, lattice[0]==0 ])
#--------------------------------------
    def test_update_2():
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
        iflag,prob = update(0,Lx,Ly,lattice,rng,prob,params)
        return all( [iflag != 0, prob <= (1e0/L)+1e-10,(lattice[0] !=(istatus|0))])
#--------------------------------------
    def test_update_3():        
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
        iflag,prob = update(0,Lx,Ly,lattice,rng,prob,params)

        return all((iflag==2,lattice[5]==3,lattice[0]==0,prob<1e-6+params[0]))
#--------------------------------------
    def test_update_4():        
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
        iflag,prob = update(0,Lx,Ly,lattice,rng,prob,params)
        return all((iflag==0,lattice[5]==1,lattice[0]==3,prob<1e-6+params[0]))
#--------------------------------------
    def test_update_5():        
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
        iflag,prob = update(11,Lx,Ly,lattice,rng,prob,params)
        return all((iflag==1,lattice[11]==0,lattice[12]==2,prob<1e-6+params[0]))
#--------------------------------------
    def test_update_6():        
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
        iflag,prob = update(0,Lx,Ly,lattice,rng,prob,params)
        return all((iflag==100,lattice[0]==1,prob<1e-6+params[2]))
#--------------------------------------
    print(">>> Test neighbours :: %s" % (test_neighbours()))
    print(">>> Test update 1   :: %s" % (test_update_1() ))
    print(">>> Test update 2   :: %s" % (test_update_2() ))
    print(">>> Test update 3   :: %s" % (test_update_3() ))
    print(">>> Test update 4   :: %s" % (test_update_4() ))
    print(">>> Test update 5   :: %s" % (test_update_5() ))
    print(">>> Test update 6   :: %s" % (test_update_6() ))
