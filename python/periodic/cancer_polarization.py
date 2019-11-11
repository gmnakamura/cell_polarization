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

def get_neighbours(k,Lx,Ly,iopen=0):
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

def get_neighbours_shared(k1,k2,Lx,Ly,iopen=0):
    """
    return the shared neighbours between sites
    k1 and k2.
    """
    tmp = np.intersect1d(get_neighbours(k1,Lx,Ly,iopen),
                         get_neighbours(k2,Lx,Ly,iopen))
    return tmp 

def select_next(k,idirection,Lx,Ly):
    """
    return the site the particle is supposed to move
    next

    idirection = 0,1,2,3,4,5 
    corresponds to e1,e2,...,e6 in the hexagonal lattice
    """
    return get_neighbours(k,Lx,Ly)[idirection]


def initial_conditions(Lx,Ly,N=1):
    n=N #int(Lx*Ly/4)
    k=0
    lattice = np.zeros(Lx*Ly,dtype=int)
    while n > 0:
        lattice[k]=1
        k += 2
        n += -1
    return lattice

def update_automata(k,Lx,Ly,lattice):
    iflag = 0
    next_site = k
    L = Lx*Ly
    available = [j for j in  get_neighbours(k,Lx,Ly) if lattice[j] < 1]
    if len(available)*lattice[k] > 0:
        next_site = np.random.choice(available)
        iflag = 1
        lattice[k] = 0
        lattice[next_site] = 1        
    return iflag,next_site


def sample_automata(Lx=10,Ly=20,N=20,steps=1):
    lattice = initial_conditions(Lx,Ly,N)
    cell = {'position':[k for k in range(Lx*Ly) if lattice[k]>0]}
    for step in range(steps):
        np.random.shuffle(cell['position'])
        for entry in range(len(cell['position'])):
            iflag,new = update_automata(cell['position'][entry],Lx,Ly,lattice)
            if iflag > 0 : cell['position'][entry] = new
    return lattice

def get_neighbourhood_frequency(Lx,Ly,lattice):
    ocup = np.zeros(Lx*Ly,dtype=int)
    for k in range(Lx*Ly):
        ocup[k] = sum(lattice[ get_neighbours(k,Lx,Ly)  ])
    hist = np.array([ sum(ocup == k) for k in range(6+1) ],dtype=np.float64)
    return hist / sum(hist)

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
