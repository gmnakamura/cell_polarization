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

