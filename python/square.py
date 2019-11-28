import numpy as np

_Lx = 64
_Ly = 64

_base = 2**(_Lx-1)
_mask = 2**(_Lx) - 1

def lshift(k):
    return ( (k<<1) + (k//_base) ) & _mask

def rshift(k):
    return ( (k>>1) + (k%2)*_base) & _mask

def compare(a,b=0,kind=0):
    """
    return bit fields in 'a' that can move to
    empty bit fields in 'b'

    'kind' specifies the type of action
     0: XOR(a,b) and then AND with a
     1: check bits in 'a' that can move to the right
    -1: check bits in 'a' that can move to the left
    """
    return { 0: (a ^        b ) & a,
             1: (a ^ lshift(a)) & a,
            -1: (a ^ rshift(a)) & a}[kind] & _mask

def pretty_string(a):
    return format(a,'b').rjust(_L,'0')
    
def count_bits(a):
    return format(a,'b').count('1')

def update(lines,rates):
    prob = 0.0
    rng = np.random.rand()
    k=0
    iflag = 0
    while ((k<_Ly) and (iflag < 1)):
        kp = (k+1+_Ly)%_Ly
        km = (k-1+_Ly)%_Ly
        stack = [rshift(lines[k]),lshift(lines[k]),
                 lines[kp],lines[km]]
        for entry in stack:
            v     = compare(lines[k],entry)
            icount= count_bits(v)
            prob += icount*rates
            if rng < prob:
                iselected = np.random.randint(icount)
                split = list(pretty_string(v))
        
        k +=1
    return lines

def init():
    return [1 & _mask for k in range(_Ly)]

def sample(steps,rates=0.1):
    lines = init()
    for k in range(steps):
        lines = update(lines,rates)
    return lines


