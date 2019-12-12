
using StatsBase
#
# Auxiliary globals
#
const _empty        = 0 
const _nonpolarized = 1

#
#  Auxiliary functions 
#

function update(k,cell,lattice,_rates,prob,rng)
    # k is 1-based index

    iflag  = 0
   
    ipos   = cell[k]
    istate = lattice[ipos]    
    # NOTE comment the double-check below
    if (istate < 1)
        return iflag,prob
    end    
    #
    # time interval
    #
    dt = 1e0/prod(size(lattice)) # detailed scale
    #    dt = 1e0/size(cell)[1] # only for particle conservation
    #
    #  ONE  BODY transitios
    #
    itmp = Int( istate == _nonpolarized )
    # using total polarization rate! i.e., 4 times the indiv. polarz
    prob += _rates["polarization"]*dt*itmp 
    if (rng < prob)
        # 4 is the coordination number. NOTE Use a variable 
        lattice[ipos] = rand((_nonpolarized+1):(_nonpolarized+4)) 
        return 1,prob
    end
    itmp = Int(istate > _nonpolarized)
    prob += _rates["depolarization"]*dt*itmp
    if (rng < prob)
        lattice[ipos] = _nonpolarized
        return 2,prob
    end
    # get_neighbours works for nonpolarized as well
    L = [i for i in size(lattice)]
    neigh = get_neighbours(ipos,L)[istate]
    # cell is polarized and the target site is empty
    itmp = itmp * Int( lattice[neigh] == _empty  )
    prob += _rates["displacement"]*itmp*dt    
    if (rng < prob)
        cell[k]        = neigh
        lattice[neigh] = istate
        lattice[ipos]  = _empty
        return 3,prob
    end
    return iflag,prob
end

function init(L,N0 = 1,inum_states=5)
    L_tot = prod(L)
    lattice = zeros(Int,L...)
    cell = StatsBase.sample(1:L_tot,N0,replace=false)
    lattice[cell] = rand(1:inum_states,N0)
    return cell,lattice 
end

function get_neighbours(k,L)
    kshift  = k - 1 # because julia uses 1-based indices
    idims   = size(L)[1]
    icoords = get_coords(kshift,L)
    neighs  = zeros(Int,2*idims+1)
    neighs[1]=kshift
    # forward neighs
    for i=1:idims
        ioutro = icoords[:]
        ioutro[i] = mod(icoords[i]+1+L[i],L[i])
        neighs[i+1] = get_index(ioutro,L) 
    end
    # back neighs
    for i=1:idims
        ioutro = icoords[:]
        ioutro[i] = mod(icoords[i]-1+L[i],L[i])
        neighs[i+idims+1] = get_index(ioutro,L) 
    end    
    return neighs + ones(Int,2*idims+1) # 1-based indices
end

function get_coords(k,L)
    # k is a 0-based index here    
    idims      = size(L)[1]
    icoords    = zeros(Int,idims)
    icoords[1] = mod(k,L[1])
    ifrac      = L[1]
    for j = 2:idims
        icoords[j] = div(k,ifrac)
        ifrac = ifrac*L[j]
    end
    return icoords # 0-based grid index
end
    
function get_index(icoords,L)
    idims  = size(L)[1]
    index = icoords[1]
    ifrac = L[1]
    for i=2:idims
        index += icoords[i]*ifrac
        ifrac = ifrac*L[i]
    end
    return index # 0-based index
end
