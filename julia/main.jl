
#  Thu Dec 12 09:49:16 CET 2019

#
# Simulation: cluster formation for cells with polarized movements
#

# auxiliary codes
include("aux.jl")
using Random
#
# PARAMETERS
#



# NOTE:: add a function to collect command line parameters
const _L = [32,32]
const _steps   = 30000
const _samples = 1
const _icoordination = 4
_rates   = Dict()
_rates["depolarization"]= 0.05e0
_rates["polarization"]  = 0.1e0  # total polariz rate ( 4 times larger )
_rates["displacement"]  = 1e0

#
# DATA INIT
#
iskip=10
data = zeros(Float64,Int(_steps/iskip),2)
L_tot = prod(_L)
#
# MONTE CARLO SAMPLES
#
for sample=1:_samples
    n0 = 64 #Int(prod(_L)*0.1)
    cell,lattice = init(_L,n0);
    idata_idx = 0
    for step=0:_steps-1
#        print(step,lattice[cell],'\n')
        N = size(cell)[1]
        # data collection
        if (mod(step,iskip) == 0)
            idata_idx += 1
            x = sum( lattice .== _nonpolarized)
            y = sum( lattice .>  _nonpolarized)
            data[idata_idx,1] += (x)/_samples
            data[idata_idx,2] += (y)/_samples
        end
        # shuffle the cell vector
        shuffle!(cell);
        # update
        iflag = 0
        k = 1
        prob=0e0
        rng = rand() # use MersenneTwister       
        while ((iflag < 1) & (k < N+1))
            iflag,prob = update(k,cell,lattice,_rates,prob,rng)           
            k += 1
        end
    end
end
