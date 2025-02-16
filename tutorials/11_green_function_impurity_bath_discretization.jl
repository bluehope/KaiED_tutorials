# ENV["PROJECT_PATH_ED"]="../envs/KED"
# include("../src/mybase.jl")
using KaiEDJ
using KaiEDJ: I , BenchmarkTools, Optimization, Plots

norb        = 1
nspin       = 2
nspinorb    = norb*nspin

nbath       = 5
@inline IndBath( ibath )            = ibath + nspinorb
@inline IndBathSpinup( ibathorb )   = 2*ibathorb + nspinorb - 1
@inline IndBathSpindn( ibathorb )   = 2*ibathorb + nspinorb

ntot    = nspinorb+nbath
dim     = 2^ntot
@show dim


ebathl     = zeros(nbath)
ebathl     = collect(LinRange( -1, 1, nbath ))
@show ebathl

Vil     = zeros(nspinorb,nbath)
for i in 1:nspinorb
    for j in 1:nbath
        Vil[i,j]    = 2 - ebathl[j]^2
    end
end
@show Vil

tij     = zeros(nspinorb,nspinorb)
tij[1,1]= 0


beta    = 32
NImFreq = 4*beta
ImFreqGridVal   = GetImFreqValGrid( beta, NImFreq )
ImFreqGrid      = ImFreqGridVal * im

G0iw= GetGzBethe.( ImFreqGrid )

Hyb11iw   = 1. / 4 * G0iw   # = t^2 G = 1/4 G


## Testing functions of flattened-parameters ##
BParam  = BathParamFlatten( ebathl, Vil )
Diw  = GetDeltaHybDiscGrid( ebathl, Vil, ImFreqGrid )
Diw2 = GetDeltaHybDiscGridFromFlat( BParam, ImFreqGrid, nspinorb, nbath )
D11iw    = GetijarrayFromVecMat( Diw, 1, 1 ) 
D11iw2   = GetijarrayFromVecMat( Diw2, 1, 1 ) 

Hybiw        = 1. / 4 * GetGzBetheDim.( ImFreqGrid, 2 )

SParam  = ( 
                Hybiw, 
                ImFreqGrid, 
                nspinorb, 
                nbath 
                )
@show typeof(SParam)
@show typeof(BParam)

# @inline function GetCostFromFlat( bathparam, systemparam ) 
#     HybwOrig    = systemparam[1]
#     wgrid       = systemparam[2]
#     dimorb      = systemparam[3]
#     dimbath     = systemparam[4]
#     HybwParam   = GetDeltaHybDiscGridFromFlat( bathparam, wgrid, dimorb, dimbath )
#     return GetCost( HybwOrig, HybwParam )
# end

# using BenchmarkTools
@time cost = GetCostFromFlat( BParam, SParam )
println( "initial cost : $(cost) " )


# using Optimization
using KaiEDJ.Optimization
prob = KaiEDJ.Optimization.OptimizationProblem(GetCostFromFlat, BParam, SParam)
# using OptimizationOptimJL
using KaiEDJ.OptimizationOptimJL
@time sol = KaiEDJ.OptimizationOptimJL.solve(prob, NelderMead(), maxiters=4000)
@show sol.original

using KaiEDJ.Optim
@time res = KaiEDJ.Optim.optimize( x -> GetCostFromFlat(x,SParam), BParam, LBFGS(), KaiEDJ.Optim.Options(iterations=4000))
@show res

# lb  = [-20.0 for i in BParam ]
# ub  = -lb
# prob = OptimizationProblem(GetCostFromFlat, BParam, SParam, lb = lb, ub = ub)
# sol = solve(prob, Fminbox(GradientDescent()))
# @show sol.original

# using ForwardDiff
# optf = OptimizationFunction( GetCostFromFlat, Optimization.AutoForwardDiff())
# prob = OptimizationProblem(optf, BParam, SParam)
# @show optf(BParam,SParam)
# sol = solve(prob, BFGS())
# @show sol.original



# BParamNew   = [ sol... ]  # From NelderMead
BParamNew   = [ res.minimizer... ]    # From BFGS
@time cost = KaiEDJ.GetCostFromFlat( BParamNew, SParam )
println( "final cost : $(cost) " )

enew, vnew  = KaiEDJ.BathParamReshape( BParamNew, nbath )
@show enew
println("vnew :")
using DelimitedFiles
writedlm(stdout, vnew)

DiwNew      = KaiEDJ.GetDeltaHybDiscGridFromFlat( BParamNew, ImFreqGrid, nspinorb, nbath )
GiwNew      = KaiEDJ.GetGreenLocalFromHybGrid( DiwNew, ImFreqGrid, collect(I(nspinorb)*0.0) )

x   = ImFreqGridVal
y1  = Hyb11iw
y2  = KaiEDJ.GetijarrayFromVecMat( DiwNew, 1, 1 )



## ReFreq spectra ##

NReFreq = 100
epsilon = 0.04
ReFreqGridVal   = collect( LinRange( -3, 3, NReFreq ) )
ReFreqGrid      = ReFreqGridVal .+ im * epsilon

Gw          = KaiEDJ.GetGzBethe.( ReFreqGrid )
DwNew       = KaiEDJ.GetDeltaHybDiscGridFromFlat( BParamNew, ReFreqGrid, nspinorb, nbath )
Gwn         = KaiEDJ.GetGreenLocalFromHybGrid( DwNew, ReFreqGrid, collect(I(nspinorb)*0.0) )

bPlot = true
if bPlot
    # using Plots
    xw  = ReFreqGridVal
    gw  = KaiEDJ.GetijarrayFromVecMat( Gw, 1, 1 )
    gwn = KaiEDJ.GetijarrayFromVecMat( Gwn, 1, 1 ) 

    Plots.plot(  xw, [ real(gw)  imag(gw)  ], label="G (Gw)" )
    Plots.plot!( xw, [ real(gwn) imag(gwn) ], label="G (Gwn)", xlabel="Frequency", ylabel="Green's function", title="Green's Function Comparison")
    Plots.savefig("Green_Function_Comparison.png")
end


