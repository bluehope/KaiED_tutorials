# ENV["PROJECT_PATH_ED"]="../envs/KED"
# include("../src/mybase.jl")
using KaiEDJ
using KaiEDJ: read_wann90, I, Plots
#using DFTforge

# hr  = read_wann90( "ex_wann_1d_2sites.dat" )
hr  = read_wann90( "ex_wann_1d_2sites.dat" )

@show hr

for hri in hr
    @show hri
end


nk      = 1000
# karr    = [ [x,0,0] for x in LinRange(0,2*pi,nk) ]
# wkarr   = [ 1. for x in LinRange(0,2*pi,nk) ]
xarr, wkarr = KaiEDJ.gausslegendre( nk )
karr    = [ [(x .+ 1) * pi, 0,0] for x in xarr ]
wkarr   = wkarr / 2.        # integrate an unit volume

norb    = 2
hk  = hr_to_hk( hr, norb, karr ; ISPIN=2 )        # Default (ISPIN=1) : dim(hk) = 2*norb   (with spin-doubling onto spinless wannier functions)
@show hk

for hki in hk
    @show hki
end

NReFreq   = 100
ReFreqGridVal   = collect(LinRange( -3, 3, NReFreq ))
epsilon         = 0.03
ReFreqGrid      = ReFreqGridVal .+ im * epsilon

chem    = collect( 0.0 * I(norb) ) 
glocal_w  = KaiEDJ.GetGreenLocalFromHkGrid( hk, wkarr, ReFreqGrid, chem )     # chemical potential convection : H_chem = - mu * N
gdiag_w = [ [ g[i,i] for g in glocal_w ] for i in 1:norb ]




p_g     = Plots.plot(  ReFreqGridVal, imag(gdiag_w[1]) )
p_g     = Plots.plot!( ReFreqGridVal, real(gdiag_w[1]) )


Eorb    = collect( 0.0 * I(norb) ) 
@time hyblatt_w  = KaiEDJ.GetHybFromGreenLocalFromHkGrid( hk, wkarr, ReFreqGrid, chem, Eorb )     # Approx. : inv(g0) = omega - Eorb + chem - Delta_hyb
println("hyblatt_w obtained.") ; flush(stdout)
@time hyblattdiag_w = [ [ g[i,i] for g in hyblatt_w ] for i in 1:norb ]
println("hyblattdiag_w obtained.") ; flush(stdout)

p_h     = Plots.plot(  ReFreqGridVal, imag(hyblattdiag_w[1]) )
p_h     = Plots.plot!( ReFreqGridVal, real(hyblattdiag_w[1]) )

Plots.plot( p_g, p_h )


#### Matsubara Green Function Setup ####
beta    = 128
NImFreq  = 1*beta
ImFreqGridVal   = KaiEDJ.GetImFreqValGrid( beta, NImFreq )
ImFreqGrid      = ImFreqGridVal * im

@time hyblatt_iw  = KaiEDJ.GetHybFromGreenLocalFromHkGrid( hk, wkarr, ImFreqGrid, chem, Eorb )
println("hyblatt_iw obtained.") ; flush(stdout)
@time hyblattdiag_iw = [ [ g[i,i] for g in hyblatt_iw ] for i in 1:norb ]
println("hyblattdiag_iw obtained.") ; flush(stdout)

# plot!( ImFreqGridVal, imag(hyblattdiag_iw[1]), marker=:circle )
