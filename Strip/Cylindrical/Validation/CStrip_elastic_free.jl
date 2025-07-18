# Strip in cylindrical coordinates
# Load packages
# Easy identity matrix and kronecker tensor notation
using LinearAlgebra, Kronecker
# Plotting 
using CairoMakie
using FileIO
# POnG
include("../../../POnG.jl")
using .POnG

## Problem definition: We look for Lamb waves propagating in direction 1 in a curved strip. 
# 
# Strip with height h, width w and curvature radius R in vacuum
#
#                           top, zero normal traction
#  e2 (= r)             ---------------------------------- 
#  |_ e3 (= θ)     left - solid: μs, λs, ρs              - right, zero displacement
#                       ---------------------------------- 
#                           bottom, zero normal traction
#
# Waves propagate in direction 1 (= x),  U(x,t) = u(r,θ)e^{i(kx-ωt)}
# normalization using ρs, μs and h

# displacement dof, coupled Lamb and SH waves
udof = 1:3 
# Solid properties 
ρs = 1000
μs = 36e3
λs = 1000^2*ρs - 2*μs # measured speed of sound in Ecoflex 1000 m/s
h = 3e-3
w = 4e-2
R = 1 # radius of curvature
# Define elastic tensor
C = C_elastic(λs, μs)
C = C ./ μs # normalization

# Discretization: We differentiate in direction 2 and 3 using SCM.
N = 8 # number of nodes in direction 2, must be even for in-plane modes and odd for out-of-plane modes
P = 14 # number of nodes in direction 3
x2, Dt1 = chebdif(N,2) # Nodes in [-1, 1]
x2 = (R .+ (h .* x2) ./ 2) ./ h
D2 = 2 .* (h/h) * Dt1[:,:,1]
D22 = 4 .* (h/h)^2 * Dt1[:,:,2]
_, Dt2 = chebdif(P,2)
D3 = 2 * (R/w) .* Dt2[:,:,1] 
D33 = 4 *(R/w)^2 .* Dt2[:,:,2]
A = [0 0 0; 0 0 -1; 0 1 0] # Derivation of the basis vectors

# Material matrices used in the elastodynamic equation
C11 = C[1,udof,udof,1]
C12 = C[1,udof,udof,2]
C21 = C[2,udof,udof,1]
C13 = C[1,udof,udof,3]
C31 = C[3,udof,udof,1]
C22 = C[2,udof,udof,2]
C23 = C[2,udof,udof,3]
C32 = C[3,udof,udof,2]
C33 = C[3,udof,udof,3]
M = collect(I(length(udof)))

# Discretized matrices to build the polynomial eigenvalue problem
# (ik)^2 Lkk + (ik) Lk + L0 + ω^2 Md = 0
Lkk = C11 ⊗ I(N*P); Lkk = collect(Lkk)
Lk = (C12 .+ C21) ⊗ I(P) ⊗ D2 .+ (C21 ⊗ I(P) .+ (C13 .+ C31) ⊗ D3 .+ (C13*A .+ A*C31) ⊗ I(P)) ⊗ Diagonal(1 ./ x2)
L0 = C22 ⊗ I(P) ⊗ D22 .+ (C22 ⊗ I(P) .+ C32 ⊗ D3 .+ A*C32 ⊗ I(P)) ⊗ (Diagonal(1 ./ x2)*D2) .+ 
    (C23 ⊗ D3 .+ C23*A ⊗ I(P)) ⊗ (D2*Diagonal(1 ./ x2) .+ Diagonal(1 ./ x2.^2)) .+
    (C33 ⊗ D33 .+ (C33*A .+ A*C33) ⊗ D3 .+ A*C33*A ⊗ I(P)) ⊗ Diagonal(1 ./ x2.^2)
Md = M ⊗ I(N*P); Md = collect(Md)

# boundary conditions
# N×P vector to N×P matrix
bc = reshape(1:N*P,N,P)
top = bc[1,:]
bottom = bc[N,:]
left = bc[:,1]
right = bc[:,P]
# boundary traction operators top
Bkt = C21 ⊗ I(N*P)[top,:]
B0t = C22 ⊗ (I(P) ⊗ D2)[top,:] .+ C23 ⊗ (D3 ⊗ Diagonal(1 ./ x2))[top,:] .+ C23*A ⊗ (I(P) ⊗ Diagonal(1 ./ x2))[top,:]
# boundary traction operators bottom
Bkb = C21 ⊗ I(N*P)[bottom,:]
B0b = C22 ⊗ (I(P) ⊗ D2)[bottom,:] .+ C23 ⊗ (D3 ⊗ Diagonal(1 ./ x2))[bottom,:] .+ C23*A ⊗ (I(P) ⊗ Diagonal(1 ./ x2))[bottom,:]
# zero displacement left
Bkl = C31 ⊗ I(N*P)[left,:]
B0l = C32 ⊗ (I(P) ⊗ D2)[left,:] .+ C33 ⊗ (D3 ⊗ Diagonal(1 ./ x2))[left,:] .+ C33*A ⊗ (I(P) ⊗ Diagonal(1 ./ x2))[left,:]
# zero displacement right
Bkr = C31 ⊗ I(N*P)[right,:]
B0r = C32 ⊗ (I(P) ⊗ D2)[right,:] .+ C33 ⊗ (D3 ⊗ Diagonal(1 ./ x2))[right,:] .+ C33*A ⊗ (I(P) ⊗ Diagonal(1 ./ x2))[right,:]

# Impose boundary conditions
topall = (0:length(udof)-1)*N*P .+ top'; topall = topall'[:]
bottomall = (0:length(udof)-1)*N*P .+ bottom'; bottomall = bottomall'[:]
leftall = (0:length(udof)-1)*N*P .+ left'; leftall = leftall'[:]
rightall = (0:length(udof)-1)*N*P .+ right'; rightall = rightall'[:]

# top
Lkk[topall,:] .= 0
Lk[topall,:] .= Bkt
L0[topall,:] .= B0t
Md[topall,:] .= 0
# bottom
Lkk[bottomall,:] .= 0
Lk[bottomall,:] .= Bkb
L0[bottomall,:] .= B0b
Md[bottomall,:] .= 0
# left
Lkk[leftall,:] .= 0
Lk[leftall,:] .= Bkl
L0[leftall,:] .= B0l
Md[leftall,:] .= 0
# right
Lkk[rightall,:] .= 0
Lk[rightall,:] .= Bkr
L0[rightall,:] .= B0r
Md[rightall,:] .= 0

## Solve polynomial eigenvalue problem
# Impose ω and get eigenvalues and eigenvectors
ω = 2*pi .* range(1e-3, 400, length = 120) 
k = Array{ComplexF64,2}(undef,length(ω),2*length(udof)*N*P)
u = Array{ComplexF64,3}(undef,length(ω),length(udof)*N*P,2*length(udof)*N*P)

@time for i = 1:length(ω)
    local temp = ω[i] * h / sqrt(μs/ρs)
    local problem = PEP([temp^2 .* Md .+ L0, Lk, Lkk])
    local ktmp, utmp = polyeig(problem)
    k[i,:] = -1im * ktmp' / h
    u[i,:,:] = utmp
end

# Store result of a 2D calculation before (and for) postprocessing
result = result2D(udof,N,P,ω,k,u)
idx = modedirection(result)

## Plots
# remove negative real(k) and purely imaginary solutions
f = repeat(ω,1,Int(2*(length(udof)*N*P))) /(2*pi)
k[(1 .- abs.(imag.(k[:]))) .< 0] .= NaN
f[isnan.(k)] .= NaN
kplot = k[:]
fplot = f[:]
idxplot = idx[:]
# Load Comsol Data
data = load(joinpath(@__DIR__,"comsol_dispersion_free.jld2"))

with_theme(theme_latexfonts()) do 
    x = 0:1500
    fig = Figure(size = (1.2 * 360, 360))

        ax = Axis(fig[1,1])
        ax.limits = (0, 300, 0, 400)
        ax.xlabel = L"k \,\, \mathrm{(1/m)}"
        ax.ylabel = L"f \,\, \mathrm{(Hz)}"
        ax.xlabelsize = 20 
        ax.ylabelsize = 20
        ax.title = "In plane modes"
        scatter!(ax,real(kplot[idxplot .!== 2]),fplot[idxplot .!== 2],
            markersize = 10, color = :red)
        scatter!(ax,data["k"][:],data["f"][:], marker = :xcross, markersize = 6, label = "Comsol")
        colsize!(fig.layout, 1, Aspect(1, 1.2))

    resize_to_layout!(fig)
    display(fig)
end