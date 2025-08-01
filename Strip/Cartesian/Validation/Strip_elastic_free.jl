# Guided wave dispersion in an elastic strip using SCM 
# Load packages
# Easy identity matrix and kronecker tensor notation
using LinearAlgebra, Kronecker
using FileIO
# Plotting 
using CairoMakie
# POnG
include("../../../POnG.jl")
using .POnG


## Problem definition: We look for waves propagating in direction 1 in a strip with 
#height h and width w in vaccuum.
#
#                               top
#  e2             ---------------------------------- 
#  |_ e3     left - solid: μs, λs, ρs              - right
#                 ---------------------------------- 
#                               bottom
# Waves propagate in direction 1, U(x,t) = u(y,z)e^{i(kx-ωt)}
# normalization using ρs, μs and h

# displacement dof
udof = 1:3 
# Solid properties
ρs = 1000
μs = 36e3
h = 3e-3
w = 4e-2
λs = 1000^2*ρs - 2*μs # measured speed of sound in Ecoflex 1000 m/s
# Define elastic tensor
C = C_elastic(λs, μs)
C = C ./ μs # normalization

# Discretization: We differentiate in direction 2 and 3 using SCM.
N = 8 # number of nodes in direction 2 (even as we look at in-plane modes)
P = 14 # number of nodes in direction 3
_, Dt1 = chebdif(N,2)
D2 = -2 * (h/h) .* Dt1[:,:,1]
D22 = 4 * (h/h)^2 .* Dt1[:,:,2]
_, Dt2 = chebdif(P,2)
D3 =  - 2 * h/w .* Dt2[:,:,1] # normalized with h
D33 = 4 *(h/w)^2 .* Dt2[:,:,2]
D23d = D3 ⊗ D2 
D2d = I(P) ⊗ D2
D22d = I(P) ⊗ D22
D3d = D3 ⊗ I(N)
D33d = D33 ⊗ I(N)

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
# (ik)^2 Lkk + (ik) Lk + L0 + ω^2 Md =0
Lkk = C11 ⊗ I(N*P); Lkk = collect(Lkk)
Lk = (C12 .+ C21) ⊗ D2d .+ (C13 .+ C31) ⊗ D3d
L0 = C22 ⊗ D22d .+ C33 ⊗ D33d .+ (C23 .+ C32) ⊗ D23d
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
B0t = C22 ⊗ D2d[top,:] .+ C23 ⊗ D3d[top,:]
# boundary traction operators bottom
Bkb = C21 ⊗ I(N*P)[bottom,:]
B0b = C22 ⊗ D2d[bottom,:] .+ C23 ⊗ D3d[bottom,:]
# boundary traction operators left
Bkl = C31 ⊗ I(N*P)[left,:]
B0l = C32 ⊗ D2d[left,:] .+ C33 ⊗ D3d[left,:]
# boundary traction operators right
Bkr = C31 ⊗ I(N*P)[right,:]
B0r = C32 ⊗ D2d[right,:] .+ C33 ⊗ D3d[right,:]

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

for i = 1:length(ω)
    local temp = ω[i] * h / sqrt(μs/ρs)
    local problem = PEP([temp^2 .* Md .+ L0, Lk, Lkk])
    local ktmp, utmp = polyeig(problem)
    k[i,:] = -1im * ktmp' / h
    u[i,:,:] = utmp
end

# Store result of a 2D calculation before (and for) postprocessing
result = result2D(udof,N,P,ω,k,u)
idx = modedirection(result)
sym = modesymmetry(result, dir = "ip")

## Plots
# remove negative real(k) and purely imaginary solutions
f = repeat(ω,1,Int(2*(length(udof)*N*P))) /(2*pi)
k[(1 .- abs.(imag.(k[:]))) .< 0] .= NaN
f[isnan.(k)] .= NaN
kplot = k[:]
fplot = f[:]
idxplot = idx[:]
symplot = sym[:]
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
        ax.title = "In-plane modes"
        scatter!(ax,real(kplot[idxplot .!== 2 .&& symplot .== 1]),fplot[idxplot .!== 2 .&& symplot .== 1], 
            color = :red, markersize = 10,  label = "symmetric")
        scatter!(ax,real(kplot[idxplot .!== 2 .&& symplot .== 0]),fplot[idxplot .!== 2 .&& symplot .== 0], 
            color = :blue, markersize = 10,  label = "antisymmetric")
        scatter!(ax,data["k"][:],data["f"][:], marker = :xcross, markersize = 6, label = "Comsol")
        colsize!(fig.layout, 1, Aspect(1, 1.2))
        axislegend(position = :rb)

    resize_to_layout!(fig)
    display(fig)
end
