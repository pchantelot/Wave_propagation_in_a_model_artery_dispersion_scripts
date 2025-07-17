# Guided wave dispersion in a tube 
# Load packages
# Easy identity matrix and kronecker tensor notation
using LinearAlgebra, Kronecker
# Plotting 
using CairoMakie
# File IO
using FileIO
# POnG
path_to_repository = "Documents/GitHub/Wave_propagation_in_a_model_artery_dispersion_scripts/"
push!(LOAD_PATH,joinpath(pwd(),path_to_repository))
using POnG

## Problem definition: We look for waves propagating in direction 1 in a cylinder with mean radius R with wall thickness h.
#
#                           top, zero normal traction
#  e2 (= r)             ---------------------------------- 
#  |_ e3 (= θ)          - solid: μs, λs, ρs              - 
#                       ---------------------------------- 
#                           bottom, zero normal traction

# Waves propagate in direction 1 (= x),  U(x,t) = u(r)e^{i(kx+nθ-ωt)}
# normalization using ρs, μs and h

# displacement dof
udof = 1:3 # displacement dof
# Solid properties 
ρs = 7932
μs = ρs*3260^2
λs = ρs*5960^2 - 2 * μs
R = 20.5 #radius of curvature
h = 1
n = 1
# Define elastic tensor
C = C_elastic(λs, μs)
C = C ./ μs # normalization

# Discretization: We differentiate in direction 2 using SCM.
N = 20 # number of nodes in direction 2
x2, D = chebdif(N, 2) # Nodes in [-1, 1]
x2 = (R .+ (h .* x2) ./ 2) ./ h
D2 = 2 .* (h/h) * D[:,:,1]
D22 = 4 .* (h/h)^2 * D[:,:,2]
A = [0 0 0; 0 0 -1; 0 1 0] # Derivation of the basis vectors
AA = (1im * n * I(length(udof)) .+ A) 

# Material matrices
C11 = C[1,udof,udof,1]
C12 = C[1,udof,udof,2]
C21 = C[2,udof,udof,1]
C13 = C[1,udof,udof,3]
C31 = C[3,udof,udof,1]
C22 = C[2,udof,udof,2]
C23 = C[2,udof,udof,3]
C32 = C[3,udof,udof,2]
C33 = C[3,udof,udof,3]
M = I(length(udof))

# Discretized matrices to build the polynomial eigenvalue problem
# (ik)^2 Lkk + (ik) Lk + L0 + ω^2 Md =0
Lkk = C11 ⊗ I(N); Lkk = collect(Lkk)
Lk = (C12 .+ C21) ⊗ D2 .+ (C21 .+ C13*AA .+ AA*C31) ⊗ Diagonal(1 ./ x2) 
L0 = C22 ⊗ D22 .+ (C22 .+ C23*AA .+ AA*C32) ⊗ (Diagonal(1 ./x2)*D2) .+ AA*C33*AA ⊗ Diagonal(1 ./ x2.^2)
Md = M ⊗ I(N); Md = collect(Md)
# boundary traction operators top
Bkt = C21 ⊗ I(N)[1,:]'
B0t = C22 ⊗ D2[1,:]' .+ C23*AA ⊗ Diagonal(1 ./x2)[1,:]'
# boundary traction operators bottom
Bkb = C21 ⊗ I(N)[N,:]'
B0b = C22 ⊗ D2[N,:]' .+ C23*AA ⊗ Diagonal(1 ./x2)[N,:]'

# Impose BC: location of boundary vectors in discretized matrix
bc =  [(0:length(udof)-1).*N .+ 1, (1:length(udof)).*N]
# top
Lkk[bc[1],:] .= 0
Lk[bc[1],:] = Bkt
L0[bc[1],:] = B0t
Md[bc[1],:] .= 0
# bottom
Lkk[bc[2],:] .= 0
Lk[bc[2],:] = Bkb
L0[bc[2],:] = B0b
Md[bc[2],:] .= 0

## Solve polynomial eigenvalue problem
# Impose ω and get eigenvalues and eigenvectors
ω = 2 * pi * sqrt(μs/ρs) / h * range(1e-6, 5e-2, length = 70) 
k = Array{ComplexF64,2}(undef, length(ω), 2 * length(udof) * N)
u = Array{ComplexF64,3}(undef, length(ω), length(udof) * N, 2 * length(udof) * N)

for i = 1:length(ω)
    local temp = ω[i] * h / sqrt(μs/ρs) 
    local problem = PEP([temp^2 .* Md .+ L0, Lk, Lkk])
    local ktmp, utmp = polyeig(problem)
    k[i,:] = -1im * ktmp' / h
    u[i,:,:] = utmp
end

# Store result of a 1D calculation before (and for) postprocessing
result = result1D(udof,N,ω,k,u)
# remove negative real(k) and purely imaginary solutions
k[(0.49*pi .< angle.(k)) .| (angle.(k) .< -0.49*pi)] .= NaN 
# Remove purely evanescent solutions
k[(1 .- abs.(imag.(k[:]))) .< 0] .= NaN
## Plots
kplot = k[:]
fplot = repeat(ω, 1, Int(2 * length(udof) * N))[:] / (2 * pi)
alpha = 1 .- abs.(imag.(kplot)) / (0.1 * maximum(abs.(imag.(filter(!isnan,kplot)))))

with_theme(theme_latexfonts()) do 
    fig = Figure()
        ax = Axis(fig[1,1])
        ax.limits = (0, 0.5, 0, 0.05)
        ax.xlabel = L"kh"
        ax.ylabel = L"fh/c_t"
        ax.xlabelsize = 20
        ax.ylabelsize = 20
        scatter!(ax, real(kplot)*h, fplot * h / sqrt(μs/ρs), color = tuple.(:red,alpha))
        colsize!(fig.layout, 1, Aspect(1, 1.2))
        resize_to_layout!(fig)
    display(fig)
end
