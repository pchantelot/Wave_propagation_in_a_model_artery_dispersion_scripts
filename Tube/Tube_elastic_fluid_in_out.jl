# Guided wave dispersion in a tube with fluid inside and outside
# Load packages
# Easy identity matrix and kronecker tensor notation
using LinearAlgebra, Kronecker
# Hankel functions for BC
using SpecialFunctions
# Plotting 
using CairoMakie
# POnG
include("../POnG.jl")
using .POnG

## Problem definition: We look for waves propagating in direction 1 in a cylinder with mean radius R with wall thickness h.
#
#                           top, fluid
#  e2 (= r)             ---------------------------------- 
#  |_ e3 (= θ)          - solid: μs, λs, ρs              - 
#                       ---------------------------------- 
#                           bottom, fluid until r = 0

# Waves propagate in direction 1 (= x),  U(x,t) = u(r)e^{i(kx+nθ-ωt)}
# normalization using λf, ρs and h

# displacement dof
udof = 1:3
# Solid properties 
ρs = 1070
μs = 23e3
λs = 1500^2 * ρs - 2 * μs
R = 4e-3 #radius of curvature
h = 1e-3
n = 1
# Fluid properties 
ρf = 1000
vf = 1500
μf = μs * 1e-6 # We model the fluid as a solid with μ → 0
λf = vf^2 * ρf - 2 * μf

# Define elastic tensor
C = C_elastic(λs, μs)
C = C ./ λf # normalization
# Define elastic "fluid" tensor
Cf = C_elastic(λf, μf)
Cf = Cf ./ λf # normalization

# Discretization in solid
N = 20 # number of nodes in direction 2
x2, D = chebdif(N, 2) # Nodes in [-1, 1]
x2 = (R .+ (h .* x2) ./ 2) ./ h
D2 = 2 .* (h/h) * D[:,:,1]
D22 = 4 .* (h/h)^2 * D[:,:,2]
A = [0 0 0; 0 0 -1; 0 1 0] # Derivation of the basis vectors
AA = (1im * n * I(length(udof)) .+ A) 

# Solid Material matrices
C11 = C[1,udof,udof,1]
C12 = C[1,udof,udof,2]
C21 = C[2,udof,udof,1]
C13 = C[1,udof,udof,3]
C31 = C[3,udof,udof,1]
C22 = C[2,udof,udof,2]
C23 = C[2,udof,udof,3]
C32 = C[3,udof,udof,2]
C33 = C[3,udof,udof,3]
M = ρs / ρs * I(length(udof))

# Discretized matrices to build the polynomial eigenvalue problem
# (ik)^2 Lkk + (ik) Lk + L0 + ω^2 Md =0
# solid
Lkk = C11 ⊗ I(N); Lkk = collect(Lkk)
Lk = (C12 .+ C21) ⊗ D2 .+ (C21 .+ C13*AA .+ AA*C31) ⊗ Diagonal(1 ./ x2) 
L0 = C22 ⊗ D22 .+ (C22 .+ C23*AA .+ AA*C32) ⊗ (Diagonal(1 ./x2)*D2) .+ AA*C33*AA ⊗ Diagonal(1 ./ x2.^2)
Md = M ⊗ I(N); Md = collect(Md); Md = convert(Array{ComplexF64},Md)
# boundary traction operators top
Bkt = C21 ⊗ I(N)[1,:]'
B0t = C22 ⊗ D2[1,:]' .+ C23*AA ⊗ Diagonal(1 ./ x2)[1,:]'
# boundary traction operators bottom
Bkb = C21 ⊗ I(N)[N,:]'
B0b = C22 ⊗ D2[N,:]' .+ C23*AA ⊗ Diagonal(1 ./ x2)[N,:]'
# displacement bottom 
D0b = I(length(udof)) ⊗ I(N)[N,:]'
# Outer fluid addition
e2 = zeros(length(udof)); e2[2] = 1
Bρft = 1im*ρf/ρs * I(3).*e2 ⊗ I(N)[1,:]'  # approximate BC to get a ploynomial eigenvalue problem

# Discretization in fluid
P = 15
xf2 = lgrpointsleft(P)
Df = lgrdiffleft(P,xf2)
xf2 = (1/2 * (R - h/2) .+ (R - h/2) .* -xf2 ./ 2 ) ./ h
Df2 = - 2 .* (h/(R - h/2)) * Df # normalization !
Df22 = 4 * (h/(R - h/2))^2 * Df*Df

# Fluid Material matrices
Cf11 = Cf[1,udof,udof,1]
Cf12 = Cf[1,udof,udof,2]
Cf21 = Cf[2,udof,udof,1]
Cf13 = Cf[1,udof,udof,3]
Cf31 = Cf[3,udof,udof,1]
Cf22 = Cf[2,udof,udof,2]
Cf23 = Cf[2,udof,udof,3]
Cf32 = Cf[3,udof,udof,2]
Cf33 = Cf[3,udof,udof,3]
Mf = ρf / ρs * I(length(udof))

# Build fluid polynomial eigenvalue problem
Lfkk = Cf11 ⊗ Diagonal(xf2.^2); Lfkk = collect(Lfkk)
Lfk = (Cf12 .+ Cf21) ⊗ (Diagonal(xf2.^2)*Df2) .+  Cf21 ⊗ Diagonal(xf2) .+ 
    Cf13*AA ⊗ Diagonal(xf2) .+ AA*Cf31 ⊗ Diagonal(xf2)
Lf0 = Cf22 ⊗ (Diagonal(xf2.^2)*Df22) .+ Cf22 ⊗ (Diagonal(xf2)*Df2)  .+ Cf23*AA ⊗ (Diagonal(xf2)*Df2) .+ 
    AA*Cf32 ⊗ (Diagonal(xf2)*Df2) .+ AA*Cf33*AA ⊗ I(P)
Mfd = Mf ⊗ Diagonal(xf2.^2); Mfd = collect(Mfd)
# boundary traction operators top
Bfkt = Cf21 ⊗ I(P)[1,:]'
Bf0t = Cf22 ⊗ Df2[1,:]' .+ Cf23*AA ⊗ Diagonal( 1 ./ xf2)[1,:]'
# displacement top 
Df0t = I(length(udof)) ⊗ I(P)[1,:]'

# Impose solid BC: 
bc =  [(0:length(udof)-1).*N .+ 1, (1:length(udof)).*N]
# top
Lkk[bc[1],:] .= Bkt
Lk[bc[1],:] .= B0t
L0[bc[1],:] .= 0
Md[bc[1],:] .= Bρft
# bottom
Lkk[bc[2],:] .= 0
Lk[bc[2],:] .= Bkb
L0[bc[2],:] .= B0b
Md[bc[2],:] .= 0

# Impose fluid BCs
bcf =  [(0:length(udof)-1).*P .+ 1, (1:length(udof)).*P]
# top
Lfkk[bcf[1],:] .= 0
Lfk[bcf[1],:] .= 0
Lf0[bcf[1],:] .= Df0t
Mfd[bcf[1],:] .= 0

# Assemble matrices
Fkk, Fk, F0, Fd = assemblelayer(Lkk,Lk,L0,Md,Lfkk,Lfk,Lf0,Mfd,bc,bcf,D0b,Bfkt,Bf0t)

## Solve polynomial eigenvalue problem
# Impose ω and get eigenvalues and eigenvectors
ω = 2 * pi * range(1e-2, 400, length = 70) 
k = Array{ComplexF64,2}(undef, length(ω), 2 * length(udof) * (N + P))
u = Array{ComplexF64,3}(undef, length(ω), length(udof) * (N+P), 2 * length(udof) *(N+P))

for i  in axes(ω,1)
    temp = ω[i] * h / sqrt(λf/ρs) 
    problem = PEP([temp^2 .* Fd .+ F0, Fk, Fkk])
    ktmp, utmp = polyeig(problem)
    k[i,:] = -1im * ktmp' / h
    u[i,:,:] = utmp
end

# Store result of a 1D calculation before (and for) postprocessing
result = result1D(udof,N,ω,k,u)
# remove negative real(k) and purely imaginary solutions
#k[(0.49*pi .< angle.(k)) .| (angle.(k) .< -0.49*pi)] .= NaN 
# Remove evanescent solutions
k[(1 .- abs.(imag.(k[:]))) .< 0] .= NaN
k[real(k) .> 1500 .|| real(k).< 0] .= NaN
## Plots
kplot = k[:]
fplot = repeat(ω, 1, Int(2 * length(udof) * (N+P)))[:] / (2 * pi)
fplot[isnan.(kplot)] .= NaN
kplot = filter(!isnan,kplot)
fplot = filter(!isnan,fplot)

# Perform an iterative process to apply the accurate BC
# Initialisation of the iterative process
kplot2 = Vector{ComplexF64}()
# Make blocks once to reassemble Fd
TRB = zeros(ComplexF64, size(Lkk)[1], size(Lfkk)[1])
BLB = zeros(ComplexF64, size(Lfkk)[1], size(Lkk)[1])
# Loop on all (ω,k) eigenvalues couples
niter = 15
for i in axes(kplot,1)
    # current (ω,k)
    ktmp = kplot[i]
    ωtmp = fplot[i] * 2*pi
    # find associated eigenvector
    idxf = argmin(abs.(ω .- ωtmp))
    uf = u[idxf,:,:]
    kf = k[idxf,:]
    uf = uf[:,isfinite.(kf)]
    kf = kf[isfinite.(kf)]
    idxk = argmin(abs.(real.(kf) .- ktmp))
    utmp = uf[:,idxk]

    #compute the new BC
    tol = 1
    tic = 1
    while tol > 1e-3 && tic < niter
        tic = tic + 1
        ky = sqrt((ω[idxf] * h / vf)^2 - (ktmp*h)^2)
        if imag(ky) < 0
            ky = conj(ky)
        end
        ro = (R + h/2) / h
        Bρftup =  I(3).*e2 ⊗ I(N)[1,:]' * (1im*ktmp*h) / ky *
            (ρf/ρs * hankelh1(n, ky*ro)/(n/(ky*ro) * hankelh1(n, ky*ro) - hankelh1(n+1, ky*ro))) 
        # update BC
        Md[bc[1],:] .= Bρftup
        # reassemble Fd
        tt = hcat(Md,TRB)
        tb = hcat(BLB,Mfd)
        Fdup = vcat(tt,tb)
        # helper matrices for inverse iteration
        A1 = (1im*ktmp*h)^2*Fkk .- F0 .- (ω[idxf] * h / sqrt(λf/ρs))^2*Fdup
        A2 = 2*1im*ktmp*h*Fkk .+ Fk

        utmp = (A1 .- 1im*ktmp*h*A2) \ A2*utmp
        sol = (A2*utmp \ A1*utmp) * -1im /h 
        tol = abs(real(ktmp) - real(sol))
        ktmp = sol
    end
    push!(kplot2,ktmp)
end

with_theme(theme_latexfonts()) do 
    x = 0:1000
    fig = Figure()
        ax = Axis(fig[1,1])
        ax.limits = (0, 1100, 0, 400)
        ax.xlabel = L"k \; \mathrm{(1/m)}"
        ax.ylabel = L"f \; \mathrm{(Hz)}"
        ax.xlabelsize = 20
        ax.ylabelsize = 20
            scatter!(ax, real(kplot), fplot, marker = :utriangle)
            scatter!(ax, real(kplot2), fplot)
            lines!(ax,x, sqrt(3*μs/ρs) * x / (2*pi), color = :black, linestyle = :dash, linewidth = 3, label = L"\sqrt{E_s/ρ_s}")
            lines!(ax,x, sqrt(μs/ρs) * x / (2*pi), color = :black, linewidth = 3, label = L"\sqrt{μ_s/ρ_s}")
            lines!(ax, x, x .* sqrt(3*μs/ρf * h / R / 2) / (2*pi), linewidth = 3, color = :green, label = L"v_b")
            lines!(ax,x, 1500 * x / (2*pi), linewidth = 3, color = :red, label = L"v_f")
        axislegend(position = :rb)
        colsize!(fig.layout, 1, Aspect(1, 1.2))
        ax2 = Axis(fig[1,2])
        ax2.limits = (0, 2, 0, 400)
        ax2.xlabel = L"k \; \mathrm{(1/m)}"
        ax2.ylabel = L"f \; \mathrm{(Hz)}"
        ax2.xlabelsize = 20
        ax2.ylabelsize = 20
            scatter!(ax2, real(kplot), fplot, marker = :utriangle)
            scatter!(ax2, real(kplot2), fplot)
            lines!(ax2,x, sqrt(3*μs/ρs) * x / (2*pi), color = :black, linestyle = :dash, linewidth = 3, label = L"\sqrt{E_s/ρ_s}")
            lines!(ax2,x, sqrt(μs/ρs) * x / (2*pi), color = :black, linewidth = 3, label = L"\sqrt{μ_s/ρ_s}")
            lines!(ax2, x, x .* sqrt(3*μs/ρf * h / R / 2) / (2*pi), linewidth = 3, color = :green, label = L"v_b")
            lines!(ax2,x, 1500 * x / (2*pi), linewidth = 3, color = :red, label = L"v_f")
        colsize!(fig.layout, 2, Aspect(1, 1.2))
    resize_to_layout!(fig)
    display(fig)
end

