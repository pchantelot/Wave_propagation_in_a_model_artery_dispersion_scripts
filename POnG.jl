# POnG module: Programme pour les Ondes Guidées
module POnG

    # Load Julia's implementation of DMSUITE
    using LinearAlgebra, ToeplitzMatrices
    include(joinpath(@__DIR__,raw"Differentiation/chebpoints.jl"))
    include(joinpath(@__DIR__,raw"Differentiation/chebdif.jl"))
    include(joinpath(@__DIR__,raw"Differentiation/lpoly.jl"))
    include(joinpath(@__DIR__,raw"Differentiation/lgrpointsleft.jl"))
    include(joinpath(@__DIR__,raw"Differentiation/lgrdiffleft.jl"))
    export chebdif, cheb1extrema, lgrpointsleft, lgrdiffleft

    # Load Elastic tensor
    include(joinpath(@__DIR__,raw"Material/createC.jl"))
    export C_viscoacoustoelastic, C_elastic
    # Layer coupling
    include(joinpath(@__DIR__,raw"Material/assemblelayer.jl"))
    export assemblelayer

    # Load NonlinearEigenproblems to use polyeig
    using NonlinearEigenproblems
    export polyeig, PEP

    # Make 1D and 2D discretization class of results for use in processing functions.
    struct result1D
        udof::UnitRange{Int64}
        N::Int64
        ω::Array{Float64}
        k::Array{ComplexF64}
        u::Array{ComplexF64}
    end
    struct result2D
        udof::UnitRange{Int64}
        N::Int64
        P::Int64
        ω::Array{Float64}
        k::Array{ComplexF64}
        u::Array{ComplexF64}
    end
    export result1D, result2D

    # Postprocessing functions
    include(joinpath(@__DIR__,raw"Processing/SCMprocessing.jl"))
    export modedirection, modesymmetry

end