# Post-processing functions for SCM
using POnG

"""
	modedirection(p::result1D)

Compute an array of integers `idx` associating each eigenvalue of a SCM 1D computation `p` 
to its principal displacement direction (i.e. 1, 2, or 3).

"""
function modedirection(p::result1D)

    # Initialize array to associate a mode displacement direction to each eigenvalue
    idx = Array{Int,2}(undef, size(p.k)[1], size(p.k)[2])
    # Discretization 
    Drange = [(0:length(p.udof)-1) * p.N .+ 1 (1:length(p.udof)) * p.N]

    for i in axes(p.k,1)

        temp = [dropdims(sum(abs.(p.u[i,Drange[j,1]:Drange[j,2],:]).^2; dims = 1); dims = 1) for j =1:size(Drange)[1]]
        temp = stack(temp)'
        _, tempm = findmax(temp, dims = 1)
        idxtmp = getindex.(tempm[:],[1])
        idx[i,:] = idxtmp

    end

    return idx
    
end

"""
	modedirection(s::result2D)

Compute an array of integers `idx` associating each eigenvalue of a SCM 2D computation `s` 
to its principal displacement direction (i.e. 1, 2, or 3).
"""
function modedirection(s::result2D)

    # Initialize array to associate a mode displacement direction to each eigenvalue
    idx = Array{Int,2}(undef, size(s.k)[1], size(s.k)[2])
    # Discretization 
    Drange = [(0:length(s.udof)-1) * s.N * s.P .+ 1 (1:length(s.udof)) * s.N * s.P]

    for i in axes(s.k,1)

        temp = [dropdims(sum(abs.(s.u[i,Drange[j,1]:Drange[j,2],:]).^2; dims = 1); dims = 1) for j =1:size(Drange)[1]]
        temp = stack(temp)'
        _, tempm = findmax(temp, dims = 1)
        idxtmp = getindex.(tempm[:],[1])
        idx[i,:] = idxtmp

    end

    return idx
    
end

"""
	modesymmetry(p::result1D)

Determine the symmetry `sym` (1 if symmetric, 0 if antisymmetric) of the displacement field associated  
to each eigenvalue of a SCM 1D computation `p`.
"""
function modesymmetry(p::result1D)

    # Initialize array to associate a mode displacement direction to each eigenvalue
    sym = Array{Int,2}(undef, size(p.k)[1], size(p.k)[2])
    # Discretization
    Drange = [(0:length(p.udof)-1) * p.N .+ 1 (1:length(p.udof)) * p.N]

    for i in axes(p.k,1)

        utmp =  p.u[i,Drange[1,1]:Drange[1,2],:] # working on u1 
        utmps = utmp .+ reverse(utmp; dims = 1)
        utmpa = utmp .- reverse(utmp; dims = 1)
        sym[i,:] = sum(abs.(utmps); dims = 1) .> sum(abs.(utmpa); dims = 1) 

    end

    return sym
    
end

"""
	modesymmetry(s::result2D; dir)

Determine the symmetry `sym` (1 if symmetric, 0 if antisymmetric) of the displacement field associated  
to each eigenvalue of a SCM 2D computation `s`. Specify the Keyword `dir="ip"` (respectively `dir="oop"`) 
when looking to determine the symmetry of in-plane (respectively out-of-plane) modes.
"""
function modesymmetry(s::result2D; dir = "oop")

    # Initialize array to associate a mode displacement direction to each eigenvalue
    sym = Array{Int,2}(undef, length(s.ω), size(s.k)[2])
    # Discretization
    Drange = [(0:length(s.udof)-1) * s.N * s.P .+ 1 (1:length(s.udof)) * s.N * s.P]

    if dir == "ip"

        for i in axes(s.ω,1)

            utmp =  s.u[i,Drange[1,1]:Drange[1,2],:] # working on u1 

            for j in axes(utmp,2)

                # select a mode
                utmpmode = utmp[:,j]
                # reshape before making sym and antisym polarisation
                utmpmode = reshape(utmpmode, s.N, s.P)
                utmps = sum(utmpmode; dims = 1) .+ reverse(sum(utmpmode; dims = 1); dims = 2)
                utmpa = sum(utmpmode; dims = 1) .- reverse(sum(utmpmode; dims = 1); dims = 2)
                sym[i,j] = (sum(abs.(utmps).^2) .> sum(abs.(utmpa).^2))[1]

            end

        end

    end

    if dir == "oop"

        for i in axes(s.ω,1)

            utmp =  s.u[i,Drange[3,1]:Drange[3,2],:] # working on u3

            for j in axes(utmp,2)

                # select a mode
                utmpmode = utmp[:,j]
                # reshape before making sym and antisym polarisation
                utmpmode = reshape(utmpmode, s.N, s.P)
                utmps = sum(utmpmode; dims = 2) .+ reverse(sum(utmpmode; dims = 2); dims = 1)
                utmpa = sum(utmpmode; dims = 2) .- reverse(sum(utmpmode; dims = 2); dims = 1)
                sym[i,j] = (sum(abs.(utmps).^2) .> sum(abs.(utmpa).^2))[1]

            end

        end

    end

    return sym
    
end

"""
	modepolarisation(p::result1D)

Compute the average displacement polarisation `pol` ∈ [0, 1] where 0 is a purely in-plane polarisation 
and 1 a purely out-of-plane polarisation for the displacement field associated  
to each eigenvalue of a SCM 1D computation `p`.
"""
function modepolarisation(p::result1D)

    # Initialize array to associate a mode displacement direction to each eigenvalue
    pol = Array{Float64,2}(undef,length(p.ω), size(p.k)[2])
    # Discretization
    Drange = [(0:length(p.udof)-1).*p.N .+ 1 (1:length(p.udof)).*p.N]

    for i in axes(p.ω,1)

        utmp = p.u[i,:,:]
        umag = zeros(Float64, p.N, size(p.k)[2])
        inmag = zeros(Float64, p.N, size(p.k)[2])

        for j = 1:size(Drange)[1]
            # L2 norm od displacements in each direction
            umag += real(conj(utmp[Drange[j,1]:Drange[j,2],:]) .* utmp[Drange[j,1]:Drange[j,2],:])
        end
        umag = sum(umag; dims = 1) ./ p.N

        # case 1: udof = 1:3
        if length(p.udof) > 2
            for j = 1:2:size(Drange)[1]
                inmag += real(conj(utmp[Drange[j,1]:Drange[j,2],:]) .* utmp[Drange[j,1]:Drange[j,2],:])
            end
        # case 2: udof = 1:2 (NB: does not make sense if udof = 3...) 
        else
            inmag = real(conj(utmp[Drange[1,1]:Drange[1,2],:]) .* utmp[Drange[1,1]:Drange[1,2],:])
        end
        inmag = sum(inmag; dims = 1) ./ p.N

        pol[i,:] =  (umag .- inmag) ./ umag 

    end

    return pol
end

"""
	modepolarisation(s::result2D)

Compute the average displacement polarisation `pol` ∈ [0, 1] where 0 is a purely in-plane polarisation 
and 1 a purely out-of-plane polarisation for the displacement field associated  
to each eigenvalue of a SCM 2D computation `s`.
"""
function modepolarisation(s::result2D)

    # Initialize array to associate a mode displacement direction to each eigenvalue
    pol = Array{Float64,2}(undef,length(s.ω), size(s.k)[2])
    # Discretization 
    Drange = [(0:length(s.udof)-1) * s.N * s.P .+ 1 (1:length(s.udof)) * s.N * s.P]

    for i in axes(s.ω,1)

        utmp = s.u[i,:,:]
        umag = zeros(Float64, s.N * s.P, size(s.k)[2])
        inmag = zeros(Float64, s.N * s.P, size(s.k)[2])

        for j = 1:size(Drange)[1]
            # L2 norm od displacements in each direction
            umag += real(conj(utmp[Drange[j,1]:Drange[j,2],:]) .* utmp[Drange[j,1]:Drange[j,2],:])
        end
        umag = sum(umag; dims = 1) ./ (s.N * s.P)

        # case 1: udof = 1:3
        if length(s.udof) > 2
            for j = 1:2:size(Drange)[1]
                inmag += real(conj(utmp[Drange[j,1]:Drange[j,2],:]) .* utmp[Drange[j,1]:Drange[j,2],:])
            end
        # case 2: udof = 1:2 (NB: does not make sense if udof = 3...) 
        else
            inmag = real(conj(utmp[Drange[1,1]:Drange[1,2],:]) .* utmp[Drange[1,1]:Drange[1,2],:])
        end
        inmag = sum(inmag; dims = 1) ./ (s.N * s.P)

        pol[i,:] =  (umag .- inmag) ./ umag 

    end

    return pol
end


"""
	modesurfacepolarisation(p::result1D)

Compute the surface displacement polarisation `pol` ∈ [0, 1] where 0 is a purely in-plane polarisation 
and 1 a purely out-of-plane polarisation for the displacement field associated  
to each eigenvalue of a SCM 1D computation `p`.
"""
function modesurfacepolarisation(p::result1D)

    # Initialize array to associate a mode displacement direction to each eigenvalue
    pol = Array{Float64,2}(undef,length(p.ω), size(p.k)[2])
    # Discretization
    Drange = [(0:length(p.udof)-1).*p.N .+ 1 (1:length(p.udof)).*p.N]

    for i in axes(p.ω,1)

        utmp = p.u[i,:,:]
        umag = zeros(Float64,  size(p.k)[2])
        inmag = zeros(Float64, size(p.k)[2])

        for j = 1:size(Drange)[1]
            # L2 norm od displacements in each direction
            umag += real(conj(utmp[Drange[j,1],:]) .* utmp[Drange[j,1],:])
        end

        # case 1: udof = 1:3
        if length(p.udof) > 2
            for j = 1:2:size(Drange)[1]
                inmag += real(conj(utmp[Drange[j,1],:]) .* utmp[Drange[j,1],:])
            end
        # case 2: udof = 1:2 (NB: does not make sense if udof = 3...) 
        else
            inmag = real(conj(utmp[Drange[1,1],:]) .* utmp[Drange[1,1],:])
        end

        pol[i,:] =  (umag .- inmag) ./ umag 

    end

    return pol

end
