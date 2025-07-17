# legsrd
# Based on legsrd.m (SPECTRAL METHODS: Algorithms, Analysis and Applications, Jie Shen et al.)
# --- compute the Legendre n Legendre Gauss Radau points with x(1) = -1.
# refer to L. N. Trefethen: Spectral Methods in Matlab(2000)

function lgrpointsleft(n)

   # Eigenmethod is used for computing nodes.    
    A = Diagonal(1 ./ ((2*(0:n-2) .+ 1) .* (2*(0:n-2) .+ 3))) # Main diagonal
    A = A .+ SymTridiagonal(zeros(n-1), sqrt.((1:n-2).*((1:n-2) .+ 1)) ./ (2*(1:n-2) .+ 1)) # Create Jacobi matrix
    x = sort(eigvals(A)) # Compute eigenvalues
    x = [-1.0; x]

    return x

end