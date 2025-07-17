# lgrdiffleft
# Based on lgrdiff.m (SPECTRAL METHODS: Algorithms, Analysis and Applications, Jie Shen et al.)
# --- compute the first order differentiation matrix of size n × n associated to the Legendre Gauss Radau points with x(1) = -1.
# refer to L. N. Trefethen: Spectral Methods in Matlab(2000)

function lgrdiffleft(n,x)
  if n == 0 

    D=[]

  else 

    nx = size(x)
    y1, dy1 = lpoly(n-1,x)
    y2, dy2 = lpoly(n,x)
    dy = dy1 .+ dy2 # Compute L_{n-1}+L_n and its 1st-order derivative
    y = y1 .+ y2
    D = x ./ dy * dy' - (1. ./ dy) * (x .* dy)' # Compute dy(x_j) (x_k-x_j)/dy(x_k), 1/d_{kj} for k not= j (see (3.204))
    D = D .+ I(n) # add the identity matrix so that 1./D can be operated  
    D = 1 ./ D
    D = D .- I(n)
    diagvec = [-(n+1)*(n-1)/4; x[2:end] ./ (1 .- x[2:end].^2) .+ n * y1[2:end] ./ ((1 .- x[2:end].^2) .* dy[2:end])]
    D = D .+ Diagonal(diagvec)

  end

  return D

end