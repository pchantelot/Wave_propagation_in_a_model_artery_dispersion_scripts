# Helper function for differentiation on the Legendre Gauss Radau points
# Based on lepoly.m (SPECTRAL METHODS: Algorithms, Analysis and Applications, Jie Shen et al.)
# --- compute the Legendre polynomial of degree n at the value given in the vectorx.
# refer to L. N. Trefethen: Spectral Methods in Matlab(2000)

function lpoly(n,x)
    
    if n == 0

        y = ones(length(x))
        dy = zeros(length(x))

    elseif n == 1

        y = x
        dy = ones(length(x))

    else
    
        polylst = ones(length(x))
        pderlst = zeros(length(x))
        poly = x 
        pder = ones(length(x))
        # initialise
        polyn = zeros(length(x))
        pdern = zeros(length(x))


            for k = 2:n # Three-term recurrence relation:  

                polyn = ((2*k - 1)* x .* poly - (k-1) * polylst) / k # kL_k(x)=(2k-1)xL_{k-1}(x)-(k-1)L_{k-2}(x)
                pdern = pderlst + (2*k-1) * poly # L_k'(x)=L_{k-2}'(x)+(2k-1)L_{k-1}(x)
                polylst = poly
                poly = polyn
                pderlst = pder 
                pder = pdern

            end

        y = polyn 
        dy = pdern

    end
    
    return y, dy

end