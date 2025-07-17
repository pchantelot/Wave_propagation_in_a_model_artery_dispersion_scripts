"""
    C_elastic(λ,μ)

Returns the 4-th order elastic tensor for a material with Lamé parameters `λ` and `μ`.
"""
function C_elastic(λ,μ)

  II = collect(I(3)) .* collect(reshape(I(3),1,1,3,3))
  c = λ .* II .+ μ .* (permutedims(II,(1, 3, 4, 2)) .+ permutedims(II,(1, 3, 2, 4)))

  return c

end

"""
    C_viscoacoustoelastic(lambda1,lambda3,vT,vL,rho,f,tau,n,alpha,betaprime)

Returns the 4-th order elastic tensor for an incompressible medium taking into acount pre-stress and viscoelastic losses.

input arguments:
- `lambda1`: stretch ratio in direction 1
- `lambda3`: stretch ratio in direction 3 (NB: lambda2 is computed from the incompressibility condition)
- `vT`: transverse wave velocity
- `vL`: longitudinal wave velocity
- `rho`: material density
- `f`: frequency
- `tau`: relaxation time
- `n`: fractional derivative order
- `alpha`: first hyperelastic parameter
- `betaprime`: second hyperelastic parameter
"""
function C_viscoacoustoelastic(lambda1,lambda3,vT,vL,rho,f,tau,n,alpha,betaprime)

    # Mooney-rivlin parameters
    C10 = (1-alpha)*rho*vT^2
    C01 = alpha*rho*vT^2
    # Fractional Kelvin-Voigt model parameters
    omega = 2*pi*f
    nu = (1-betaprime)*rho*vT^2*tau^n
    beta = betaprime*rho*vT^2*tau^n

    c = Array{ComplexF64}(undef,3,3,3,3)

    C0_1111=(1/9).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(
  3*vL.^2+(-4).*vT.^2).^(-1)).^(-7/3).*(14*C01.*lambda1.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(
  -1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+9*rho.*vL.^2*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+
  lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(10/3)+(-12).*rho.*vT.^2*(1+lambda1.^(-2).*
  lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(
  -1)).^(10/3)+C10.*(5*(lambda1.*lambda3).^(2/3).*(lambda1.^(-3).*lambda3.^(-3).*(lambda1.^2*lambda3.^2+(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*
  lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1))).^(8/3)+5*lambda3.^(8/3).*(lambda3.^(-1)+lambda1.^(-2).*lambda3.^(-3).*(C10.*((
  -2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(2/3)+2*
  lambda1.^(8/3).*(lambda1.^(-1)+lambda1.^(-3).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*
  lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(2/3))+(-1).*C01.*lambda1.^2*(lambda3.^2+lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+
  lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2)+18*((sqrt(complex(-1))*(-1)
  ).*omega).^n.*(beta.*lambda1.^(13/3).*(lambda1.^(-1)+lambda1.^(-3).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*
  lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(7/3)+nu.*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*
  lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(7/3)))

    C0_1122=(1/9).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(
    3*vL.^2+(-4).*vT.^2).^(-1)).^(-7/3).*(2*C10.*lambda3.^(8/3).*(lambda3.^(-1)+lambda1.^(-2).*lambda3.^(-3).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*
    lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(2/3)+(-4).*C01.*lambda1.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+ 
    lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+(-3).*rho.*(3*
    vL.^2+(-4).*vT.^2).*((-1)+lambda1.^(-2).*lambda3.^(-2).*((-2).*C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+2*C01.*(lambda1.^2+lambda3.^2+(-2).*lambda1.^4*
    lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*
    lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(7/3)+(-4).*C10.*((lambda1.*lambda3).^(2/3).*(lambda1.^(-3).*lambda3.^(-3).*(lambda1.^2*
    lambda3.^2+(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)))
    .^(8/3)+lambda1.^(8/3).*(lambda1.^(-1)+lambda1.^(-3).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*
    lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(2/3))+2*C01.*lambda1.^2*((-2).*lambda3.^2+lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(
    C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2))

    C0_1133=(1/9).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(
    3*vL.^2+(-4).*vT.^2).^(-1)).^(-7/3).*((-4).*C01.*lambda1.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+ 
    (-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+(-3).*rho.*(3*vL.^2+(-4).*vT.^2).*((-1)+lambda1.^(-2).*lambda3.^(-2).*((-2).*C10.*((-2)+
    lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+2*C01.*(lambda1.^2+lambda3.^2+(-2).*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).*(1+lambda1.^(-2).*lambda3.^(
    -2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(7/3)
    +C10.*(2*(lambda1.*lambda3).^(2/3).*(lambda1.^(-3).*lambda3.^(-3).*(lambda1.^2*lambda3.^2+(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*
    lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1))).^(8/3)+(-4).*lambda3.^(8/3).*(lambda3.^(-1)+lambda1.^(-2).*lambda3.^(-3).*(C10.*((-2)+lambda1.^4*
    lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(2/3)+(-4).*lambda1.^(8/3).*(
    lambda1.^(-1)+lambda1.^(-3).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(
    -1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(2/3))+C01.*lambda1.^2*(2*lambda3.^2+(-4).*lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+
    lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2))

    C0_1212=C10.*(lambda1.^6*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*
    rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-5)).^(1/3)+C01.*(lambda1.^(-6).*lambda3.^(-6).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*
    lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^7).^(-1/3)+(1/2).*((sqrt(complex(-1))*(-1)).*omega).^n.*(2*nu+
    beta.*(lambda1.^2+lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+
    2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2))

    C0_1221=(1/6).*(3*beta.*lambda1.^2*((sqrt(complex(-1))*(-1)).*omega).^n+6*nu.*((sqrt(complex(-1))*(-1)).*omega).^n+(-8).*rho.*vT.^2+(-6).*lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+
    lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*vL.^2*(3*vL.^2+(-4).*vT.^2).^(-1)+8*vT.^2*(rho+lambda1.^(-2).*lambda3.^(-2).*(C10.*
    ((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*(3*vL.^2+(-4).*vT.^2).^(-1))+3*beta.*lambda1.^(-2).*
    lambda3.^(-2).*((sqrt(complex(-1))*(-1)).*omega).^n.*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*
    lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+2*C10.*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1)
    .*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-5/3).*(lambda1.^2+lambda3.^2+lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2).*
    lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(
    -1)).^2)+C01.*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*
    rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-7/3).*(4*lambda1.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*
    lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+(-2).*lambda1.^2*((-2).*lambda3.^2+lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(
    -2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*
    vT.^2).^(-1)).^2)))

    C0_1313=(1/2).*(beta.*(lambda1.^2+lambda3.^2)+2*nu).*((sqrt(complex(-1))*(-1)).*omega).^n+C10.*(lambda1.^6*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*
    lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-5)).^(1/3)+C01.*(lambda3.^6+lambda1.^(-2).*lambda3.^4*
    (C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-1/3)

    C0_1331=(1/6).*(3*beta.*lambda1.^2*((sqrt(complex(-1))*(-1)).*omega).^n+3*beta.*lambda3.^2*((sqrt(complex(-1))*(-1)).*omega).^n+6*nu.*((sqrt(complex(-1))*(-1)).*omega).^n+(-8).*rho.*vT.^2+(-6).*lambda1.^(-2).*
    lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*vL.^2*(3*vL.^2+(-4).*vT.^2).^(-1)+
    8*vT.^2*(rho+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*(3*
    vL.^2+(-4).*vT.^2).^(-1))+2*C10.*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*
    lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-5/3).*(lambda1.^2+lambda3.^2+lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+
    lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2)+2*C01.*(1+
    lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(
    -4).*vT.^2).^(-1)).^(-7/3).*(2*lambda1.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+
    2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+lambda1.^2*((-1).*lambda3.^2+2*lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)
    +lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2)))

    C0_2121=C10.*(lambda1.^(-8).*lambda3.^(-8).*(lambda1.^2*lambda3.^2+(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*
    lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1))).^(1/3)+C01.*(lambda1.^6+lambda1.^4*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*
    lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-1/3)+(1/2).*((sqrt(complex(-1))*(-1)).*omega).^n.*(2*nu+beta.*(lambda1.^2+lambda1.^(-2)
    .*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*
    rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2))

    C0_2222=(1/9).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(
    3*vL.^2+(-4).*vT.^2).^(-1)).^(-7/3).*(C10.*(2*(lambda1.*lambda3).^(2/3).*(lambda1.^(-3).*lambda3.^(-3).*(lambda1.^2*lambda3.^2+(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*
    lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1))).^(8/3)+5*lambda3.^(8/3).*(lambda3.^(-1)+lambda1.^(-2).*
    lambda3.^(-3).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(
    -1)).^(2/3)+5*lambda1.^(8/3).*(lambda1.^(-1)+lambda1.^(-3).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*
    lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(2/3))+(-1).*C01.*(lambda1.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+
    lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+lambda1.^2*((-14).*lambda3.^2+lambda1.^(
    -2).*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4))
    .*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2))+3*(3*rho.*vL.^2*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*
    lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(10/3)+(-4).*rho.*vT.^2*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+
    lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(10/3)+6*((sqrt(complex(-1))*(
    -1)).*omega).^n.*(beta.*(lambda1.*lambda3).^(7/3).*(lambda1.^(-3).*lambda3.^(-3).*(lambda1.^2*lambda3.^2+(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*
    lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1))).^(13/3)+nu.*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+
    lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(7/3))))

    C0_2233=(1/9).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(
    3*vL.^2+(-4).*vT.^2).^(-1)).^(-7/3).*(2*C01.*lambda1.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(
    -1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+(-3).*rho.*(3*vL.^2+(-4).*vT.^2).*((-1)+lambda1.^(-2).*lambda3.^(-2).*((-2).*C10.*((-2)+
    lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+2*C01.*(lambda1.^2+lambda3.^2+(-2).*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).*(1+lambda1.^(-2).*lambda3.^(
    -2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(7/3)
    +2*C10.*((-2).*(lambda1.*lambda3).^(2/3).*(lambda1.^(-3).*lambda3.^(-3).*(lambda1.^2*lambda3.^2+(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+
    (-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1))).^(8/3)+(-2).*lambda3.^(8/3).*(lambda3.^(-1)+lambda1.^(-2).*lambda3.^(-3).*(C10.*((-2)+
    lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(2/3)+lambda1.^(8/3).*
    (lambda1.^(-1)+lambda1.^(-3).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(
    -1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(2/3))+(-4).*C01.*lambda1.^2*(lambda3.^2+lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+
    lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2))

    C0_2323=C10.*(lambda1.^(-8).*lambda3.^(-8).*(lambda1.^2*lambda3.^2+(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*
    lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1))).^(1/3)+C01.*(lambda3.^6+lambda1.^(-2).*lambda3.^4*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*
    lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-1/3)+(1/2).*((sqrt(complex(-1))*(-1)).*omega).^n.*(2*nu+beta.*(lambda3.^2+lambda1.^(-2)
    .*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*
    rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2))

    C0_2332=(1/6).*(3*beta.*lambda3.^2*((sqrt(complex(-1))*(-1)).*omega).^n+6*nu.*((sqrt(complex(-1))*(-1)).*omega).^n+(-8).*rho.*vT.^2+(-6).*lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+
    lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*vL.^2*(3*vL.^2+(-4).*vT.^2).^(-1)+8*vT.^2*(rho+lambda1.^(-2).*lambda3.^(-2).*(C10.*
    ((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*(3*vL.^2+(-4).*vT.^2).^(-1))+3*beta.*lambda1.^(-2).*
    lambda3.^(-2).*((sqrt(complex(-1))*(-1)).*omega).^n.*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*
    lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+2*C10.*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1)
    .*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-5/3).*(lambda1.^2+lambda3.^2+lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2).*
    lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(
    -1)).^2)+(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1)
    .*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-7/3).*((-2).*C01.*lambda1.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*
    lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+4*C01.*lambda1.^2*(lambda3.^2+lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2)
    .*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2)
    .^(-1)).^2)))

    C0_3131=(1/2).*(beta.*(lambda1.^2+lambda3.^2)+2*nu).*((sqrt(complex(-1))*(-1)).*omega).^n+C10.*(lambda3.^6*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*
    lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-5)).^(1/3)+C01.*(lambda1.^6+lambda1.^4*lambda3.^(-2).*
    (C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-1/3)

    C0_3232=C10.*(lambda3.^6*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*
    rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(-5)).^(1/3)+C01.*(lambda1.^(-6).*lambda3.^(-6).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*
    lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^7).^(-1/3)+(1/2).*((sqrt(complex(-1))*(-1)).*omega).^n.*(2*nu+
    beta.*(lambda3.^2+lambda1.^(-2).*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+
    2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2))

    C0_3333=(1/9).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(
    3*vL.^2+(-4).*vT.^2).^(-1)).^(-7/3).*(C10.*(5*(lambda1.*lambda3).^(2/3).*(lambda1.^(-3).*lambda3.^(-3).*(lambda1.^2*lambda3.^2+(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*
    lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1))).^(8/3)+2*lambda3.^(8/3).*(lambda3.^(-1)+lambda1.^(-2).*
    lambda3.^(-3).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(
    -1)).^(2/3)+5*lambda1.^(8/3).*(lambda1.^(-1)+lambda1.^(-3).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*
    lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(2/3))+(-1).*C01.*(lambda1.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+
    lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2+lambda1.^2*(lambda3.^2+(-14).*lambda1.^(
    -2).*lambda3.^(-2).*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4))
    .*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^2))+3*(3*rho.*vL.^2*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*
    lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(10/3)+(-4).*rho.*vT.^2*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+
    lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(10/3)+6*((sqrt(complex(-1))*(
    -1)).*omega).^n.*(beta.*lambda3.^(13/3).*(lambda3.^(-1)+lambda1.^(-2).*lambda3.^(-3).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*
    lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(7/3)+nu.*(1+lambda1.^(-2).*lambda3.^(-2).*(C10.*((-2)+lambda1.^4*lambda3.^2+lambda1.^2*
    lambda3.^4)+C01.*((-1).*lambda1.^2+(-1).*lambda3.^2+2*lambda1.^4*lambda3.^4)).*rho.^(-1).*(3*vL.^2+(-4).*vT.^2).^(-1)).^(7/3))))

    c[1,1,1,1]=C0_1111
    c[1,1,2,2]=C0_1122
    c[2,2,1,1]=C0_1122
    c[1,1,3,3]=C0_1133
    c[3,3,1,1]=C0_1133
    c[1,2,1,2]=C0_1212
    c[1,2,2,1]=C0_1221
    c[2,1,1,2]=C0_1221
    c[1,3,1,3]=C0_1313
    c[1,3,3,1]=C0_1331
    c[3,1,1,3]=C0_1331
    c[2,1,2,1]=C0_2121
    c[2,2,2,2]=C0_2222
    c[2,2,3,3]=C0_2233
    c[3,3,2,2]=C0_2233
    c[2,3,2,3]=C0_2323
    c[2,3,3,2]=C0_2332
    c[3,2,2,3]=C0_2332
    c[3,1,3,1]=C0_3131
    c[3,2,3,2]=C0_3232
    c[3,3,3,3]=C0_3333

    c = permutedims(c,(1, 2, 4, 3))

    return c

end