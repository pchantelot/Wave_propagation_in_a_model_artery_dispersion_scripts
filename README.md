# Wave propagation in a model artery: dispersion scripts [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16037558.svg)](https://doi.org/10.5281/zenodo.16037558)

This repository contains the scripts used to compute the dispersion relation of axial guided elastic waves in infinite tubes and (curved) strips reported in:
> P. Chantelot, A. Delory, C. Prada, and F. Lemoult, Wave propagation in a model artery, [10.48550/arXiv.2507.17698](https://doi.org/10.48550/arXiv.2507.17698) (2025).

We provide elasticity tensors to obtain the dispersion curves for elastic media (described by their Lamé parameters $\lambda$ and $\mu$) and for almost incompressible viscoelastic media subjected to an homogeneous finite deformation, which we account for using a Mooney-Rivlin model and a fractional Kelvin-Voigt model, within the acoustoelastic framework.

This work extends the codes developed in:
- [GEW dispersion script](https://github.com/dakiefer/GEW_dispersion_script), doi: [10.5281/zenodo.7010603](https://doi.org/10.5281/zenodo.7010603) 
- [GEW soft strip](https://github.com/dakiefer/GEW_soft_strip), doi: [10.5281/zenodo.10372576](http://doi.org/10.5281/zenodo.10372576)

## Usage
1. Download the repository.
2. Install the following [Julia](https://github.com/JuliaLang) packages: `FileIO`, `SpecialFunctions`, `LinearAlgebra`, `ToeplitzMatrices`, `Kronecker`,`NonlinearEigenproblems`, and `CairoMakie`. 
3. Execute the desired script, *e.g.* `Tube_elastic.jl`.

We include scripts to compute the dispersion relation in the following cases:
- elastic tube in vacuum, `Tube_elastic.jl`.
- elastic tube filled and surrounded by a semi-infinite fluid medium, `Tube_elastic_fluid_in_out.jl`.
- elastic free strip, `Strip_elastic_free.jl`.
- elastic clamped strip, `Strip_elastic_clamped.jl`.
- viscoelastic clamped strip subjected to a large homogeneous deformation, `Strip_viscoacoustoelastic_clamped.jl`.
- curved elastic clamped strip, `CStrip_elastic_clamped.jl`.
- curved viscoelastic clamped strip subjected to a large homogeneous deformation, `CStrip_viscoacoustoelastic_clamped.jl`.
- curved viscoelastic clamped strip subjected to a large homogeneous deformation with approximate coupling to a semi-infinite fluid medium, `CStrip_viscoacoustoelastic_clamped_fluidapprox.jl`.

> [!TIP]
> If you want to implement a different elasticity tensor (*i.e.* different hyperelastic model, viscoelastic model) have a look at the Mathematica files provided in the [GEW soft strip](https://github.com/dakiefer/GEW_soft_strip) repository. 

## Methods
We are interested in the dispersion relation of an homogeneous elastic waveguide, infinite in the direction of propagation, with a two-dimensional cross-section. 
We solve the three-dimensional elastodynamics equation,  
```math
    \boldsymbol{\nabla} \cdot (\boldsymbol{C} : \boldsymbol{\nabla} \boldsymbol{u}) = \rho \ddot{\boldsymbol{u}} \tag{1}
```
where $\boldsymbol{C}$ is the elasticity tensor, $\rho$ is the density, and $\boldsymbol{u}$ is the displacement field, subjected to the appropriate boundary conditions at the edges of the waveguide.
After inserting a propagative displacement ansatz, the dispersion relation is given by the non-trivial solutions of this boundary value problem.

We detail below how to obtain these non-zero solutions using the spectral collocation method that relies on discretizing the cross-section of the waveguide, transforming the boundary value problem into an polynomial eigenvalue problem.
This method has already been presented for plates[^1] and strips with a rectangular cross-section [^2]. We focus here on its implementation in cylindrical coordinates for tubes and curved strips, as well as on the implementation of fluid coupling.

We recall that in cylindrical coordinates $` \boldsymbol{\nabla} = \boldsymbol{e}_x \partial_x + \boldsymbol{e}_r \partial_r  + \boldsymbol{e}_\theta \frac{1}{r} \partial_\theta `$, so that,
```math
    \boldsymbol{\nabla u} = \boldsymbol{e}_x \partial_x \otimes \boldsymbol{u} + \boldsymbol{e}_r \partial_r \otimes \boldsymbol{u} + \boldsymbol{e}_\theta \frac{1}{r} \otimes \boldsymbol{A} \cdot \boldsymbol{u}, \tag{2}
```
where $`\boldsymbol{A} = \boldsymbol{I}\partial_\theta + \boldsymbol{e}_\theta \otimes \boldsymbol{e}_r - \boldsymbol{e}_r \otimes \boldsymbol{e}_\theta `$ with $`\boldsymbol{I} `$ the three-dimensional identity tensor.
> [!NOTE]
> The convention used in equation (2), that reads in Cartesian coordinates $` \boldsymbol{\nabla u} = u_{i,j} \, \boldsymbol{e}_j \otimes \boldsymbol{e}_i `$ matches our implementation of the gradient. In many textbooks, the gradient of the vector field is defined as $`\boldsymbol{\nabla u} = u_{j,i} \, \boldsymbol{e}_j \otimes \boldsymbol{e}_i `$.

### Tube
We consider a tube with mean radius $R$ and thickness $h$.
We work in cylindrical coordinates $(x,r,\theta)$ and the propagative displacement ansatz reads,
```math
    \boldsymbol{u} = \boldsymbol{u}(k_x,r,m,\omega)e^{i(k_x x + m\theta - \omega t)}, \tag{3}
```
evidencing that a one-dimensional discretization, in the radial direction, allows us to compute axial guided waves in a tube. 

Inserting $(3)$ in $(1)$, we get,
```math
\begin{split}  
    \left[ (ik_x)^2 \boldsymbol{c}_{xx}  + ik_x \left( (\boldsymbol{c}_{xr} + \boldsymbol{c}_{rx}) \partial_r + \frac{1}{r} \left(\boldsymbol{c}_{rx} + \boldsymbol{c}_{x\theta} \cdot \boldsymbol{A} + \boldsymbol{A}\cdot\boldsymbol{c}_{\theta x}\right)\right)\right. + \\
    \left. \boldsymbol{c}_{rr} \partial^2_r + \frac{1}{r} \left(\boldsymbol{c}_{rr} + \boldsymbol{c}_{r\theta} \cdot \boldsymbol{A} + \boldsymbol{A} \cdot \boldsymbol{c}_{\theta r} \right) \partial_r + \frac{1}{r^2} \boldsymbol{A} \cdot \boldsymbol{c}_{\theta\theta} \cdot \boldsymbol{A}  + \rho \omega^2 \boldsymbol{I}\right] \cdot \boldsymbol{u} &= 0,
\end{split} 
\tag{4}
```
where $`\boldsymbol{c}_{ij} = \boldsymbol{e}_i \cdot \boldsymbol{C} \cdot \boldsymbol{e}_j `$ and $`\boldsymbol{A} `$ now reads $` \boldsymbol{A} = i m \boldsymbol{I} + \boldsymbol{e}_\theta \otimes \boldsymbol{e}_r - \boldsymbol{e}_r \otimes \boldsymbol{e}_\theta `$.

The elastodynamics equation is complemented by boundary conditions at $r = R - h/2$ and $r = R + h/2$. Here, we describe the case of a tube in vacuum where we impose homogeneous Neumann boundary conditions,
```math
    \left.\boldsymbol{e}_r \cdot \boldsymbol{C} : \boldsymbol{\nabla u}\right|_{r = R \pm h/2} = \left[ik_x \boldsymbol{c}_{rx} + \boldsymbol{c}_{rr} \partial_r + \frac{1}{r} \boldsymbol{c}_{r\theta} \cdot \boldsymbol{A}\right] \cdot \boldsymbol{u} = 0. \tag{5}
```

The boundary value problem defined by equations (4) and (5) is now discretized using the spectral collocation method.
The displacement vector is approximated at the $N$ Chebyshev-Gauss-Lobatto collocation points $r_i, i \in \{1\ldots N\}$ located in the interval $[R-h/2, R+h/2]$, becoming a vector of size $3N$, and the first and second order differentiation matrices along the radial direction, $D_r$ and $D_{rr}$ of size $N \times N$, are computed using `chebdif.jl`[^3].

The discrete form of equations (4) and (5) is,
<!--
%\begin{align}
%    \begin{split}
%     \left[ (ik_x)^2 \boldsymbol{c}_{xx} \otimes I_N + ik_x \left( (\boldsymbol{c}_{xr} + \boldsymbol{c}_{rx}) \otimes D_2 + \left(\boldsymbol{c}_{rx} + \boldsymbol{c}_{x\theta} \cdot \boldsymbol{A} + \boldsymbol{A}\cdot\boldsymbol{c}_{\theta x}\right) \otimes R_N \right)  \right. + \\
%     \left. \boldsymbol{c}_{rr} \otimes D_{22} + \left(\boldsymbol{c}_{rr} + \boldsymbol{c}_{r\theta} \cdot \boldsymbol{A} + \boldsymbol{A} \cdot \boldsymbol{c}_{\theta r} \right) \otimes R_ND_2 + \boldsymbol{A} \cdot \boldsymbol{c}_{\theta\theta} \cdot \boldsymbol{A} \otimes R_N^2  + \rho \omega^2 \boldsymbol{I} \otimes I_N \right] u &= 0 ,
%     \end{split} \quad \mathrm{for} \; r \in [R-h/2,R+h/2] \\
%    \left[ik_x \boldsymbol{c}_{rx} \otimes I_N + \boldsymbol{c}_{rr} \otimes D_2 + \boldsymbol{c}_{r\theta} \cdot \boldsymbol{A} \otimes R_N \right] u = 0 \quad \mathrm{for} \; r = R \pm h/2
%\end{align}
-->

```math
\begin{aligned}
    \left[(ik_x)^2L_{kk} + ik_x L_k + L_0 + \omega^2 M\right]u &= 0, \quad \mathrm{for} \; r \in [R-h/2,R+h/2]\\ 
    \left[ik_x B_k + B_0 \right]u &= 0, \quad \mathrm{for} \; r = R \pm h/2 
\end{aligned}
\tag{6}
```
with,
```math
\begin{aligned}
    &L_{kk} = \boldsymbol{c}_{xx} \otimes I_N, \quad L_k = (\boldsymbol{c}_{xr} + \boldsymbol{c}_{rx}) \otimes D_r + \left(\boldsymbol{c}_{rx} + \boldsymbol{c}_{x\theta} \cdot \boldsymbol{A} + \boldsymbol{A}\cdot\boldsymbol{c}_{\theta x}\right) \otimes R_N^{-1}, \\
    &L_0 = \boldsymbol{c}_{rr} \otimes D_{rr} + \left(\boldsymbol{c}_{rr} + \boldsymbol{c}_{r\theta} \cdot \boldsymbol{A} + \boldsymbol{A} \cdot \boldsymbol{c}_{\theta r} \right) \otimes R_N^{-1}D_r + \boldsymbol{A} \cdot \boldsymbol{c}_{\theta\theta} \cdot \boldsymbol{A} \otimes R_N^{-2}, \quad M = \rho \boldsymbol{I} \otimes I_N, \\
    &B_k = \boldsymbol{c}_{rx} \otimes I_N, \quad B_0 = \boldsymbol{c}_{rr} \otimes D_r + \boldsymbol{c}_{r\theta} \cdot \boldsymbol{A} \otimes R_N^{-1},
\end{aligned}
```
where $\otimes$ denotes the Kronecker product between two matrices, $I_N$ is the identity matrix of size $N \times N$, $R_N^{-1}$ is the matrix $\mathrm{Diagonal}(1/r_i, i =1 \ldots N)$, and $u$ is the discretized displacement vector.

Finally, we obtain a polynomial eigenvalue problem of size $3N \times 3N$ by replacing the rows of the elastodynamics equation in (6) corresponding to boundary collocation points (*i.e.* the rows $1$, $N$, $N+1$, $2N$, $2N+1$, and $3N$) by the corresponding rows of the boundary condition.
We solve this problem for the eigenpair $(k_x,u)$ parametrized in $\omega$ using the `polyeig` function.

We test this implementation against predictions obtained using the commercial software DISPERSE showing excellent agreement (see the [/Tube/Validation/](/Tube/Validation/) folder).

### Tube with inner fluid
We consider that the elastic tube is filled with fluid. As the domain $r \in [0, R - h/2]$ is finite, we choose to discretize it.
> [!NOTE]
> We model the fluid as an elastic material with a shear modulus $\mu_f \rightarrow 0$, allowing us to reuse the methods described in the previous section.

The only difference lies in the choice of the collocation points. Indeed, the Chebyshev-Gauss-Lobatto points include the domain boundaries, leading to a singularity at $r = 0$ which makes them unsuitable to discretize the fluid domain. Instead, we choose to use $P$ Legendre-Gauss-Radau collocation points which lie in the $]0, R - h/2]$ interval. These collocation points are computed using the function `lgrpointsleft` and the differentiation matrices $Df_r$ and $Df_{rr}$ (with size $P \times P$) using the function `lgrdiffleft`[^4].

We couple the two domains at $r = R - h/2$ by expressing the continuity of the normal traction and of the normal displacement,
```math
\begin{aligned}
    \boldsymbol{e}_r \cdot \boldsymbol{C}_s : \boldsymbol{\nabla} \boldsymbol{u}_s - \boldsymbol{e}_r \cdot \boldsymbol{C}_f : \boldsymbol{\nabla} \boldsymbol{u}_f = 0, \\
    \boldsymbol{e}_r \cdot \boldsymbol{u}_s - \boldsymbol{e}_r \cdot \boldsymbol{u}_f = 0, 
\end{aligned}
\tag{7}
```
where the subscripts $s$ and $f$ stand for solid and fluid. 

Practically, this amounts to constructing two polynomial eigenvalue problems, one for the solid with size $3N \times 3N$ and one for the fluid with size $3P \times 3P$, as described in the section above.
The matrices $L^s_{kk}$, $L^s_{k}$, $L^s_{0}$, $M^s$, and $L^f_{kk}$, $L^f_{k}$, $L^f_{0}$, $M^f$ (for the solid and fluid problem, respectively) are then assembled into four matrices of size $3(N+P)$ so that the elastodynamics equation reads,
```math
\begin{aligned}
    \left[(ik_x)^2\begin{pmatrix}L^s_{kk} & 0_{3N\times 3P} \\ 0_{3P\times 3N} & L^f_{kk} \end{pmatrix} + (ik_x)\begin{pmatrix}L^s_{k} & 0 \\ 0 & L^f_{k}\end{pmatrix} + 
    \begin{pmatrix}L^s_{0} & 0 \\ 0 & L^f_{0}\end{pmatrix}  + \omega^2\begin{pmatrix}M^s & 0 \\ 0 & M^f\end{pmatrix}\right] \begin{pmatrix} u_s \\ u_f \end{pmatrix} = 0. 
\end{aligned}
\tag{8}
```
The Neumann boundary condition at $r = R + h/2$ is then incorporated by replacing the rows $1$, $N+1$, and $2N+1$ of equation (8) by the corresponding rows of of the boundary condition in (6) taking care of adding $3P$ trailing zeros to accommodate the size increase of the coupled displacement vector. 

The coupling boundary conditions given by the equations (7) are discretized similarly as above and injected in the rows $N$, $2N$, $3N$, $3N+1$, $3N+P+1$, and $3N+2P+1$ of the assembled matrices.
 
Finally, we obtain a polynomial eigenvalue problem of size $3(N+P)$ which can be solved for the eigenpair $(k_x,u)$ parametrized by $\omega$.
### Tube with outer fluid
We consider an elastic tube surrounded by an infinite fluid medium. This time it is not convenient to discretize the fluid domain, and we choose to take into account its presence in the boundary condition at $r = R + h/2$ which now reads:
```math
    \boldsymbol{e}_r\cdot(\boldsymbol{C} : \boldsymbol{\nabla} \boldsymbol{u})|_{r=R+h/2} =  \frac{\rho_f\omega^2}{\gamma}\frac{H_m\left(\gamma r\right)}{H'_m\left(\gamma r\right)}  u_r \boldsymbol{e}_r, \tag{9}
```
where $\rho_f$ is the fluid density, $H_m$ is the Hankel function of the first kind of order $m$, and $\gamma = \sqrt{k_f^2 - k_x^2}$  with $k_f$ the wavenumber in the fluid.
Directly discretizing equation (9) and incorporating it in equation (6) transforms the polynomial eigenvalue problem into a non-linear eigenvalue problem.

As solving such non-linear eigenvalue problems is challenging, we turn to the iterative approach proposed by Gravenkamp *et al.*[^5].
We first use the following approximate boundary condition,
```math
    \boldsymbol{e}_r\cdot(\boldsymbol{C} : \boldsymbol{\nabla} \boldsymbol{u})|_{r=R+h/2} = -\frac{\rho_f\omega^2}{k_x} u_r \boldsymbol{e}_r, \tag{10}
```
which can be discretized as,
```math
    \left[(ik_x)^2 \boldsymbol{c}_{rx} \otimes I_N + ik_x\left(\boldsymbol{c}_{rr} \otimes D_r + \boldsymbol{c}_{r\theta} \cdot \boldsymbol{A} \otimes R_N^{-1} \right)+ \rho_f \omega^2 \boldsymbol{I}e_r^T \otimes I_N \right] u = 0, \tag{11}
```
where $e_r$ is the vector $(0,1,0)$, and included in the relevant rows of equation (6).

We can solve the eigenvalue problem,
```math
\begin{split}
    \left[(ik_x)^2L_{kk} + ik_x L_k + L_0 + \omega^2 M\right]u &= 0, \quad \mathrm{for} \; r \in [R-h/2,R+h/2]\\
    \left[(ik_x)^2 B_k + ikB_0 + \omega^2 B_M\right]u &= 0, \quad \mathrm{for} \; r = R + h/2 \\
    \left[ik_x B_k + B_0 \right]u &= 0, \quad \mathrm{for} \; r = R - h/2 
\end{split}
\tag{12}
```
with $B_M = i\rho_f \boldsymbol{I}e_r^T \otimes I_N$, and obtain approximate solutions that we denote $(k_x^+, u^+)$.

> [!NOTE]
> For clarity, we treat here the case of a Neumann boundary condition at the tube inner wall. This approach is also valid for other boundary conditions, as in `Tube_elastic_fluid_in_out.jl` where we consider both inner and surrounding fluid.

We now consider the exact eigenvalue problem, 
```math
    \left[(ik_x)^2L_{kk} + ik_x L_k + L_0 + \omega^2 M(k_x)\right]u = 0, \tag{13}
```
where the appropriate rows of the matrices have been replaced to accommodate the Neumann boundary condition at $r=R-h/2$ and the exact boundary condition given by equation (9) at $r=R+h/2$. We highlight the fact that the eigenvalue problem becomes non-linear by making the dependence of $M$ on $k_x$ explicit.

We perform a Taylor expansion of this eigenvalue problem in the vicinity of $k_x^+$,
```math
\begin{aligned}
    (ik_x)^2L_{kk} + ik_x L_k + L_0 + \omega^2 M(k_x) &=  ik_x\left(2ik_x^+L_{kk}+L_k\right) - \left((ik_x^+)^2L_{kk} - L_0 - \omega^2M(k_x^+)\right)  \\
    & = ik_x A_2(k_x^+) - A_1(k_x^+), 
\end{aligned}
\tag{14}
```
which we inject in equation (10) to obtain the following generalized eigenvalue problem,
```math
    ik_xA_2(k_x^+) u = A_1(k_x^+) u. \tag{15}
```
Inverse iteration is then used to improve the approximate eigenpairs $(k_x^+,u^+)$, starting from their initial value $k_x^{+(0)}$ and $u^{+(0)}$ by computing $u^{+(1)}$
```math
    \left(A_1(k_x^{+(0)}) - ik_x^{+(0)}A_2(k_x^{+(0)}) \right)u^{+(1)} = A_2 u^{+(0)}, \tag{16}
```
and $k_x^{+(1)}$,
```math
    ik_x^{+(1)}A_2(k_x^{+(0)})u^{+(1)} = A_1(k_x^{+(0)})u^{+(1)}, \tag{17}
```
using the backslash operator to invert the matrices in equations (16) and (17).
This procedure is iterated until convergence is reached. 

### Curved strip
We now focus on describing the propagation of axial waves in curved strip with radius $R$, thickness $h$, and curvilinear width $w$. 
We work in cylindrical coordinates $(x,r,\theta)$ and the propagative displacement ansatz reads,
```math
    \boldsymbol{u} = \boldsymbol{u}(k_x,r,\theta,\omega)e^{i(k_x x - \omega t)}, \tag{18}
```
indicating that we perform a two-dimensional discretization, unlike in the case of a tube.

Inserting this ansatz in equation (1) still yields equation (4),
```math
    \begin{split}
     \left[ (ik_x)^2 \boldsymbol{c}_{xx}  + ik_x \left( (\boldsymbol{c}_{xr} + \boldsymbol{c}_{rx}) \partial_r + \frac{1}{r} \left(\boldsymbol{c}_{rx} + \boldsymbol{c}_{x\theta} \cdot \boldsymbol{A} + \boldsymbol{A}\cdot\boldsymbol{c}_{\theta x}\right)\right)\right. + \\
     \left. \boldsymbol{c}_{rr} \partial^2_r + \frac{1}{r} \left(\boldsymbol{c}_{rr} + \boldsymbol{c}_{r\theta} \cdot \boldsymbol{A} + \boldsymbol{A} \cdot \boldsymbol{c}_{\theta r} \right) \partial_r + \frac{1}{r^2} \boldsymbol{A} \cdot \boldsymbol{c}_{\theta\theta} \cdot \boldsymbol{A}  + \rho \omega^2 \boldsymbol{I}\right] \cdot \boldsymbol{u} &= 0,
     \end{split} 
```
but this time the matrix $\boldsymbol{A}$ still involves an azimuthal derivative as $`\boldsymbol{A} = \boldsymbol{I}\partial_\theta + \boldsymbol{e}_\theta \otimes \boldsymbol{e}_r - \boldsymbol{e}_r \otimes \boldsymbol{e}_\theta `$.

Boundary conditions are required at $r = R \pm h/2$ and at $\theta = \pm w/(2R)$. 
We show here the case of a strip clamped on its lateral edges in vacuum,
```math
\begin{aligned}
    \left.\boldsymbol{e}_r \cdot \boldsymbol{C} : \boldsymbol{\nabla u}\right|_{r = R \pm h/2} = \left[ik_x \boldsymbol{c}_{rx} + \boldsymbol{c}_{rr} \partial_r + \frac{1}{r} \boldsymbol{c}_{r\theta} \cdot \boldsymbol{A}\right] \cdot \boldsymbol{u} &= 0, \\
    \left.\boldsymbol{u}\right|_{\theta = \pm w/(2R)} = \boldsymbol{I} \cdot \boldsymbol{u} &= 0. 
\end{aligned}
\tag{19}
```
We discretize this boundary value problem on a grid of $N \times P$ Chebyshev-Gauss-Lobatto collocation points and calculate the radial $N \times N$ differentiation matrices $D_r$ and $D_{rr}$ as well as the azimuthal $P \times P$ differentiation matrices $D_\theta$ and $D_{\theta\theta}$ using `chebdif`. 
The discrete form of the boundary value problem given by equations (4) and (19) reads,
```math
\begin{split}
    \left[(ik_x)^2L_{kk} + ik_x L_k + L_0 + \omega^2 M\right]u &= 0, \quad \mathrm{for} \; r \in [R-h/2,R+h/2], \;\theta \in [-w/(2R), w/(2R)]\\
    \left[(ik_x) B_k + B_0 \right]u &= 0, \quad \mathrm{for} \; r = R \pm h/2 \\
    \left[\boldsymbol{I} \otimes I_{NP} \right]u &= 0, \quad \mathrm{for} \; \theta = \pm w/(2R), 
\end{split}
\tag{20}
```
with,
```math
\begin{aligned}
    &L_{kk} = \boldsymbol{c}_{xx} \otimes I_{NP}, \quad L_k = (\boldsymbol{c}_{xr} + \boldsymbol{c}_{rx}) \otimes I_P \otimes D_r + \left(\boldsymbol{c}_{rx} \otimes I_P + (\boldsymbol{c}_{x\theta} + \boldsymbol{c}_{\theta x}) \otimes D_\theta + (\boldsymbol{c}_{x\theta}A + A\boldsymbol{c}_{\theta x}) \otimes I_P \right) \otimes R_N^{-1}, \\
    &L_0 = \boldsymbol{c}_{rr} \otimes I_P \otimes D_{rr} + \left(\boldsymbol{c}_{rr} \otimes I_P + \boldsymbol{c}_{\theta r} \otimes D_\theta + A\boldsymbol{c}_{\theta r} \otimes I_P \right) \otimes R_N^{-1}D_r +  (\boldsymbol{c}_{r\theta} \otimes D_\theta + \boldsymbol{c}_{r\theta}A \otimes I_P ) \otimes (D_rR_N^{-1} + R_N^{-2}) + \\
    &\quad \left(\boldsymbol{c}_{\theta\theta} \otimes D_{\theta\theta}+ (\boldsymbol{c}_{\theta\theta}A + A\boldsymbol{c}_{\theta\theta}) \otimes D_\theta + A\boldsymbol{c}_{\theta\theta}A \otimes I_P\right) \otimes R_N^{-2}, \quad M = \rho \boldsymbol{I} \otimes I_{NP}, \\
    &B_k = \boldsymbol{c}_{rx} \otimes I_{NP}, \quad B_0 = \boldsymbol{c}_{rr} \otimes I_P \otimes D_r + (\boldsymbol{c}_{r\theta} \otimes D_\theta + \boldsymbol{c}_{r\theta}A \otimes I_P) \otimes R_N^{-1},
\end{aligned}
```
where $` A = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & -1 \\ 0 & 1 & 0 \end{pmatrix} `$,

and all the matrices of the discretized polynomial eigenvalue problem are of size $3NP \times 3NP$.

>[!NOTE]
> In the implementation of the $L_0$ term in equation (20), we found that we have to write $\frac{\partial}{\partial r}\left(\frac{1}{r}-\right) + \frac{1}{r^2}-$ instead of $\frac{1}{r}\frac{\partial}{\partial r}-$.

We validate this implementation by comparing predictions for *almost* flat strips (*i.e.* with $R\gg w$) to Comsol simulations of the in-plane guided modes for free and clamped flat strips in [/Strip/Validation/](/Tube/Validation/). We also check the influence of $R$ against Comsol simulations, focusing on the first out-of-plane mode in [/Strip/Validation/Courbure](/Tube/Validation/Courbure), obtaining excellent agreement in both test cases.

> [!WARNING]
> The mixed boundary conditions used in the case of a clamped strip lead to issues at the corners of the domain. Use an even number number of collocation points to compute the in-plane waves and an odd number of collocation points to compute the out-of-plane wave dispersion relation.

### Curved strip with approximate coupling to a semi-infinite fluid domain
We discuss the coupling of a curved strip to a semi-infinite fluid medium. 
We introduce an approximate boundary condition that reads,
```math
    \left.\boldsymbol{e}_r \cdot \boldsymbol{C} : \boldsymbol{\nabla u}\right|_{r = R + h/2} = \rho_f \omega^2 \frac{\boldsymbol{u}\cdot \boldsymbol{e}_r}{i\boldsymbol{k}_f\cdot\boldsymbol{e}_r}, \tag{21}
```
where we used the continuity of the normal displacement at the interface, and where $`\boldsymbol{k}_f `$ is the wavevector of an inhomogeneous plane wave radiated into the fluid. 
We further simplify this boundary condition by considering that only evanescent waves propagate in the fluid, and by neglecting the transverse wavenumber in the $` \boldsymbol{e}_\theta `$ direction. Under these assumptions, the approximate boundary condition at $r=R+h/2$ can be discretized as,
```math
    \left[(ik_x)^2 B_k + ik_x B_0 + i\rho_f \omega^2 \boldsymbol{I}e_r^T \otimes I_{NP} \right] \cdot \boldsymbol{u} = 0. \tag{22}
```
After inserting the rows of equation (22) in the rows corresponding to collocation points at the top boundary in the discretized eigenvalue problem defined in (20), we can determine the eigenpairs $(k_x,u)$ corresponding to the problem of a curved strip coupled to a semi-infinite fluid medium.
### Viscoacoustoelasticity
We also provide scripts to compute the dispersion relation of small amplitude waves in homogeneously pre-stressed viscoelastic media. 
In the framework of acoustoelatic theory, the propagation of incremental waves on top of a finitely deformed configuration is still described by the wave equation given in (1) where the elasticity tensor $\boldsymbol{C}$ is replaced by a modified elasticity tensor $\boldsymbol{C^0}$, which takes into account the intertwined contributions of pre-stress and viscoleasticity, and where $\boldsymbol{u}$ now stands for the incremental displacement.

This similarity allows us to use the methods described for a stress-free curved strip after:

- computing the modified elasticity tensor $\boldsymbol{C^0}$ using the function `C_viscoacoustoelastic` which takes as input the stretch ratios in the $x$ and $\theta$ direction, the frequency $f$ and the relevant material parameters used in the compressible Mooney-Rivlin model and in the fractional Kelvin-Voigt model. We refer the reader to the repository [GEW soft strip](https://github.com/dakiefer/GEW_soft_strip) for details[^2].
- adapting the geometry to the considered level of pre-stress as $w = \lambda_\theta w_0$, and $h = h_0/(\lambda_\theta \lambda_x)$.

## Dependencies
The functions `chebdif.jl` and `chebpoints.jl` are bundled from [DMSuite.jl](https://github.com/l90lpa/DMSuite.jl) with their license file.

The functions `lpoly.jl`, `lgrpointsleft.jl`, and `lgrdiffleft.jl` are adapted from the book:
> J. Shen, T. Tang and L. Wang, Spectral Methods: Algorithms, Analysis and Applications, Springer Science & Business Media (2011) [https://doi.org/10.1007/978-3-540-71041-7](https://doi.org/10.1007/978-3-540-71041-7)

The package [NonlinearEigenProblems.jl](https://github.com/nep-pack/NonlinearEigenproblems.jl) is developed by E. Jarlebring, M. Bennedich, G. Mele, E. Ringh and  P. Upadhyaya:
> E. Jarlebring, M. Bennedich, G. Mele, E. Ringh and  P. Upadhyaya, NEP-PACK: A Julia package for nonlinear eigenproblems, *[arXiv:1811.09592](https://arxiv.org/abs/1811.09592)* (2018)

## Authors
Code created in 2024-2025 by

P. Chantelot, Institut Langevin, ESPCI Paris, Université PSL - [Google Scholar](https://scholar.google.fr/citations?user=BQWXUKYAAAAJ&hl=fr&oi=ao) 

D. A. Kiefer, Institut Langevin, ESPCI Paris, Université PSL - [Google Scholar](https://scholar.google.fr/citations?hl=fr&user=odSy3v4AAAAJ)

[^1]: D. A. Kiefer, _Elastodynamic quasi-guided waves for transit-time ultrasonic flow metering_, ser. FAU Forschungen, Reihe B, Medizin, Naturwissenschaft, Technik, vol. 42. Erlangen: FAU University Press (2022), doi: [10.25593/978-3-96147-550-6](http://doi.org/10.25593/978-3-96147-550-6).
[^2]: A. Delory, D. A. Kiefer, M. Lanoy, A. Eddi, C. Prada, and F. Lemoult, Viscoelastic dynamics of a soft strip subject to a large deformation, *Soft Matter*, (2024), doi: [10.1039/D3SM01485A](https://doi.org/10.1039/D3SM01485A).
[^3]: J. A. Weideman and S. C. Reddy, “A MATLAB Differentiation Matrix Suite,” ACM Trans. Math. Softw., vol. 26, no. 4, pp. 465–519 (2000), doi: [10.1145/365723.365727](http://doi.org/10.1145/365723.365727).
[^4]: J. Shen, T. Tang and L. Wang, Spectral Methods: Algorithms, Analysis and Applications, Springer Science & Business Media (2011) [10.1007/978-3-540-71041-7](https://doi.org/10.1007/978-3-540-71041-7)
[^5]: H. Gravenkamp, C. Birk, and C. Song, Numerical modeling of elastic waveguides coupled to infinite fluid media using exact boundary conditions, Comput. Struct. 141,36 (2014) [10.1016/j.compstruc.2014.05.010](https://doi.org/10.1016/j.compstruc.2014.05.010) 
