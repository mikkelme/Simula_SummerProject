# Simula_SummerProject


<!-- ## To do list
- [x] Wiring for directly importing mesh model (https://github.com/gridap/GridapGmsh.jl/blob/master/src/GmshDiscreteModels.jl)
- [ ] Basbuška -->


## Introduction 


## Theory 

### Equations 

$u_S, p_S$ is velocity and pressure in stokes domain, and $p_D$ is pressure in Darcy domain.

$$

\begin{align}
    - \nabla \cdot \sigma(u_S, p_S) &= f_s \quad \text{in} \ \Omega_S \\
    \nabla \cdot u_S &= 0  \quad \text{in} \ \Omega_S \\
    u_S &= g  \quad \text{on} \ \Lambda_S
\end{align}
\\
\begin{align}
   \sigma(u_S, p_S) &= 2\mu \varepsilon(u_S) - p_S \mathbb{I} \\
   \varepsilon(u_S) &= \frac{1}{2}(\nabla u_S + \nabla^Tu_S)
\end{align}


$$
---

$$

- \nabla \cdot \sigma(u_S, p_S) = f_s, \quad \text{in} \ \Omega_S \\
x = 4
$$



- Coupled equations + Nitsche extension

- Domain
- Boundary conditions


<!-- ### Stokes 

$$ -\nabla \cdot \sigma(u,p) = f, \quad -\nabla \cdot u = 0, \quad \in \Omega, \quad \partial \Omega = \Gamma_D \cup \Gamma_N $$

$$ \sigma(u,p) = 2\mu \varepsilon(u) - P \mathbb{I}, \quad \varepsilon(u) = \frac{1}{2}\left(\nabla u + (\nabla u)^T \right) $$

$$ u = g \ \text{on} \ \Gamma_D, \quad \sigma \cdot \nu = h \ \text{on} \ \Gamma_N $$



### Babuska Stokes

Start from 

$$ -2 \mu \nabla \cdot(\mathbf{D}(\mathbf{u}))+\nabla p =\rho \mathbf{f}, \quad \text { in } \Omega $$

$$ \nabla \cdot \mathbf{u} =0, \quad \text { in } \Omega $$

$$ \mathbf{u} =\mathbf{0 ,} \text { on } \Gamma_{1} $$

$$ \mathbf{u} \times \mathbf{n} =\mathbf{0}, \quad \text { on } \Gamma_{2}, \text { and } $$

$$ p =p_{0}, \text { on } \Gamma_{2} $$ -->



