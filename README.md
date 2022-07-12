# Simula_SummerProject


## To do list
- [x] Wiring for directly importing mesh model (https://github.com/gridap/GridapGmsh.jl/blob/master/src/GmshDiscreteModels.jl)
- [ ] Basbuška


## Introduction 


## Theory 

### Stokes 

$$ -\nabla \cdot \sigma(u,p) = f, \quad -\nabla \cdot u = 0, \quad \in \Omega, \quad \partial \Omega = \Gamma_D \cup \Gamma_N $$

$$ \sigma(u,p) = 2\mu \varepsilon(u) - P \mathbb{I}, \quad \varepsilon(u) = \frac{1}{2}\left(\nabla u + (\nabla u)^T \right) $$

$$ u = g \ \text{on} \ \Gamma_D, \quad \sigma \cdot v = h \ \text{on} \ \Gamma_N $$



### Babuska Stokes

Start from 

$$ -2 \mu \nabla \cdot(\mathbf{D}(\mathbf{u}))+\nabla p =\rho \mathbf{f}, \quad \text { in } \Omega $$
$$ \nabla \cdot \mathbf{u} =0, \quad \text { in } \Omega $$
$$ \mathbf{u} =\mathbf{0 ,} \text { on } \Gamma_{1} $$
$$ \mathbf{u} \times \mathbf{n} =\mathbf{0}, \quad \text { on } \Gamma_{2}, \text { and } $$
$$ p =p_{0}, \text { on } \Gamma_{2} $$



