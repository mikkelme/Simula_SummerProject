# Simula_SummerProject (FIND better name)

This repo contains the work done as a summer intern at Simula during a six week period during the summer of 2022. The project was guided by my supervisor Miroslav Kutcha (Github link).

## Introduction 

We are going to build a framwork for simulating the flow of the cerebrospinal which flows on the outside of the brain. We want to investigate whether one can use dimension reduction to increase the efficiency of such simulations. That is, modelling the outer shell in 3D as a surface without thickness in 2D. Or in the simpler case we can consider going from the 2D model, where the fluid flows across a 2D surface and investigating whether this can be redyced to a 1D line. In this project we are goind to limit ourself to the dimension reduction starting from the 2D case. However, the framework for making the 3D brain geometry is already included here such that an extension to the 3D study is not to demanding.

We are going to use Julia as the programming language for this project. We work with Gridap for the finite element part and gmsh for creating the geometry and the mesh.

## Method & Theory 

### Domain

We are going to model our brain as a composition of two domains: The *Stokes* domain and the *Darcy* domain corresponding the equations that governs the fluid flow in these domains. In the *Stokes* domain the fluid flows in a free path on the outside of the brain tissude, where the motion is described as Stokes flow (low Reynolds number). In the *Darcy* domain the fluid flows though the pores of the brain tissue where the fluid motion is described as percolation. 

<!-- Insert model image -->

### Equations 

We denote $u_S, p_S$ as velocity and pressure respectively in the Stokes domain $S$, and $p_D$ as pressure in the Darcy domain $D$. We define the problem by the following equations.

#### Stokes domain

$$
\begin{align}
    - \nabla \cdot \sigma(u_S, p_S) &= f_s  &\text{in} \ \Omega_S \\
    \nabla \cdot u_S &= 0   &\text{in} \ \Omega_S \\
    u_S &= u_{S,0}   &\text{on} \ \Lambda_S \\
    u_S \times \hat{n}_S &= 0   &\text{on} \ \Gamma_S \\
    p_S &= p_{S,0}   &\text{on} \ \Gamma_S \\
\end{align}
$$

where

$$
\begin{align}
   \sigma(u_S, p_S) &= 2\mu \varepsilon(u_S) - p_S \mathbb{I} \\
   \varepsilon(u_S) &= \frac{1}{2}(\nabla u_S + \nabla^Tu_S)
\end{align}
$$

#### Darcy domain

$$
\begin{align}
    \nabla \cdot (-\frac{\kappa}{\mu}\nabla p_D) &= f_D  &\text{in} \ \Omega_D \\
     P_D &= P_{D,0}  &\text{on} \ \partial \Omega_S \setminus \Gamma 
\end{align}
$$



#### Interface conditions
$$
\begin{align}
    u_S\cdot\hat{n}_S + (-\frac{\kappa}{\mu}\nabla p_D \cdot \hat{n}_D) &= g\Gamma  &\text{on} \ \Gamma \\ 
    -[\sigma(u_S, p_S)\cdot\hat{n}_S]\cdot\hat{n}_S &= P_D  &\text{on} \ \Gamma \\
    -[\sigma(u_S, p_S)\cdot\hat{n}_S]\cdot\hat{\tau}_S &= \alpha u_S \cdot\hat{\tau}_S  &\text{on} \ \Gamma \\
\end{align}
$$


### Weak formulation

Check consistency with $dx$ or $dV, dS, dL$.

#### Stokes
From equation (...) we get

$$
\begin{align}
    \int_{\Omega_S} f_s \cdot v_S \ dx &=  \int_{\Omega_S} 2\mu \ \varepsilon(u_S) \odot \varepsilon(v_S) \ dx - \int_{\Omega_S} p_S \nabla\cdot v_S \ dx - \int_{\partial\Omega_S} \big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot v_S \ dx \\ 
    0 &= - \int_{\Omega_S} (\nabla \cdot u_S) \cdot q_s
\end{align}
$$

In the first equation (put number) Decompose the last term in normal $\hat{n}_S$ and tangential $\hat{\tau}_S$ direction

$$
\begin{align}
    - \int_{\partial\Omega_S} \big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot v_S \ dS   &=  \int_{\partial\Omega_S} \underbrace{-\Big[\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{n}_S \Big]\Big[\hat{n}_S \cdot v_S \Big]}_{P_D \ \text{on} \ \ \Gamma} -\underbrace{\Big[\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{\tau}_S \Big]\Big[\hat{\tau}_S \cdot v_S \Big]}_{\alpha u_S \cdot\hat{\tau}_S \ \ \text{on} \ \Gamma} \\
    &= \int_{\partial\Omega_S\setminus\Gamma} -\Big[\underbrace{\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{n}_S }_{- p_S} \Big]\Big[\hat{n}_S \cdot v_S \Big] -\Big[\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{\tau}_S \Big]\Big[\hat{\tau}_S \cdot v_S \Big]  +  \int_{\Gamma} P_D - \alpha u_S \cdot\hat{\tau}_S \ dL
\end{align}
$$

where we used the transistion


$$
\begin{align}
    \big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{n}_S &= \big(2\mu \ \varepsilon(u_S) - p_S \mathbb{I} \big)\cdot\hat{n}_S \\
    &= 2\mu \ \underbrace{\hat{n}_S \cdot \varepsilon(u_S) \cdot \hat{n}_S}_{0} - p_S \ \underbrace{\hat{n}_S \cdot \hat{n}_S}_{1} = - p_S
\end{align}
$$

with the zero in line (...) comes from theorem something in that paper (REFER).

We can then handle the remaining tangential component of the  $\partial\Omega_S\setminus\Gamma$ boundary using the Nitsche method (theorem?).  


$$
\begin{align}
    - \int_{\partial\Omega_S\setminus\Gamma} \Big[\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{\tau}_S \Big]\Big[\hat{\tau}_S \cdot v_S \Big] = 
    &- \int_{\partial\Omega_S\setminus\Gamma}\Big[\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{\tau}_S \Big]\Big[\hat{\tau}_S \cdot v_S \Big] \\
    &- \int_{\partial\Omega_S\setminus\Gamma}\Big[\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{\tau}_S \Big]\Big[\hat{\tau}_S \cdot u_S - u_{S,\text{tan}} \Big] \\
    &+ \int_{\partial\Omega_S\setminus\Gamma}\frac{\gamma}{h} \Big[ \hat{\tau}_S \cdot u_S - u_{S,\text{tan}} \Big]\Big[\hat{\tau}_S \cdot u_S \Big]
\end{align}
$$

where $u_{S,\text{tan}}$ is the condition for the tangential part of the stokes velocity on the $\partial\Omega_S\setminus\Gamma$. We only want the normal component and thus we set $u_{S,\text{tan}}$ = 0



#### Darcy 

$$
\begin{align}
    \int_{\Omega_D} f_D \cdot q_d \ d\vec{x} &= \int_{\Omega_D} \frac{\kappa}{\mu} \nabla p_D \cdot \nabla q_D \ dx + \int_{\partial\Omega_D} \underbrace{-\hat{n}_D \cdot \frac{\kappa}{\mu} \nabla p_D}_{g\Gamma - u_S \cdot \hat{n}_S \ \text{on} \ \Gamma} \cdot q_D \ dx \\
    &=\int_{\Omega_D} \frac{\kappa}{\mu} \nabla p_D \cdot \nabla q_D \ dx + \int_{\Gamma} (g\Gamma - u_S \cdot \hat{n}_S) \cdot q_D - \int_{\partial\Omega\setminus\Gamma} \hat{n}_D \cdot \kappa \nabla p_D \cdot q_D
\end{align}
$$

where we handle the last term as neuman condition.


## Modelling the brain

We are going to use 
$$
\begin{align}
    u_{S,0} &= 0 \\

\end{align}
$$

$\nabla \cdot x$ is the same as $Div(x)$.


## Results and Discussion


<!-- <p align="center">
  <img src="figures/solution_convergence.png" alt="Text is not showing?"/>
</p> -->


