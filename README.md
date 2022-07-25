# Simula_SummerProject


## To do list
- [ ] Consider closing outfile and appending during simulation in order to save the data if it crashes
- [ ] Make struct for PDE parameters for printing of function body


## Introduction 


## Theory 
### Equations 

$u_S, p_S$ is velocity and pressure in stokes domain $S$, and $p_D$ is pressure in Darcy domain $D$.

#### Stokes

$$
\begin{align}
    - \nabla \cdot \sigma(u_S, p_S) &= f_s \quad \text{in} \ \Omega_S \\
    \nabla \cdot u_S &= 0  \quad \text{in} \ \Omega_S \\
    u_S &= u_{S,0}  \quad \text{on} \ \Lambda_S \\
    u_S \times \hat{n}_S &= 0  \quad \text{on} \ \Gamma_S \\
    p_S &= p_{S,0}  \quad \text{on} \ \Gamma_S \\
\end{align}
$$

where

$$
\begin{align}
   \sigma(u_S, p_S) &= 2\mu \varepsilon(u_S) - p_S \mathbb{I} \\
   \varepsilon(u_S) &= \frac{1}{2}(\nabla u_S + \nabla^Tu_S)
\end{align}
$$

#### Darcy 

$$
\begin{align}
    \nabla \cdot (-\frac{\kappa}{\mu}\nabla p_D) &= f_D \quad \text{in} \ \Omega_D \\
     P_D &= P_{D,0}  \quad \text{on} \ \partial \Omega_S \setminus \Gamma 
\end{align}
$$


#### Interface 
$$
\begin{align}
    u_S\cdot\hat{n}_S + (-\frac{\kappa}{\mu}\nabla p_D \cdot \hat{n}_D) &= g\Gamma \\
    -[\sigma(u_S, p_S)\cdot\hat{n}_S]\cdot\hat{n}_S &= P_D \\
    -[\sigma(u_S, p_S)\cdot\hat{n}_S]\cdot\hat{\tau}_S &= \alpha u_S \cdot\hat{\tau}_S\\
\end{align}
$$
___

### Conditons for the brain simulations

We are going to use 
$$
\begin{align}
    u_{S,0} &= 0 \\

\end{align}
$$
___

### Weak formulation

#### Stokes

$
\begin{align}
    \int_{\Omega_S} f_s \cdot v_S =  \int_{\Omega_S} 2\mu \varepsilon(u_S) \odot \varepsilon(v_S) \ dx - \int_{\Omega_S} p_S \nabla\cdot v_S \ dx + \int_{\partial\Omega_S} [-\sigma(u_S, p_S)\cdot\hat{n}_S] \cdot v_S \ dS  \ldots
\end{align}
$

#### Darcy 



## Results and Discussion

![Alt text](/figures/solution_convergence.png "Optional title")

<img src="/figures/solution_convergence.png" alt="Alt text" title="Optional title">
