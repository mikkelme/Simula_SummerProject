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
    u_S\cdot\hat{n}_S + (-\frac{\kappa}{\mu}\nabla p_D \cdot \hat{n}_D) &= g\Gamma \quad \text{on} \ \Gamma \\ 
    -[\sigma(u_S, p_S)\cdot\hat{n}_S]\cdot\hat{n}_S &= P_D \quad \text{on} \ \Gamma \\
    -[\sigma(u_S, p_S)\cdot\hat{n}_S]\cdot\hat{\tau}_S &= \alpha u_S \cdot\hat{\tau}_S \quad \text{on} \ \Gamma \\
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

$\nabla \cdot x$ is the same as $Div(x)$.
___

### Weak formulation

#### Stokes
From equation (...) we get

$$
\begin{align}
    \int_{\Omega_S} f_s \cdot v_S &=  \int_{\Omega_S} 2\mu \ \varepsilon(u_S) \odot \varepsilon(v_S) \ dx - \int_{\Omega_S} p_S \nabla\cdot v_S \ dx - \int_{\partial\Omega_S} \big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot v_S \ dS \\ 
    0 &= - \int_{\Omega_S} (\nabla \cdot u_S) \cdot q_s
\end{align}
$$

In the first equation (put number) Decompose the last term in normal $\hat{n}_S$ and tangential $\hat{\tau}_S$ direction

$$
\begin{align}
    - \int_{\partial\Omega_S} \big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot v_S \ dS   &=  \int_{\partial\Omega_S} \underbrace{-\Big[\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{n}_S \Big]\Big[\hat{n}_S \cdot v_S \Big]}_{P_D \ \text{on} \ \ \Gamma} -\underbrace{\Big[\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{\tau}_S \Big]\Big[\hat{\tau}_S \cdot v_S \Big]}_{\alpha u_S \cdot\hat{\tau}_S \ \ \text{on} \ \Gamma} \\
    &= \int_{\partial\Omega_S\setminus\Gamma} -\Big[\underbrace{\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{n}_S }_{- p_S} \Big]\Big[\hat{n}_S \cdot v_S \Big] -\Big[\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{\tau}_S \Big]\Big[\hat{\tau}_S \cdot v_S \Big]  +  \int_{\Gamma} P_D - \alpha u_S \cdot\hat{\tau}_S \ dL
    % &= \int_{\partial\Omega_S\setminus\Gamma} -\Big[\underbrace{\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{n}_S }_{
    %     (2\mu \varepsilon(u_S) - p_S \mathbb{I})\cdot\hat{n}_S

    % } \Big]\Big[\hat{n}_S \cdot v_S \Big] -\Big[\big(\sigma(u_S, p_S)\cdot\hat{n}_S\big) \cdot \hat{\tau}_S \Big]\Big[\hat{\tau}_S \cdot v_S \Big]  +  \int_{\Gamma} P_D - \alpha u_S \cdot\hat{\tau}_S \ dL
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

We can then handle the remaining $\partial\Omega_S\setminus\Gamma$ boundary using the Nitsche method (theorem?).  again 



#### Darcy 



## Results and Discussi


<!-- <p align="center">
  <img src="figures/solution_convergence.png" alt="Text is not showing?"/>
</p> -->


