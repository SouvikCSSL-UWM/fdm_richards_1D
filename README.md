# fdm_richards_1D
Finite difference solution of 1D Richards' equation (head-based) using Picard iteration. Backward Euler for time, central difference for space, and Crank–Nicolson scheme for diffusion and advection. Includes multiple initial conditions and supports unsaturated soil flow simulations.

This repository provides a Python-based solver for the **1D Richards Equation** in unsaturated soils, employing the **Van Genuchten model** for soil hydraulic properties. The solver uses **finite difference methods** with **Picard iteration** and supports **two types of boundary conditions**: constant head and constant flux.

## 🔬 Governing Equation

The water content form of the 1D Richards equation solved is:

\[
C(h) \frac{\partial h}{\partial t} = \frac{\partial}{\partial z} \left[ K(h) \left( \frac{\partial h}{\partial z} + 1 \right) \right]
\]

Where:

- \( h \): Pressure head [m]  
- \( \theta \): Volumetric water content  
- \( K(h) \): Unsaturated hydraulic conductivity [m/s]  
- \( C(h) \): Specific moisture capacity [1/m]

The water retention and hydraulic conductivity functions follow the **Van Genuchten** model.

## ⚙️ Features

- Nonlinear soil moisture dynamics simulation.
- Van Genuchten parameterization:
  - \( \alpha = 7.5 \)
  - \( n = 1.89 \)
  - \( m = 1 - 1/n \)
  - \( \theta_s = 0.41 \)
  - \( \theta_r = 0.065 \)
  - \( K_s = 1.06 \times 10^{-4} \, \text{m/s} \)
- Time integration using Picard iterations.
- Supports:
  - **Constant Head (Dirichlet)** boundary condition.
  - **Constant Flux (Neumann)** boundary condition.
- Animated plot of θ (water content) over time.

## 📌 Boundary Conditions

You can toggle between two boundary conditions by modifying the `Picard_iteration()` function in the code:

### 1. Constant Head (Dirichlet)

Uncomment the following lines to apply a fixed head boundary:

```python
# A[0,0] = 1
# A[0,1] = 0
# b[0] = h[0]
# A[-1, -1] = 1
