# Dynamic Predictive Control for Networked Control Systems (Delay + Packet Loss)

This repository contains MATLAB code accompanying the paper:

**The design of dynamic predictive control for networked control systems subject to latency and packet loss**  
Ali Safi, Mostafa Nasiri, Reza Farasat

The project implements a **dynamic predictive control** scheme for **networked control systems (NCSs)** operating over a **UDP-like communication channel** subject to:
- **random measurement latency** (time-varying delay), and
- **packet dropout** (data loss without acknowledgment).

The controller is based on an **LQTI (Linear Quadratic Tracker with Integrator)** structure and extends classical networked predictive control by predicting **both**:
1) the plant state evolution, and  
2) the controller internal dynamics,  
to maintain stable regulation and tracking under delay and packet loss.

---

## Repository contents (expected)

Typical files used by the simulation include:

- MATLAB scripts/functions (`*.m`) implementing:
  - plant/controller discretization
  - predictive estimation logic under delay/dropout
  - simulation and plotting
- `StateSpace.mat` (plant matrices)
- `Inverted_Pendulum.m` (nonlinear pendulum dynamics)
- `RungeKutta.m` (Runge–Kutta Solver)

If your repository uses different filenames, update the references above accordingly.

---

## Method summary

### Plant model (discrete time)

$x_p(k + 1) = A_p x_p(k) + B_p u_p(k)$

$y_p(k) = C_p x_p(k)$

where:
- $x_p \in \mathbb{R}^{n}$ is the plant state,
- $u_p \in \mathbb{R}^{m}$ is the plant input,
- $y_p \in \mathbb{R}^{p_p}$ is the plant output.

### Dynamic controller (LQTI form)

$x_c(k + 1) = A_c x_c(k) + B_c u_c(k)$

$y_c(k) = C_c x_c(k) + D_c u_c(k)$

The formulation used in the paper applies the following assumptions:

$y_{ref} = 0$

$u_p(k) = y_c(k)$

$u_c(k) = -\hat{x}_p(k)$

---

## Network model (UDP-like)

The measurement channel is modeled as UDP-like:
- packets experience random delay, and
- packets may be dropped without retransmission/acknowledgment.

A single parameter $M$ is used as a proxy for network quality (as described in the paper):

- Minimum delay: $M$ samples  
- Maximum delay: $2M$ samples  
- Packet loss probability: $10M\%$

This repository is set up to run simulations for:

$M \in \{1,2,3,4\}$

---

## Prediction logic (core equations)

Let $d_\tau$ denote the effective number of consecutive delayed/missing samples used by the predictor.

### Predicted plant state

$\hat{x}_p(k+1) = A_p^{d_\tau+1} x_p(k-d_\tau)$

$x_{p}(k+1) = A_p^{d_\tau+1} x_p(k-d_\tau) + \sum_{n=1}^{d_\tau+1} A_p^{d_\tau+1-n} B_p u_p \left(k-d_\tau+n-1\right)$

### Predicted controller state

$x_c(k+1) = A_c^{d_\tau+1} x_c(k-d_\tau) - \sum_{n=1}^{d_\tau+1} A_c^{d_\tau+1-n} B_c \hat{x}_p \left(k-d_\tau+n-1\right)$

### Control signal

$u_p(k + 1) = C_c x_c(k + 1) - D_c \hat{x}_p(k + 1)$


## Quick start

1. Put all MATLAB files in one folder.
2. Ensure `StateSpace.mat` exists in the same folder and contains:
   - `A` (continuous-time state matrix)
   - `B` (continuous-time input matrix)
3. Run the main simulation script from MATLAB.

If you are using the refactored script that sweeps $M$ automatically, it will run:

$M = 1,2,3,4$

and generate:
- dropout + delay profile plots for each $M$
- response comparison plots (predictive vs. baseline)

## License

This project is licensed under the MIT License. See [`LICENSE`](LICENSE).
