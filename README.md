# Efficient Uncertainty Quantification for Emergency Well Planning

This repository contains the numerical implementation of the methodology described in the paper:  
**"Efficient approach to uncertainty quantification for emergency well planning under data scarcity"**.

## Overview
We propose an efficient strategy, based on the **Method of Distributions (MoD)**, to evaluate groundwater availability under data scarcity. The methodology utilizes a **Meshless Generalized Finite-Difference (mGFD)** technique to solve deterministic equations for the cumulative distribution function (CDF) of hydraulic head.

### Key Highlights
* **High Efficiency**: Our uncertainty-quantification methodology is orders of magnitude faster than Monte Carlo simulations.
* **Rare Event Estimation**: It enables the estimation of rare events, furnishing valuable guidance for risk assessment.
* **Field Application**: The method is illustrated on a field application dealing with emergency well planning in the Taichung Basin, Taiwan.

## Repository Structure
The repository includes five main MATLAB scripts corresponding to the alternative methods evaluated in the study:

1. **`uniformFlow_MCS.m`**: Monte Carlo simulations (MCS) providing the "ground truth" reference solution.
2. **`uniformFlow_MCSCDF.m`**: Hybrid MC-CDF method using a reduced number of MC realizations.
3. **`uniformFlow_MEsCDF.m`**: Deterministic ME-CDF method solved with the Finite Difference Method (FDM).
4. **`uniformFlow_mMEsCDF.m`**: Deterministic mME-CDF method using the meshless mGFD scheme on regular nodes.
5. **`uniformFlow_iMEsCDF.m`**: Deterministic iME-CDF method using **Advancing Front Node (AFN)** generation for irregular and adaptive node placement.

## Prerequisites
* **MATLAB**: Developed and tested on MATLAB R2024b.
* **Parallel Computing Toolbox**: Recommended for running Monte Carlo simulations.
* **S-GeMS (Stanford Geostatistical Modeling Software)**: Required for the Sequential Gaussian Simulation (SGS) used in MCS scripts.

## Required Helper Functions
Ensure the following supporting functions (available in this repository) are in your MATLAB path:
* `fct_GFDM_Coef`: Computes the meshless derivative weights.
* `fct_2DStencil_FDM` / `fct_CDF3DStencil_FDM`: Generates finite difference stencils.
* `fct_Steady_DeterHead_HLG_FDM`: Solves the steady-state groundwater flow equation.
* `fct_Steady_Closure_FDM`: Calculates closure variables for the MoD equations.
* `fct_Steady_CdfHead`: Solves the final CDF governing equation.
* `node_drop_2d_ctps_tol`: Implementation of the AFN node generation algorithm.

## Citation
If you use this code or methodology in your research, please cite:
> Shang-Ying Chen, Kuo-Chin Hsu, Chien-Chung Ke, Nai-Chin Chen, and Daniel M. Tartakovsky (2026). Efficient approach to uncertainty quantification for emergency well planning under data scarcity. 

## Acknowledgments
This research was supported by the National Science and Technology Council (NSTC) of Taiwan, Sinotech Engineering Consultants, the U.S. Strategic Environmental Research and Development Program (SERDP), and the National Science Foundation (NSF).
