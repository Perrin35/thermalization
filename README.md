# Repository Overview

This repository includes the notebook and data used in [[1]](#references) and can be cited with the DOI: [10.5281/zenodo.14639657](https://doi.org/10.5281/zenodo.14639656). It is organized as follows:

## Directory Structure

### `data/`
Contains all data used in the project, organized into subdirectories:

- **`circuits/`**:
  - Stores quantum circuit files used for simulations.
  - Includes 100 disorder realization folders **`disorder_realization_<i>/`**, each containing:
    - `random_interaction_<i>.npy` (1 ≤ i ≤ 100): NumPy file storing the interaction matrix \( J \), where \( J[k,l] \) represents the coupling strength between qubits \( k \) and \( l \), and \( J[l,k] \) is its conjugate.
    - Subfolders for different simulation times: **`T=2.5/`**, **`T=5/`**, **`T=10/`**:
        - **Recompiled Circuits**: Two circuits `recompiled_circuit_<initial_state>_disorder_<i>_T=<T>.qasm` in QASM format per disorder realization for the two different initial states: ground state and highest energy state.
        - **`RC/`**: Contains 100 crosstalk randomized compiling (cRC, see [[2]](#references)) versions of the recompiled circuits.
        - **`ZNE3/`**: Contains ZNE circuits with tripled CNOT gates and 100 cRC versions of each.

- **`quantum data/`**:
  - Raw data from six IBMQ sessions in subfolders:
    - **`session thermalization <i>/`** (1 ≤ i ≤ 6): Simulations at \( T = 2.5, 5, 10 \) for both initial states. The 3 first sessions correspond to the ground state with session 1 for T=2.5, session 2 for T=5 and session 3 for T=10. The 3 last sessions correspond to the highest energy state with the same ordering for the time.
      - **`measurement_backend=ibm_hanoi.json`**: `measurement_backend=ibm_hanoi.json`: A 16x16 matrix representing the calibration of the measurement for that session. The 
2^4=16 rows correspond to different initial bitstring states that are prepared and directly measured. In a perfect measurement scheme, this matrix would be the identity matrix. However, due to noise, other states may become populated. This file can be used to mitigate the effects of noisy measurements by applying the Iterative Bayesian Unfolding (IBU, see [[3]](#references)) method, as done in the notebook `compute populations.ipynb`
      - **`result_disorder_realization_<i>_backend=ibm_hanoi.json`**: Results from circuits in the `RC/` folder. Each file contains a list of 100 dictionaries. Each dictionary represents the quantum run of a cRC circuit version, with keys corresponding to bitstring configurations (in base 10) and values representing their measured probabilities.
      - **`result_ZNE3_disorder_realization_<i>_backend=ibm_hanoi.json`**: Results from circuits in the `ZNE3/` folder, with the same structure as above.

- **`populations/`**:
  - Stores population data derived from raw quantum data using the `compute populations.ipynb` notebook. Includes:
    - `population_exact_dynamics_<initial_state>.npy`: The theoretical  population expected from the exact hamiltonian dynamics. A NumPy array \( 100 x 1500 x 16 \). Axes correspond to:
      1. Disorder realizations (100).
      2. Simulation times (1500, sampled at intervals of 0.01 from 0 to 15).
      3. Bitstring configurations (16).
    - `population_noiseless_<initial_state>.npy`: the population expected from the recompiled circuit but on a perfect quantum computer. It is very close to the 	exact dynamics and therefore not shown in the paper. A NumPy array \( 3 x 16 x 100 \). Axes correspond to:
      1. Simulation times (3: \( T = 2.5, 5, 10 \)).
      2. Bitstring configurations (16).
      3. Disorder realizations (100).
    - `population_RC_<initial_state>.npy`: Population data calculated from the files `result_disorder_realization_<i>_backend=ibm_hanoi.json` using IBU method, combined with the `measurement_backend=ibm_hanoi.json` calibration file to correct for the effects of noisy measurements. A NumPy array \( 3 x 16 x 100 x 100 \). Axes correspond to:
      1. Simulation sessions (3: \( T = 2.5, 5, 10 \)).
      2. Bitstring configurations (16).
      3. Disorder realizations (100).
      4. Randomized circuit versions (100).
    - `population_ZNE3_<initial_state>.npy`: Same structure as `population_RC_<initial_state>.npy`, but derived from ZNE circuits.

### `figure/`
Contains Figures 2 and 3 of the article, generated using `plot article.ipynb`.

## Python Files
Three Python scripts containing essential functions:

- **`recompilation_functions.py`**: For recompiling circuits.
- **`randomized_compiling_functions.py`**: Implements Crosstalk Randomized Compiling (cRC) as described in [[2]](#references).
- **`ZNE.py`**: Implements Zero Noise Extrapolation (ZNE) by tripling the CNOT gates in the circuits.

## Jupyter Notebooks
Three notebooks for circuit generation, data processing, and analysis:

- **`circuits creation.ipynb`**: Generates quantum circuits and saves them in `data/circuits/`.
- **`compute populations.ipynb`**: Processes raw quantum data to compute populations and saves the results in `data/populations/`.
- **`plot article.ipynb`**: Analyzes population data and generates Figures 2 and 3 for the associated article.


## Dependencies
The project has been run using the following packages (with specific versions):

- Python 3.9.7
- qiskit 0.22.4
- quimb 1.5.1
- NumPy 1.24.0
- matplotlib 3.8.0
- scipy 1.11.4
- json 2.0.9
  

## Usage Instructions

1. **Generate Circuits**:
   - Run the `circuits creation.ipynb` notebook to generate quantum circuits. (Warning if you don't change the filename you will erase circuits already stored)

2. **Compute Populations**:
   - Process raw quantum data using `compute populations.ipynb`.

3. **Analyze Results**:
   - Use `plot article.ipynb` to generate Figures 2 and 3 for publication.

## References

1. Dynamic thermalization on noisy quantum hardware, H Perrin, T Scoquart, AI Pavlov, NV Gnezdilov, arXiv preprint [arXiv:2407.04770](https://arxiv.org/abs/2407.04770) (accepted in Communication Physics).

2.  Crosstalk Randomized Compiling (cRC): H. Perrin, T. Scoquart, A. Shnirman, J. Schmalian, and K. Snizhko, Mitigating crosstalk errors by randomized compiling: Simulation of the BCS model on a superconducting quantum computer, [Phys. Rev. Res. 6, 013142](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.013142) (2024).

3. Iterative Bayesian Unfolding (IBU): B. Nachman, M. Urbanek, W. A. de Jong, and C. W. Bauer, Unfolding quantum computer readout noise, [npj Quantum Information 6](https://www.nature.com/articles/s41534-020-00309-7) (2020).
