# 1D Heisenberg Model Simulation using iTEBD

This repository contains a Python implementation of an **infinite Time-Evolving Block Decimation (iTEBD)**-inspired workflow for simulating a one-dimensional spin-$\frac{1}{2}$ Heisenberg system under an external magnetic field.

The project demonstrates tensor-network-based numerical computation in clean, modular Python using standard scientific libraries.

## Features

- Object-oriented implementation of an iTEBD-inspired simulation workflow
- Construction of a two-site Heisenberg Hamiltonian
- Imaginary-time evolution using Trotter decomposition
- Tensor contractions with NumPy
- Singular Value Decomposition (SVD) for truncation and state updates
- Magnetization measurement as a function of external magnetic field
- Scientific visualization with Matplotlib
- Structured logging for monitoring simulation progress
- SLURM/Bash batch script for execution on Linux HPC clusters

## Technologies Used

- Python
- NumPy
- SciPy
- Matplotlib
- Bash
- SLURM
- Logging
- Object-Oriented Programming (OOP)

## Project Structure
```text
.
├── .gitignore
├── README.md
├── itebd_simulation.py
└── run_simulation.sh
Installation
Clone the repository:

bash
git clone https://github.com/Amirsalimi96/itebd-heisenberg-simulation.git
cd itebd-heisenberg-simulation
Install the required packages:

bash
pip install numpy scipy matplotlib
Usage
Run the simulation locally with:

bash
python itebd_simulation.py
Running on an HPC Cluster
You can also run the simulation using the provided SLURM batch script:

bash
sbatch run_simulation.sh
Before submission, make sure to adapt the environment setup and resource requests to your cluster configuration.

Why This Project?
This project demonstrates:

Scientific programming in Python
Numerical linear algebra and SVD-based truncation
Tensor-style array operations
Modular and object-oriented code design
Familiarity with Linux, Bash, and HPC/SLURM workflows
Author
Amir Salimi

PhD Student in Condensed Matter Physics

GitHub: @Amirsalimi96
