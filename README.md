# 1D Heisenberg Model Simulation using iTEBD

This repository contains a Python implementation of the **infinite Time-Evolving Block Decimation (iTEBD)** algorithm for simulating 1D quantum spin systems (Heisenberg antiferromagnet).

## Features
- Implementation of Matrix Product States (MPS) and tensor contractions.
- Ground state energy calculation using imaginary time evolution.
- Numerical simulation and data visualization using NumPy, SciPy, and Matplotlib.
- Bash scripts for batch execution on Linux/HPC clusters.

## Requirements
- Python 3.8+
- NumPy
- SciPy
- Matplotlib

## Usage
To run the ground state energy simulation:
```bash
python itebd_simulation.py
