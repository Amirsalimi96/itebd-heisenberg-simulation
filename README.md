# 1D Heisenberg Model Simulation using iTEBD

This repository contains a Python implementation of the **infinite Time-Evolving Block Decimation (iTEBD)** algorithm for simulating a 1D quantum spin-\(\frac{1}{2}\) Heisenberg system under an external magnetic field.

The project demonstrates how tensor-network-inspired numerical methods can be implemented in **clean, modular Python** using standard scientific libraries.

---

## Features

- Object-oriented implementation of an iTEBD-inspired simulation workflow
- Construction of the two-site Heisenberg Hamiltonian
- Imaginary-time evolution using **Trotter decomposition**
- Tensor contractions with **NumPy**
- Singular Value Decomposition (**SVD**) for truncation and state update
- Magnetization measurement as a function of external magnetic field
- Scientific plotting with **Matplotlib**
- Structured logging for monitoring simulation progress

---

## Technologies Used

- **Python**
- **NumPy**
- **SciPy**
- **Matplotlib**
- **Logging**
- **Object-Oriented Programming (OOP)**

---

## Project Structure
```bash
.
├── .gitignore
├── README.md
└── itebd_tensor_simulation.py

## Installation

Clone the repository:
```bash
git clone https://github.com/Amirsalimi96/itebd-heisenberg-simulation.git
cd itebd-heisenberg-simulation

pip install numpy scipy matplotlib

