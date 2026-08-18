"""
Infinite Time-Evolving Block Decimation (iTEBD) Algorithm
Demonstrating: OOP, NumPy Tensor Operations, SVD Truncation, and Scientific Visualization.
"""

import logging
from typing import Tuple, List
import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt

# Configure structured logging
logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.INFO,
    datefmt="%H:%M:%S",
)

class QuantumSpinChainSimulation:
    """
    Simulates 1D/Quasi-1D Quantum Spin systems using the iTEBD tensor network method.
    """
    def __init__(self, bond_dim: int = 10, delta: float = 0.01, max_steps: int = 500, tol: float = 1e-8):
        self.chi = bond_dim      # Virtual bond dimension (D)
        self.delta = delta       # Imaginary time step
        self.max_steps = max_steps
        self.tol = tol
        self.d = 2               # Local physical dimension (Spin-1/2)

        # Basic Spin-1/2 Operators (Pauli Matrices)
        self.sx = 0.5 * np.array([[0, 1], [1, 0]], dtype=complex)
        self.sy = 0.5 * np.array([[0, -1j], [1j, 0]], dtype=complex)
        self.sz = 0.5 * np.array([[1, 0], [0, -1]], dtype=complex)
        self.eye = np.eye(2, dtype=complex)

    def _build_two_site_hamiltonian(self, j_coupling: float, h_field: float) -> np.ndarray:
        """Constructs a two-site Heisenberg Hamiltonian with external field hz."""
        # Exchange interaction: J * (Sx.Sx + Sy.Sy + Sz.Sz)
        h_int = j_coupling * (
            np.kron(self.sx, self.sx) +
            np.kron(self.sy, self.sy) +
            np.kron(self.sz, self.sz)
        )
        # On-site Zeeman field term: -hz * Sz (split equally on 2 sites)
        h_field_term = -0.5 * h_field * (np.kron(self.sz, self.eye) + np.kron(self.eye, self.sz))
        
        h_total = h_int + h_field_term
        return h_total.real

    def _get_trotter_gate(self, h_two_site: np.ndarray) -> np.ndarray:
        """Calculates imaginary-time Trotter evolution operator: U = exp(-delta * H)."""
        dim = h_two_site.shape[0]
        gate_matrix = expm(-self.delta * h_two_site)
        return gate_matrix.reshape(self.d, self.d, self.d, self.d)

    def init_mps(self) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        """Initializes random Gamma tensors and normalized Lambda vectors for A/B sites."""
        np.random.seed(42)
        # Gamma tensors: shape (bond_left, physical_dim, bond_right)
        ga = np.random.rand(self.chi, self.d, self.chi)
        gb = np.random.rand(self.chi, self.d, self.chi)
        ga /= np.max(np.abs(ga))
        gb /= np.max(np.abs(gb))

        # Lambda singular value vectors
        la = np.ones(self.chi) / np.sqrt(self.chi)
        lb = np.ones(self.chi) / np.sqrt(self.chi)
        return [ga, gb], [la, lb]

    def _contract_theta(self, gamma: List[np.ndarray], lambda_: List[np.ndarray], site_idx: int) -> np.ndarray:
        """Constructs the two-site wave function tensor Theta."""
        a = site_idx % 2
        b = (site_idx + 1) % 2
        
        # Theta contraction using NumPy tensordot
        theta = np.tensordot(np.diag(lambda_[b]), gamma[a], axes=(1, 0))
        theta = np.tensordot(theta, np.diag(lambda_[a]), axes=(2, 0))
        theta = np.tensordot(theta, gamma[b], axes=(2, 0))
        theta = np.tensordot(theta, np.diag(lambda_[b]), axes=(3, 0))
        return theta  # Shape: (chi_L, d_A, d_B, chi_R)

    def itebd_step(self, gamma: List[np.ndarray], lambda_: List[np.ndarray], u_gate: np.ndarray, site_idx: int):
        """Applies local Trotter gate and performs SVD truncation."""
        a = site_idx % 2
        b = (site_idx + 1) % 2

        theta = self._contract_theta(gamma, lambda_, site_idx)

        # Apply 2-site gate U: shape (d, d, d, d)
        theta_u = np.tensordot(u_gate, theta, axes=([2, 3], [1, 2]))
        theta_u = np.transpose(theta_u, (2, 0, 1, 3))  # (chi_L, d_A, d_B, chi_R)

        chi_l, d_a, d_b, chi_r = theta_u.shape
        mat = theta_u.reshape(chi_l * d_a, d_b * chi_r)

        # SVD Truncation
        u, s, vd = np.linalg.svd(mat, full_matrices=False)

        # Keep leading singular values up to bond dimension chi
        s_trunc = s[:self.chi] / np.linalg.norm(s[:self.chi])
        u_trunc = u[:, :self.chi]
        vd_trunc = vd[:self.chi, :]

        # Update MPS tensors
        inv_lb = np.diag(1.0 / (lambda_[b] + 1e-12))
        gamma[a] = np.tensordot(inv_lb, u_trunc.reshape(chi_l, d_a, self.chi), axes=(1, 0))
        gamma[b] = np.tensordot(vd_trunc.reshape(self.chi, d_b, chi_r), inv_lb, axes=(2, 0))
        lambda_[a] = s_trunc

    def measure_magnetization(self, gamma: List[np.ndarray], lambda_: List[np.ndarray]) -> float:
        """Measures the expectation value of Sz magnetization."""
        theta = self._contract_theta(gamma, lambda_, 0)
        # Apply Sz operator on the first physical index
        op_sz = np.kron(self.sz, self.eye).reshape(self.d, self.d, self.d, self.d)
        theta_op = np.tensordot(op_sz, theta, axes=([2, 3], [1, 2]))
        theta_op = np.transpose(theta_op, (2, 0, 1, 3))

        exp_val = np.sum(theta.conj() * theta_op)
        return float(np.real(exp_val))

    def run_field_sweep(self, h_values: np.ndarray, j_coupling: float = -1.0) -> List[float]:
        """Runs the complete ground state simulation over an array of external magnetic fields."""
        magnetizations = []
        gamma, lambda_ = self.init_mps()

        logging.info("Starting magnetic field sweep...")
        for hz in h_values:
            h_two = self._build_two_site_hamiltonian(j_coupling, hz)
            u_gate = self._get_trotter_gate(h_two)

            # Convergence loop
            for step in range(self.max_steps):
                self.itebd_step(gamma, lambda_, u_gate, site_idx=0)
                self.itebd_step(gamma, lambda_, u_gate, site_idx=1)

            mz = self.measure_magnetization(gamma, lambda_)
            magnetizations.append(mz)
            logging.info(f"Field hz = {hz:4.2f} | Magnetization <Sz> = {mz:.6f}")

        logging.info("Sweep finished.")
        return magnetizations


if __name__ == "__main__":
    # 1. Instantiate the simulation class
    sim = QuantumSpinChainSimulation(bond_dim=8, delta=0.02, max_steps=150)

    # 2. Define magnetic field range
    hz_list = np.arange(0.0, 2.1, 0.2)

    # 3. Execute field sweep
    mz_results = sim.run_field_sweep(h_values=hz_list, j_coupling=-1.0)

    # 4. Plot magnetization curve
    plt.figure(figsize=(7, 4.5))
    plt.plot(hz_list, mz_results, "s-", color="#1f77b4", linewidth=2, label=r"$\langle S^z \rangle$")
    plt.title("Magnetization Curve vs External Field ($h_z$)")
    plt.xlabel(r"Magnetic Field ($h_z$)")
    plt.ylabel(r"Magnetization $\langle S^z \rangle$")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()

