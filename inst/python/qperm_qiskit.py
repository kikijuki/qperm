"""Real Qiskit backend for qperm's quantum permutation test.

This actually builds and runs the quantum circuits (uniform Hadamard state preparation, a
phase oracle marking the exceeding sign-flips, Grover amplitude amplification, and a
maximum-likelihood readout over a schedule of Grover powers) on a statevector simulator.
It is faithful but only tractable for small problems: the sign-flip register uses one
qubit per subject, so keep the subject count small (<= ~12). It is called from R through
reticulate; it can also be run directly for testing.

Nothing here beats the classical engine in cost: constructing the oracle enumerates the
sign-flip null classically. The point is to run genuine quantum circuits that reproduce
the corrected p-value, as a demonstration of the method.
"""
import numpy as np


def _marked(Y, threshold, alternative):
    """Boolean array over the 2^n basis states: which sign-flips exceed the threshold."""
    Y = np.asarray(Y, dtype=float)
    n, V = Y.shape
    if n > 16:
        raise ValueError("The Qiskit backend enumerates 2^n oracle phases; "
                         "n_subjects must be <= 16 (ideally <= 12).")
    N = 1 << n
    idx = np.arange(N)
    bits = (idx[:, None] >> np.arange(n)[None, :]) & 1        # N x n, bit i of state
    S = np.where(bits == 1, -1.0, 1.0)                        # subject i flipped if bit set
    T = (S @ Y) / n                                           # N x V per-voxel statistic
    if alternative == "two.sided":
        T = np.abs(T)
    elif alternative == "less":
        T = -T
    elif alternative != "greater":
        raise ValueError("alternative must be two.sided, greater, or less")
    return T.max(axis=1) >= threshold


def _mle_theta(powers, goods, shots, grid=2000):
    th = np.linspace(0.0, np.pi / 2, grid + 1)[1:]
    powers = np.asarray(powers); goods = np.asarray(goods)
    best, best_ll = th[0], -np.inf
    for t in th:
        s2 = np.clip(np.sin((2 * powers + 1) * t) ** 2, 1e-12, 1 - 1e-12)
        ll = float(np.sum(goods * np.log(s2) + (shots - goods) * np.log(1 - s2)))
        if ll > best_ll:
            best_ll, best = ll, t
    return best


def estimate_p(Y, threshold, alternative="two.sided",
               powers=(0, 1, 2, 3, 4, 6, 8, 12, 16, 24), shots=100, seed=1234):
    """Run the real Grover/MLQAE circuits and return the estimated corrected p-value.

    Returns a dict with the estimate, the exact marked fraction (for reference), and the
    per-power success counts.
    """
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import grover_operator
    from qiskit.primitives import StatevectorSampler

    Y = np.asarray(Y, dtype=float)
    n = Y.shape[0]
    powers = [int(m) for m in powers]
    shots = int(shots)

    marked = _marked(Y, float(threshold), alternative)
    exact_fraction = float(marked.mean())

    # Phase oracle: a diagonal unitary that flips the sign of marked basis states.
    phases = np.where(marked, -1.0, 1.0).astype(complex)
    oracle = QuantumCircuit(n)
    oracle.unitary(np.diag(phases), range(n), label="oracle")
    grover = grover_operator(oracle)   # default prep = H^n, reflection about |0>

    sampler = StatevectorSampler(seed=np.random.default_rng(seed))
    goods = []
    for m in powers:
        qc = QuantumCircuit(n)
        qc.h(range(n))
        for _ in range(m):
            qc = qc.compose(grover)
        qc.measure_all()
        result = sampler.run([qc], shots=shots).result()
        counts = result[0].data.meas.get_counts()
        good = sum(c for bit, c in counts.items() if marked[int(bit, 2)])
        goods.append(int(good))

    theta = _mle_theta(powers, goods, shots)
    return {
        "estimate": float(np.sin(theta) ** 2),
        "exact_fraction": exact_fraction,
        "powers": powers,
        "good": goods,
        "shots": shots,
    }


if __name__ == "__main__":
    # Quick self-test on the bundled synthetic dataset (10 subjects, 6 voxels).
    Y = np.array([
        [1.0512, 0.2987, -0.2741, -0.8906, -0.4547, -0.9916],
        [1.1101, 1.3402, -0.4922, -0.6205, 0.4898, 0.3569],
        [1.1554, -0.9305, -0.0293, 0.6953, -1.3442, -0.4576],
        [-0.8512, -1.2895, -1.8417, -0.2351, -1.2674, 0.2713],
        [1.2068, -0.1869, -2.5168, -0.5387, -0.0485, 0.1133],
        [-0.4801, -0.4778, -0.9785, -0.8088, 1.0609, -0.8075],
        [1.0175, 0.8844, -0.5836, -0.1117, 0.1105, 0.0638],
        [-0.1751, 0.0761, 1.3588, -1.5471, 0.8594, 0.1194],
        [0.4085, 2.0004, 0.7623, -1.1993, 0.0745, 0.5767],
        [0.8612, 0.6829, -0.0665, 0.6672, 1.4385, -0.6757]])
    # threshold chosen to give a smallish marked fraction
    thr = np.quantile(
        np.abs(((np.where(((np.arange(1 << 10)[:, None] >> np.arange(10)[None, :]) & 1) == 1,
                           -1.0, 1.0)) @ Y) / 10).max(axis=1), 0.95)
    out = estimate_p(Y, thr, "two.sided", shots=200)
    print("exact marked fraction:", round(out["exact_fraction"], 4))
    print("MLQAE estimate       :", round(out["estimate"], 4))
    print("per-power good counts :", out["good"])
