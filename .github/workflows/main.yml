
import matplotlib.pyplot as plt
from qutip import fock, thermal_dm, wigner, Qobj
N = 10

r = 0.5
v = np.sqrt(1 - r**2)
alpha = 0.5
n = 0.01
kappa_tau = np.linspace(0, 0.35, 200)
alpha_values = [1]
grid_points = 100
x = np.linspace(-2, 2, grid_points)
p = np.linspace(-2, 2, grid_points)
X, P = np.meshgrid(x, p)

a = destroy(N)

def calculate_delta(alpha, kappa_tau):
    initial_state =  (( qutip.coherent(N, alpha)).unit()).proj()
    decay_factor = np.exp(-kappa_tau)
    decayed_state = decay_factor * initial_state + (1 - decay_factor) * thermal_dm(10, 0)
    W = wigner(decayed_state, x, p)

    delta = decayed_state.tr()
    return delta

delta_values = {alpha: [] for alpha in alpha_values}

for alpha in alpha_values:
    for k_tau in kappa_tau:
        delta_values[alpha].append(calculate_delta(alpha, k_tau))

plt.figure(figsize=(8, 6))
for alpha in alpha_values:
    plt.plot(kappa_tau, delta_values[alpha], label=f"alpha={alpha}")
plt.xlabel(r"$\kappa \tau$")
plt.ylabel(r"$\mu(t)$")
plt.legend()
plt.xlim(0, 0.35)
plt.grid()
plt.show()
