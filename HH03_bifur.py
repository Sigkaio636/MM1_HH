import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.integrate import solve_ivp
from scipy import stats
from scipy import optimize
from matplotlib.backend_bases import MouseButton

"""
Display the evolution of the limiting behavior 
according to values of the bifurcation parameter (I),
creating a bifurcation diagram for the reduced HH-model
"""

# Shifted Nernst equilibrium potentials at mV
EK = -12
ENa = 120
EL = 10.6
# Maximal conductances at mS/cm^2
gK = 36
gNa = 120
gL = 0.3
# Opening and closing rates
an = lambda V: 0.01 * (10 - V) / (np.exp(1 - V / 10) - 1)
bn = lambda V: 0.125 * np.exp(-V / 80)
am = lambda V: 0.1 * (25 - V) / (np.exp(2.5 - V / 10) - 1)
bm = lambda V: 4 * np.exp(-V / 18)
ah = lambda V: 0.07 * np.exp(-V / 20)
bh = lambda V: 1 / (np.exp(3 - V / 10) + 1)
# Voltage-sensitive steady-state parameters
ninf = lambda V: an(V) / (an(V) + bn(V))
minf = lambda V: am(V) / (am(V) + bm(V))
hinf = lambda V: ah(V) / (ah(V) + bh(V))
tn = lambda V: 1 / (an(V) + bn(V))
tm = lambda V: 1 / (am(V) + bm(V))
th = lambda V: 1 / (ah(V) + bh(V))
# Membrane capacitance and currents
C = 1  # µF/cm^2
I0 = 0  # µA/cm^2


# Visualization of phase portrait
nd = lambda n, v: an(v) * (1 - n) - bn(v) * n
Vd = (
    lambda n, v, i: (
        i
        - gK * n**4 * (v - EK)
        - gNa * minf(v) ** 3 * (0.89 - 1.1 * n) * (v - ENa)
        - gL * (v - EL)
    )
    / C
)


def reducHH(t, state):
    V, n, I = state
    return [
        (
            I
            - gK * n**4 * (V - EK)
            - gNa * minf(V) ** 3 * (0.89 - 1.1 * n) * (V - ENa)
            - gL * (V - EL)
        )
        / C,
        an(V) * (1 - n) - bn(V) * n,
        0,
    ]


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, projection="3d")

for i in range(0, 100, 1):
    I = i / 5
    s0 = [20, 0.2, I]
    sol = solve_ivp(reducHH, [0, 1000], s0, dense_output=True)
    t = np.linspace(950, 1000, 500)
    s = sol.sol(t)
    ax.plot(s[0], s[1], [i] * len(t), c="k", alpha=0.3)

ax.set_xlabel("V")
ax.set_ylabel("n")
ax.set_zlabel("I")

plt.show()
