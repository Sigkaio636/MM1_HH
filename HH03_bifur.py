import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

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


def reducHH(t, state):
    V, n, I = state
    return [
        (
            I
            - gK * n**4 * (V - EK)
            - gNa * minf(V) ** 3 * (0.8882 - 1.041 * n) * (V - ENa)
            - gL * (V - EL)
        )
        / C,
        an(V) * (1 - n) - bn(V) * n,
        0,
    ]


fig = plt.figure(figsize=(16, 7))
axdyn = fig.add_subplot(1, 2, 1, projection="3d")
axFre = fig.add_subplot(1, 2, 2)
Fre_list = []
I_ = np.arange(-5, 20, 0.25)

for i in I_:
    s0 = [20, 0.2, i]
    sol = solve_ivp(reducHH, [0, 1000], s0, dense_output=True)
    # t = np.linspace(950, 1000, 1000)
    t = np.linspace(500, 1000, 5000)
    s = sol.sol(t)
    axdyn.plot(s[0], s[1], [i] * len(t), c="k", alpha=0.3)

    spectrum = np.abs(np.fft.fft(s[0][1:]))
    peaks, properties = find_peaks(spectrum, height=np.max(spectrum) / 2)

    k = 0 if len(peaks) == 0 else peaks[0]
    Fs = k / (t[-1] - t[0])
    Fre_list.append(Fs)

axdyn.set_xlabel("V", fontsize=14)
axdyn.set_ylabel("n", fontsize=14)
axdyn.set_zlabel("I", fontsize=14)

axFre.plot(I_, Fre_list, "--o")
axFre.set_xlabel("Intensity density µA/cm^2")
axFre.set_ylabel("Fundamental frequency Hz")

plt.show()
"""

I = 10 / 5
s0 = [20, 0.2, I]
sol = solve_ivp(reducHH, [0, 1000], s0, dense_output=True)
t = np.linspace(600, 1000, 4000)
s = sol.sol(t)

spectrum = np.abs(np.fft.fft(s[0][1:]))

plt.plot(s[0])
plt.show()
plt.plot(spectrum)
plt.show()

peaks, properties = find_peaks(spectrum, height=np.max(spectrum) / 2)
print(peaks, properties)

k = 0 if not peaks else peaks[0]
Fs = k / (t[-1] - t[0])
print(k, Fs)
"""
