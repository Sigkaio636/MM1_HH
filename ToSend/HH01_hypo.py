import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.integrate import solve_ivp
from scipy import stats

"""
Simulate the complete Hodgkin Huxley model to study 
the assumptions that allow to obtain an approximate 
2D reduced system. In particular:

-> The variables n & h follow aproximately a linear relation
    h = 0.89 - 1.1 n 
-> The time constant of m gates is much smaller than for n and h
    Consider intantaneous transition m = m_inf
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
I0 = 15  # µA/cm^2


def fullHH(t, state):
    V, n, m, h, I = state
    return [
        (I - gK * n**4 * (V - EK) - gNa * m**3 * h * (V - ENa) - gL * (V - EL)) / C,
        an(V) * (1 - n) - bn(V) * n,
        am(V) * (1 - m) - bm(V) * m,
        ah(V) * (1 - h) - bh(V) * h,
        0,
    ]


# Visualization of whole system an reduction hypothesis check
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(2, 2, 1)
ax1.set_ylim(0, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax3.set_xlim(0, 1)
ax3.set_ylim(0, 1)
ax4 = fig.add_subplot(2, 2, 4)

ax1.set_title("Time evolution of open gates probability")
ax2.set_title("Time evolution of action potential")
ax3.set_title("n-h plane and linear tendencies")
ax4.set_title("Time constant of each gate")

s0 = np.array([0, ninf(0), minf(0), hinf(0), I0])
sol = solve_ivp(fullHH, [0, 100], s0, dense_output=True)
t = np.linspace(0, 100, 10000)
s = sol.sol(t)  # variables separated

(linen,) = ax1.plot(t, s[1], label="n")
(linem,) = ax1.plot(t, s[2], label="m")
(lineh,) = ax1.plot(t, s[3], label="h")
ax1.legend()
(lineV,) = ax2.plot(t, s[0], label="V")
ax2.legend()

(line_nh,) = ax3.plot(s[1], s[3])
ax3.set_xlabel("n")
ax3.set_ylabel("h")
n_gride = np.linspace(0, 1, 10)
(line_aprox,) = ax3.plot(n_gride, 0.89 - 1.1 * n_gride, "--g", label="aprox tendency")
nh_slope, nh_ntercept, r, p, std_err = stats.linregress(s[1], s[3])
(line_regre,) = ax3.plot(
    n_gride,
    nh_ntercept + nh_slope * n_gride,
    "--r",
    label=f"regression @ {r:.4f}",
)
ax3.legend()

(lineTaun,) = ax4.semilogy(s[0], tn(s[0]), label="n")
(lineTaum,) = ax4.semilogy(s[0], tm(s[0]), label="m")
(lineTauh,) = ax4.semilogy(s[0], th(s[0]), label="h")
ax4.set_xlabel("Observed V")
ax4.legend()

# Slider axis: [left, bottom, width, height]
ax_slider = plt.axes([0.25, 0.02, 0.65, 0.03])  # type: ignore
I_slider = Slider(ax_slider, "Applied current", 0.0, 30.0, valinit=I0)


def updateFull(val):
    I = I_slider.val
    s0 = np.array([0, ninf(0), minf(0), hinf(0), I])
    sol = solve_ivp(fullHH, [0, 100], s0, dense_output=True)
    s = sol.sol(t)  # variables separated

    linen.set_ydata(s[1])
    linem.set_ydata(s[2])
    lineh.set_ydata(s[3])
    lineV.set_ydata(s[0])
    line_nh.set_xdata(s[1])
    line_nh.set_ydata(s[3])

    nh_slope, nh_ntercept, r, p, std_err = stats.linregress(s[1], s[3])
    line_regre.set_ydata(
        nh_ntercept + nh_slope * n_gride,
    )
    line_regre.set_label(f"regression @ {r:.4f}")
    ax3.legend()

    lineTaun.set_xdata(s[0])
    lineTaun.set_ydata(tn(s[0]))
    lineTaum.set_xdata(s[0])
    lineTaum.set_ydata(tm(s[0]))
    lineTauh.set_xdata(s[0])
    lineTauh.set_ydata(th(s[0]))

    fig.canvas.draw_idle()


# Connect slider to update function
I_slider.on_changed(updateFull)

plt.show()
