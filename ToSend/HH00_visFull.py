import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.integrate import solve_ivp

"""
To set up the complete system of the Hodgkin Huxley model, 
and given an initial position to obtain the trajectory for 
different values of intensity (the bifurcation parameter).
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


def system(t, state):
    V, n, m, h, I = state
    return [
        (I - gK * n**4 * (V - EK) - gNa * m**3 * h * (V - ENa) - gL * (V - EL)) / C,
        an(V) * (1 - n) - bn(V) * n,
        am(V) * (1 - m) - bm(V) * m,
        ah(V) * (1 - h) - bh(V) * h,
        0,
    ]


# Initialize plots
fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(1, 3, 1)
ax1.set_ylim(0, 1)
ax1.set_title("Evolution of Fraction of open gates")
ax2 = fig.add_subplot(1, 3, 2)
ax3 = fig.add_subplot(1, 3, 3, projection="3d")
ax2.set_title("Evolution of Action potential")
ax3.set_xlim(0, 1)
ax3.set_ylim(0, 1)
ax3.set_zlim(0, 1)
ax3.set_title("Visualization of 4D trajectory")

# Initialize numeric integration
V0 = 20
s0 = np.array([V0, ninf(V0), minf(V0), hinf(V0), I0])  # initial state
sol = solve_ivp(system, [0, 100], s0, dense_output=True)
t = np.linspace(0, 100, 10000)  # time interval
s = sol.sol(t)  # variables separated

(linen,) = ax1.plot(t, s[1], label="n")
(linem,) = ax1.plot(t, s[2], label="m")
(lineh,) = ax1.plot(t, s[3], label="h")
ax1.legend()
(lineV,) = ax2.plot(t, s[0], label="V")
ax2.legend()
sc3 = ax3.scatter(s[1], s[2], s[3], c=s[0], cmap="hot")
ax3.set_xlabel("n")
ax3.set_ylabel("m")
ax3.set_zlabel("h")

# Add Slider axis: [left, bottom, width, height]
ax_slider = plt.axes([0.25, 0.02, 0.65, 0.03])  # type: ignore
I_slider = Slider(ax_slider, "Applied current", 0.0, 30.0, valinit=I0)


# Update function for interactive window
def update(val):
    # Recalculate the trajectory with new I value
    I = I_slider.val
    s0 = np.array([V0, ninf(V0), minf(V0), hinf(V0), I])
    sol = solve_ivp(system, [0, 100], s0, dense_output=True)
    s = sol.sol(t)

    # Update plots
    linen.set_ydata(s[1])
    linem.set_ydata(s[2])
    lineh.set_ydata(s[3])
    lineV.set_ydata(s[0])

    ax3.clear()
    ax3.scatter(s[1], s[2], s[3], c=s[0], cmap="hot")
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    ax3.set_zlim(0, 1)
    ax3.set_xlabel("n")
    ax3.set_ylabel("m")
    ax3.set_zlabel("h")

    fig.canvas.draw_idle()


# Connect slider to update function
I_slider.on_changed(update)

plt.show()

"""
Check the number of zeros of V' plotting the graph
"""

DtV = (
    lambda V: (
        I0
        - gK * ninf(V) * 4 * (V - EK)
        - gNa * minf(V) ** 3 * hinf(V) * (V - ENa)
        - gL * (V - EL)
    )
    / C
)


V_ = np.linspace(-1000, 1000)
plt.plot(V_, DtV(V_))
plt.xlabel("V")
plt.ylabel("V'")
# plt.show()
plt.plot(V_, np.sign(DtV(V_)))
plt.xlabel("V")
plt.ylabel("sign(V')")
# plt.show()
