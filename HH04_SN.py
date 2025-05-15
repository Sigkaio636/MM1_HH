import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.integrate import solve_ivp
from scipy import stats
from scipy import optimize
from matplotlib.backend_bases import MouseButton

"""
Visualization of Phase portrait for reduced HH model.

The parameters and rates have been changed to produced 
another bifurcation, Saddle Node.

Note that the Poincaré-Bendixson Theorem can also be applied.
"""


# Shifted Nernst equilibrium potentials at mV, supercritical Andronov-Hopf bifurcation
EK = -70
ENa = 50
EL = -81
# Maximal conductances at mS/cm^2
gK = 60
gNa = 30
gL = 0.4
# Opening and closing rates
an = lambda V: 0.0088 * (V + 40) / (1 - np.exp(-(V + 40) / 7))
bn = lambda V: 0.037 * np.exp(-(V + 40) / 40)
am = lambda V: 0.08 * (V + 56) / (1 - np.exp(-(V + 56) / 6.8))
bm = lambda V: 0.8 * np.exp(-(V + 56) / 18)
ah = lambda V: 0.006 * np.exp(-(V + 41) / 14.7)
bh = lambda V: 1.3 / (1 + np.exp(-(V + 41) / 7.6))
# Voltage-sensitive steady-state parameters
ninf = lambda V: an(V) / (an(V) + bn(V))
minf = lambda V: am(V) / (am(V) + bm(V))
hinf = lambda V: ah(V) / (ah(V) + bh(V))
tn = lambda V: 1 / (an(V) + bn(V))
tm = lambda V: 1 / (am(V) + bm(V))
th = lambda V: 1 / (ah(V) + bh(V))
# Membrane capacitance and currents
C = 1.9  # µF/cm^2
I0 = 2  # µA/cm^2


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


# Initialize the plots
fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(1, 2, 1)
Vmin, Vmax = (-100, 70)
#   Create the vector field for the phase portrait
N, V = np.mgrid[0:1:100j, Vmin:Vmax:100j]
Vdot = Vd(N, V, I0)
Ndot = nd(N, V)
speed = np.sqrt(Vdot**2 + Ndot**2)
strm = ax1.streamplot(
    V, N, Vdot, Ndot, density=2, color=speed, linewidth=1, cmap="viridis"
)

# Obtein the n-nullcline, where simply n = ninf(V)
Vgride = np.linspace(Vmin, Vmax, 200)
ax1.plot(Vgride, ninf(Vgride), "r-.", label="n-null")
# Obtain the V-nullcline, solving V'=0 via newton method for each value of V
nopt_Vnull = optimize.newton(Vd, [0.8] * 200, args=(Vgride, [I0] * 200), maxiter=200)
nopt_Vnull = [x if 0 <= x <= 1 else np.nan for x in nopt_Vnull]
ax1.plot(Vgride, nopt_Vnull, "r--", label="V-null")

ax1.legend()
ax1.set_xlabel("V")
ax1.set_ylabel("n")
ax1.set_xlim(Vmin, Vmax)
ax1.set_ylim(0, 1)

# Initialize the numerical integration for the trajectory
s0 = [-25, 0.2, I0]
sol = solve_ivp(reducHH, [0, 100], s0, dense_output=True)
t = np.linspace(0, 100, 10000)
s = sol.sol(t)
ax2 = fig.add_subplot(1, 2, 2)
(line_pote,) = ax2.plot(t, s[0])
ax2.set_xlabel("t")
ax2.set_ylabel("V")
ax1.plot(s0[0], s0[1], "bo")
ax1.plot(s[0], s[1], "b")

# Add Slider axis: [left, bottom, width, height]
ax_slider = plt.axes([0.25, 0.02, 0.65, 0.03])  # type: ignore
I_slider = Slider(ax_slider, "Applied current", 0.0, 30.0, valinit=I0)


# Update function for interactive window
def update(val):
    ax1.cla()  # Clean the axis of the streamplot to repaint with new I value

    I = I_slider.val
    Vdot = Vd(N, V, I)
    speed = np.sqrt(Vdot**2 + Ndot**2)
    strm = ax1.streamplot(
        V, N, Vdot, Ndot, density=2, color=speed, linewidth=1, cmap="viridis"
    )

    ax1.plot(Vgride, ninf(Vgride), "r-.", label="n-null")

    nopt_Vnull = optimize.newton(Vd, [0.8] * 200, args=(Vgride, [I] * 200), maxiter=200)
    nopt_Vnull = [x if 0 <= x <= 1 else np.nan for x in nopt_Vnull]
    ax1.plot(Vgride, nopt_Vnull, "r--", label="V-null")

    ax1.legend()
    ax1.set_xlabel("V")
    ax1.set_ylabel("n")

    s0 = [-25, 0.2, I]
    sol = solve_ivp(reducHH, [0, 100], s0, dense_output=True)
    s = sol.sol(t)
    line_pote.set_ydata(s[0])  # Update the trajectory of Action potential in ax2
    ax1.plot(s0[0], s0[1], "bo")
    ax1.plot(s[0], s[1], "b")

    fig.canvas.draw_idle()


# Connect slider to update function
I_slider.on_changed(update)

plt.show()
