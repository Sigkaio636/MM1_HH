# MM1_HH
> Source: Dynamical Systems in Neuroscience, Eugene M. Izhikevich

### Motivation 
_Our goal is to motivate the reader to think of a neuron not only
in terms of ions and channels, as many biologists do, and not only in terms of an input/
output relationship, as many theoreticians do, but also as a nonlinear dynamical
system that looks at the input through the prism of its own intrinsic dynamics._

## Introduction 
_If somebody were to put a gun to the head of the author of this book and ask him to
name the single most important concept in brain science, he would say it is the concept
of a neuron._

There are only about 10^11 neurons in the brain and they are so important because they can transmit electrical signals over long distances.

![image](https://github.com/user-attachments/assets/992ee22b-863c-41b3-bd2a-6e55d203f23d)

Electrical activity in neurons is sustained and propagated via ionic currents through neuron membranes.
Most of these transmembrane currents involve one of four ionic species: sodium (Na+), potassium (K+), calcium (Ca2+), or chloride (Cl−).

The concentrations of these ions are different on the inside and the outside of a cell, which creates electrochemical gradients – the major driving forces of neural activity. The
extracellular medium has a high concentration of Na+ and Cl− (salty, like seawater) and a relatively high concentration of Ca2+. The intracellular medium has high concentrations of K+ and negatively charged molecules (denoted by A−).

![HHdata](https://github.com/user-attachments/assets/ceefe3d0-1b76-4054-941b-5f157ab7378a)

The cell membrane has large protein molecules forming channels through which ions (but not A−) can flow according to their electrochemical gradients.

• **Passive redistribution**. The impermeable anions A− attract more K+ into the cell (opposites attract) and repel more Cl− out of the cell, thereby creating concentration gradients.

• **Active transport**. Ions are pumped in and out of the cell via ionic pumps. For example, the Na+-K+ pump depicted in Fig.2.1 pumps out three Na+ ions for every two K+ ions pumped in, thereby maintaining concentration gradients.

The positive and negative charges accumulate on the opposite sides of the membrane surface, creating an electric potential gradient across the membrane – transmembrane potential or membrane voltage. At some point the concentration gradient and the electric potential gradient exert equal and opposite
forces that counterbalance each other, and the net cross-membrane current is zero. An equilibrium is achieved. The value of such an equilibrium potential depends on the ionic species, and it is given by the Nernst equation.

### What is a Spike?
A typical neuron receives inputs from more than 10, 000 other neurons through the contacts on its dendritic tree called synapses. The inputs produce electrical transmembrane currents that change the membrane potential of the neuron.

Synaptic currents produce changes, called postsynaptic potentials (PSPs). Small currents produce small PSPs; larger currents produce significant PSPs that can be amplified by the voltage-sensitive channels embedded in the neuronal membrane and lead to the generation of an action potential or spike – an abrupt and transient change of membrane voltage that propagates to other neurons via a long protrusion called an axon.

Such spikes are the main means of communication between neurons. In general, neurons do not fire on their own; they fire as a result of incoming spikes from other neurons.

Most introductory neuroscience books describe neurons as integrators with a threshold: neurons sum incoming PSPs and “compare” the integrated PSP with a certain voltage value, called the firing threshold. If it is below the threshold, the neuron remains quiescent; when it is above the threshold, the neuron fires an all-or-none spike.
As you might expect, how this comparison arises and is made has much greater depth. 

![image](https://github.com/user-attachments/assets/bf7aa7bf-0e76-413c-8d24-516eaa993a3c)

A neuron can fire under different experiments. Highlights:

• Transient current of different intensity, where a certain amperage is needed for firing.

• Persistent intensity, where the delay to the spike depends on the amperage. 

• Sequencial pulses, that resonates with the frequency and evokes a spike response.

If excitatory inputs depolarize the membrane potential (i.e., bring it closer to the “firing threshold”), and inhibitory inputs hyperpolarize the potential and move it
away from the threshold, then how can the neuron in Fig.1.6 fire in response to the inhibitory input?

As we can understand, the firing threshold, if it exists, must be somewhere in the shaded region.

![image](https://github.com/user-attachments/assets/ca9ae9a2-f629-4b01-8073-3202d48438f6)

![image](https://github.com/user-attachments/assets/756af101-f961-4247-be4f-860dfefa460f)

### Why Do We Care?

Why would two neurons respond completely differently to the same input? A biologist would say that the response of a neuron depends on many factors, such as the type of voltage- and Ca2+-gated channels expressed by the neuron, the morphology of its dendritic tree, the location of the input, and other factors. These factors are indeed important, but they do not determine the neuronal response per se. Rather they determine the rules that govern dynamics of the neuron. Different conductances and currents can result in the same rules, and hence in the same responses; conversely, similar currents can result in different rules and in different responses. The currents define what kind of dynamical system the neuron is, and so **neurons are dynamical systems**:

We divide all currents into two major classes: amplifying and resonant, with the persistent Na+ current INa,p and the persistent K+ current IK being the typical examples of the former and the latter, respectively.

A dynamical system consists of a set of variables that describe its state and a law that describes the evolution of the state variables with time. Typically, all variables describing neuronal dynamics can be classified into four classes, according to their function and the time scale.

1. Membrane potential.
2. Excitation variables, such as activation of a Na+ current. These variables are responsible for the upstroke of the spike.
3. Recovery variables, such as inactivation of a Na+ current and activation of a fast K+ current. These variables are responsible for the repolarization (downstroke) of the spike.
4. Adaptation variables, such as activation of slow voltage- or Ca2+-dependent currents. These variables build up during prolonged spiking and can affect excitability in the long run.

## Phase space



## The model, Ionic Currents

In the rest of the book V denotes the membrane potential and ENa, ECa, EK, and ECl denote the Nernst equilibrium potentials. The major ionic currents are:
$I_K=g_K(V−E_K)$,  $I_{Na}=g_{Na}(V−E_{Na})$,  $I_{Ca}=g_{Ca}(V−E_{Ca})$,  $I_{Cl}=g_{Cl}(V−E_{Cl})$, 

When the conductance is constant, the current is said to be Ohmic. In general, ionic currents in neurons are not Ohmic, since the conductances may depend on time, membrane potential, and pharmacological agents, e.g., neurotransmitters, neuromodulators, second-messengers, etc. It is the time-dependent variation in conductances that allows a neuron to generate an action potential, or spike.

It is traditional to represent electrical properties of membranes in terms of equivalent
circuits. 
According to Kirchhoff’s law, the total current, I, flowing across a patch of a cell membrane is the sum of the membrane capacitive current $C\dot{V}$ (the capacitance C ≈ 1.0 μF/cm2 in the squid axon) and all the ionic currents $I=C\dot{V}+I_{Na}+I_{Ca}+I_K+I_{Cl}$

Or equivalently: $C\dot{V}=I−g_{Na}(V−E_{Na})−g_{Ca}(V−E_{Ca})−g_K(V−E_K)−g_{Cl}(V−E_{Cl})$

![image](https://github.com/user-attachments/assets/7c2418cc-8da2-4fbc-8888-ed3d7a46bed4)

Let's analyse the conductances. Ionic channels are large transmembrane proteins having aqueous pores through which ions can flow down their electrochemical gradients. Despite the stochastic nature of transitions between open and closed states in individual channels, the net current generated by a large population or ensemble of identical channels can reasonably be described by the equation $I=\bar{g}p(V−E)$ where $p$ is the average proportion of channels in the open state, $\bar{g}$ is the maximal conductance of the population, and E is the reverse potential of the current, i.e., the potential at which the current reverses its direction. If the channels are selective for a single ionic species, then the reverse potential E equals the Nernst equilibrium potential.

When the gating particles are sensitive to the membrane potential, the channels are said to be voltage-gated. The gates are divided into two types: those that activate or
open the channels, and those that inactivate or close them. The probability of an activation gate being in the open state is denoted by the variable $m$ and the probability of an inactivation gate being in the open state is denoted by the variable $h$. Therefore, the proportion of open channels in a large population is $p=m^ah^b$ where $a$ is the number of activation gates and $b$ is the number of inactivation gates per channel.

Lastly, the time evolution of gates is given by the rates of transtion from open gates to close ones $\alpha$ (closing velocity) and from close to open $\beta$ (openning velocity). Therefore, for a generic gate $x$, the evolution law is given by $\dot{x}=\alpha_x(V)(1−x)−\beta_x(V)$. This rates are fitted empirically 

### The Hodgkin-Huxley model 
One of the most important models in computational neuroscience is the Hodgkin-Huxley model of the squid giant axon. Using pioneering experimental techniques of that time, Hodgkin and Huxley (1952) determined that the squid axon carries three major currents: voltage-gated persistent K+ current with four activation gates (resulting in the term $n^4$ in the equation below, where $n$ is the activation variable for K+); voltagegated transient Na+ current with three activation gates and one inactivation gate (the term $m^3h$ below), and Ohmic leak current, IL, which is carried mostly by Cl− ions. 

This model is a four-dimensional dynamical system because its state is uniquely determined by the membrane potential, $V$ , and socalled gating variables $n$, $m$, and $h$ for persistent K+ and transient Na+ currents. Note that The Hodgkin-Huxley model does not have variables of the fourth type (adaptation variables), but many neuronal models do, especially those exhibiting bursting dynamics.

The complete set of space-clamped Hodgkin-Huxley equations is

![image](https://github.com/user-attachments/assets/f641f6b8-c811-4dc3-bb79-daddbad9309a)

![image](https://github.com/user-attachments/assets/ac597262-88c1-41b7-8f30-47aaf42a8042)
