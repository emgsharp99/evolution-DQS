from IPython.display import HTML
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import numpy.random as random
import os
import pandas as pd
from scipy import linalg as ln
from scipy import sparse as sparse

epsilon = 3.0
spacing = 160
V_type = 'QHO'
seed = 1
max_steps = 750

class Wave_Packet:
    def __init__(self, epsilon, spacing, V_type, seed, dt = 0.25, x_range = 40, no_steps = 300, sigma0 = 1.5, k0 = 5.5):
        # instantiating constants
        self.N = no_steps
        self.epsilon = epsilon
        self.spacing = spacing
        self.x, self.dx = np.linspace(-1*x_range/2, x_range/2, self.N, retstep = True)
        self.dt = dt
        self.seed = seed
        
        # generation of Gaussian wavepacket 
        norm = (2.0 * np.pi * sigma0**2)**(-0.25)
        self.psi = np.exp(-(self.x - 0)**2 / (4.0 * sigma0**2))
        self.psi = self.psi*np.exp(1.0j * k0 * self.x)
        self.psi *= norm
        
        ## psi for perturbed
        self.psi_perturbed = self.psi
        
        # generation of the Hamiltonian matrix and introduction of disorder
        V = np.zeros(self.N)
        self.V = V
        # free particle
        if V_type == 'Free':
            pass
        # QHO
        elif V_type == 'QHO':
            for i in range(50,250):
                V[i] = self.x[i]**2 / self.N
#         # morse potential
#         elif V_typle == 'Morse':
#             pass 
        
        # disorder with set seed
        perturbation = epsilon*random.uniform(-1, 1, size = self.N)
        # creating desired spacin in disorder
        disorder = np.zeros(self.N)
        disorder[::self.spacing + 1] = perturbation[::self.spacing + 1]
        
        self.V_perturbed = self.V + disorder
    
        # arrays for diagonal and off-diagonal values in Hamiltonian
        diag = np.zeros(self.N)
        offdiag = np.zeros(self.N)
        I = np.identity(self.N)
        
        ## normal wavepacket
        # creating the Hamiltonian for unperturbed wavepacket
        for i in range(self.N):
            diag[i] = 1 + V[i]
        for i in range(self.N):
            offdiag[i] = - 0.5
        H = np.zeros(shape=(self.N, self.N))
        for j in range(self.N):
            per = (j+1) % self.N  # periodic
            H[j, j] = diag[j]
            H[j, per] = offdiag[j]
            H[per, j] = offdiag[j]
        hamiltonian = H
        
        # Crank-Nicolson method to calculate time evolution of normal wavepacket
        backward = (I - self.dt / 2.0j * hamiltonian)
        forward = (I + self.dt / 2.0j * hamiltonian)
        self.evolution_matrix = ln.inv(backward).dot(forward)

        ## perturbed wavepacket
        # creating the Hamiltonian for perturbed wavepacket
        for i in range(self.N):
            diag[i] = 1 + V[i] + disorder[i]
        for i in range(self.N):
            offdiag[i] = - 0.5
        H = np.zeros(shape=(self.N, self.N))
        for j in range(self.N):
            per = (j+1) % self.N  # periodic
            H[j, j] = diag[j]
            H[j, per] = offdiag[j]
            H[per, j] = offdiag[j]
        hamiltonian_perturbed = H
        
        # Crank-Nicolson method to calculate time evolution of perturbed wavepacket
        
        backward = (I - self.dt / 2.0j * hamiltonian_perturbed)
        forward = (I + self.dt / 2.0j * hamiltonian_perturbed)
        self.evolution_matrix_perturbed = ln.inv(backward).dot(forward)        
      
    def evolve(self):
        ## normal wavepacket
        # evolves wave packet
        self.psi = self.evolution_matrix.dot(self.psi)
        # probability density of wave
        self.prob = np.real(self.psi * np.conjugate(self.psi))
        
        norm = sum(self.prob) # normalises to begin with, then should be 1
        self.prob /= norm
        self.psi /= norm**0.5
        
        ## perturbed wavepacket
        # evolves wave packet
        self.psi_perturbed = self.evolution_matrix_perturbed.dot(self.psi_perturbed)
        # probability density of wave
        self.prob_perturbed = np.real(self.psi_perturbed * np.conjugate(self.psi_perturbed))
        
        norm_perturbed = sum(self.prob_perturbed) # normalises to begin with, then should be 1
        self.prob_perturbed /= norm_perturbed
        self.psi_perturbed /= norm_perturbed**0.5
        
        ## inner product between psi and perturbed psi
        inner = abs(self.psi * np.conjugate(self.psi_perturbed))
        self.norm_inner = sum(inner)
        
        return self.prob, self.prob_perturbed, self.norm_inner

class Evolution_Generator:
    def __init__(self, wave_packet, max_steps = max_steps):
        self.wave_packet = wave_packet
        self.max_steps = max_steps
        self.wave_data = []
        self.inner_data = []
        
    def data_output(self):    
        for i in range(self.max_steps):
            self.wave_data.append(self.wave_packet.evolve())
            self.inner_data.append(self.wave_packet.evolve()[2])
       
        return self.wave_data, self.inner_data


wavep = Wave_Packet(epsilon, spacing, V_type, seed)
wave_data = Evolution_Generator(Wave_Packet(epsilon, spacing, V_type, seed)).data_output()[0]

fig, ax1 = plt.subplots()

axtext = fig.add_axes([0.0,0.95,0.1,0.05])
axtext.axis("off")

time = axtext.text(0.5, 0.5, str(0), ha="left", va="top")

psi, = ax1.plot([],[])
psi_perturbed, = ax1.plot([],[])

def animate(j):
    print(str(round(j*100/max_steps, 2)) +'% complete.')
    ax1.clear()
    ax1.set_title('Time evolution of perturbed and unperturbed wavepackets in {}'.format(V_type))
    ax1.set_xlabel('Position, a$_0$')
    ax1.set_ylabel('Probability density, $|Ψ(x)|^2$')
    ax1.set_ylim([0,0.3])
    ax1.plot(Wave_Packet(epsilon, spacing, V_type, seed).x, Wave_Packet(epsilon, spacing, V_type, seed).V/5, color = 'b', label = 'Potential (schematic)')
    ax1.fill_between(Wave_Packet(epsilon, spacing, V_type, seed).x, Wave_Packet(epsilon, spacing, V_type, seed).V/5, facecolor="none", hatch="/", edgecolor="b", linewidth=0.0)
    ax1.plot(Wave_Packet(epsilon, spacing, V_type, seed).x, wave_data[j][0], label = 'Wavepacket (unperturbed)')
    ax1.plot(Wave_Packet(epsilon, spacing, V_type, seed).x, wave_data[j][1], label = 'Wavepacket (perturbed, ε = {})'.format(epsilon))
    ax1.legend()
    psi.set_data(Wave_Packet(epsilon, spacing, V_type, seed).x, wave_data[j][0])
    psi_perturbed.set_data(Wave_Packet(epsilon, spacing, V_type, seed).x, wave_data[j][1])
    
    time.set_text(('Elapsed time: {:6.2f} fs').format(j * Wave_Packet(epsilon, spacing, V_type, seed).dt * 2.419e-2))
    
    return psi, psi_perturbed, time,

wave_animation = animation.FuncAnimation(fig, animate, frames = len(wave_data), interval = 5)
#plt.show()
#HTML(wave_animation.to_html5_video())
wave_animation.save(
    './animations/s_{}_e_{}.gif'.format(spacing, epsilon),
    writer = 'pillow',
    fps = 30
    )