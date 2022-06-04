from IPython.display import HTML
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import numpy.random as random
import os
import pandas as pd
from scipy import linalg as ln
from scipy import sparse as sparse
import sympy as sympy
from tqdm import tqdm


class Wave_Packet:
    def __init__(self, epsilon, spacing, dt = 0.25, x_0=0, x_range=40, resolution=100, sigma0=1.5, k0=3.0):
        # instantiating constants
        self.N = resolution
        self.epsilon = epsilon
        self.spacing = spacing
        self.x, self.dx = np.linspace(-1*x_range/2, x_range/2, self.N, retstep = True)
        self.dt = dt
        self.seed = 1
        self.V = np.zeros(resolution)
        random.seed(self.seed)

        # disorder with set seed
        perturbation = epsilon*random.uniform(-1, 1, size=resolution)
        # creating desired spacin in disorder
        self.disorder = np.zeros(resolution)
        self.disorder[::spacing + 1] = perturbation[::spacing + 1]
        
        # generation of Gaussian wavepacket 
        norm = (2.0 * np.pi * sigma0**2)**(-0.25)
        self.psi = np.exp(-(self.x - x_0)**2 / (4.0 * sigma0**2))
        self.psi = self.psi*np.exp(1.0j * k0 * self.x)
        self.psi *= norm
        
        ## psi for perturbed
        self.psi_perturbed = self.psi    
        
    def add_barrier(self, loc, curve=0, preview=False):
        x = sympy.symbols("x")
        equation = curve
        equation = sympy.lambdify(x, equation, "numpy")

        loc_start=int(loc.split(':')[0][1:])
        loc_end=int(loc.split(':')[1][:-1])
        loc=range(loc_start, loc_end)
        
        x_cent = round(len(self.x)/2)
        barr_cent = round(len(loc)/2)

        self.V[loc]=equation(self.x[x_cent-barr_cent:x_cent+barr_cent])
        self.V_perturbed = self.V + self.disorder
        self.V_perturbed_graphed = self.V + (self.disorder/15)

        if preview==True:
            fig, ax = plt.subplots(figsize=(10,10))
            ax.plot(self.x, self.V_perturbed_graphed/5)
            ax.set_ylim([0,0.3])
            ax.fill_between(self.x, self.V_perturbed_graphed/5)
            ax.set_title("Preview of potential")
            ax.grid()
            plt.show()
            accept = input("Accept? (Y/N) \n")
            if accept=="N".lower():
                print('Exiting')
                exit()
              

    def evolve(self):
        # arrays for diagonal and off-diagonal values in Hamiltonian
        diag = np.zeros(self.N)
        offdiag = np.zeros(self.N)
        I = np.identity(self.N)
        
        ## normal wavepacket
        # creating the Hamiltonian for unperturbed wavepacket
        for i in range(self.N):
            diag[i] = 1 + self.V[i]
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
            diag[i] = 1 + self.V_perturbed[i]
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
    def __init__(self, wave_packet, max_steps):
        self.wave_packet = wave_packet
        self.max_steps = max_steps
        self.wave_data = []
        self.inner_data = []
        
    def data_output(self):
        print('Simulating...')
        for i in tqdm(range(self.max_steps)):
            # print(str(round(i*100/self.max_steps, 2)) +'% complete.')
            self.wave_data.append(self.wave_packet.evolve())
            self.inner_data.append(self.wave_packet.evolve()[2])
        return self.wave_data, self.inner_data
       

def viewer(wave_packet, max_steps=300, save=True, PATH="./test.gif"):
    fig, ax1 = plt.subplots(figsize=(10,10))
    axtext = fig.add_axes([0.0,0.95,0.1,0.05])
    axtext.axis("off")
    time = axtext.text(0.5, 0.5, str(0), ha="left", va="top")
    psi, = ax1.plot([],[])
    psi_perturbed, = ax1.plot([],[])

    wave_data = Evolution_Generator(wave_packet, max_steps).data_output()[0]

    def animate(j):
        ax1.clear()
        ax1.set_title('Time evolution of perturbed and\nunperturbed wavepackets')
        ax1.set_xlabel('Position, a$_0$')
        ax1.set_ylabel('Probability density, $|Ψ(x)|^2$')
        ax1.set_ylim([0,0.3])

        ax1.plot(
            wave_packet.x, 
            wave_packet.V_perturbed_graphed/5, 
            color='b', 
            label='Potential (schematic)'
            )
        ax1.fill_between(
            wave_packet.x, 
            wave_packet.V_perturbed_graphed/5, 
            facecolor="blue",
            alpha=0.3, 
            edgecolor="b", 
            linewidth=0.0
            )
        ax1.plot(
            wave_packet.x, 
            wave_data[j][0], 
            label='Wavepacket (unperturbed)', 
            linewidth=3,
            color='orange'
            )
        ax1.plot(
            wave_packet.x, 
            wave_data[j][1], 
            label='Wavepacket (perturbed)', 
            linewidth=3
            )
        ax1.legend()
        ax1.grid()
        psi.set_data(wave_packet.x, wave_data[j][0])
        psi_perturbed.set_data(wave_packet.x, wave_data[j][1])
        
        time.set_text(('Elapsed time: {:6.2f} fs').format(j * wave_packet.dt * 2.419e-2))
        
        return psi, psi_perturbed, time,

    wave_animation = animation.FuncAnimation(fig, animate, frames=tqdm(range(len(wave_data))), interval=0.5)
 
    if save==True:
        print('Saving...')
        wave_animation.save(
        PATH,
        writer = 'pillow',
        fps = 30
        )
        print('Saved to {}.'.format(PATH))

    if save==False:
        plt.show()
        HTML(wave_animation.to_html5_video())
    return
