import numpy as np
import numpy.random as random
import os
import pandas as pd
import pickle
from scipy import linalg as ln
from scipy import sparse as sparse

seed_list = [50]
iter_list = ['spacing', 'epsilon']
V_type = 'QHO'

for seed in seed_list:
    for iter in iter_list:
        seeds = np.arange(seed)
        iter_over = iter#'epsilon'
        path = os.path.join(os.getcwd(), *['data_for_analysis_{}'.format(V_type), iter_over, str(len(seeds))+'.4'])

        try:
            os.makedirs(path)
        except FileExistsError:
            pass

        epsilonvals = []
        spacingvals = []
        for i in range(31):
            epsilonvals.append(i/10)
            
        for i in range(20,21):
            spacingvals.append(i)

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
            def __init__(self, wave_packet, max_steps = 750):
                self.wave_packet = wave_packet
                self.max_steps = max_steps
                self.wave_data = []
                self.inner_data = []
                
            def data_output(self):    
                for i in range(self.max_steps):
                    self.wave_data.append(self.wave_packet.evolve())
                    self.inner_data.append(self.wave_packet.evolve()[2])
            
                return self.wave_data, self.inner_data
            

        def data_toexport(epsilonvals, spacingvals, V_type, seeds, path, iter_over, fresh): 
            ## ensuring that it iterates over multiple runs of one variable 
            if iter_over == 'epsilon':
                epsilonvals = [epsilonvals]
                runstr = 'e_{}'.format(str(epsilonvals[0]))
            
            elif iter_over  == 'spacing':
                spacingvals = [spacingvals]
                runstr = 's_{}'.format(str(spacingvals[0]))
                
            ## if fresh data required, clears containing folder
            if fresh == True:
                print('CURRENT RUN: {}'.format(runstr))
                for folder in os.listdir(path):
                    folderdir = os.path.join(path,folder)
                    for subfolder in os.listdir(folderdir):#path + '/' + folder):
                        if subfolder.startswith('epsilon') or subfolder.startswith('spacing'):
                            filedir =  os.path.join(folderdir,subfolder)
                            for file in os.listdir(filedir):
                                os.remove(os.path.join(filedir,file))
                            os.rmdir(filedir)
                        elif subfolder.endswith('.csv'):
                            os.remove(folderdir)
                            
                    os.rmdir(folderdir)
            else:
                print('CURRENT RUN: {}'.format(runstr))
            
            ## system parameters  
            seeds = seeds
            newdir1 = os.path.join(path,'QTP_SYSTEM_RUN_{}'.format(str(runstr)))
            os.mkdir(newdir1)
            param_init = Wave_Packet(0, 0, V_type, 0)
            dt = param_init.dt
            no_steps = len(Evolution_Generator(param_init).data_output()[1])
            total_runtime = dt * no_steps * 2.419e-2
            sys_params = os.path.join(newdir1,'sys_params.csv')

            if iter_over == 'epsilon':
                epsilon = epsilonvals[0]
                for spacing in spacingvals:
                    print('--- SPACING: {}'.format(spacing))
                    newdir2 = os.path.join(newdir1,'spacing_{}'.format(str(spacing)))
                    os.mkdir(newdir2)
                    for seed in seeds:
                        wavegen = Wave_Packet(epsilon, spacing, V_type, seed)
                        inner_data = Evolution_Generator(wavegen).data_output()[1]
                        print('------ SEED: {}'.format(seed))
                        np.random.seed(seed)

                        # write parameters to file
                        parameters = {'epsilonvals': epsilonvals, 'spacingvals': spacingvals, 'dt': dt,
                                    'no_steps': no_steps, 'total_runtime': total_runtime, 'iter_over': iter_over}
                        df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in parameters.items()]))
                        df.to_csv(sys_params, index = False)

                        #write data to pickle file
                        data = {'inner_data': inner_data}
                        outbox = os.path.join(newdir2, 'SEED_{}.p'.format(str(seed)))
                        pickle.dump(data,open(outbox,'wb'))

                return

            elif iter_over == 'spacing':
                spacing = spacingvals[0]
                for epsilon in epsilonvals:
                    print('--- EPSILON: {}'.format(epsilon))
                    newdir2 = os.path.join(newdir1,'epsilon_{}'.format(str(epsilon)))
                    os.mkdir(newdir2)
                    for seed in seeds:
                        # setting data to generate and choosing save location
                        wavegen = Wave_Packet(epsilon, spacing, V_type, seed)
                        inner_data = Evolution_Generator(wavegen).data_output()[1]
                        print('------ SEED: {}'.format(seed))
                        np.random.seed(seed)

                        # write parameters to file
                        parameters = {'epsilonvals': epsilonvals, 'spacingvals': spacingvals, 'dt': dt,
                                    'no_steps': no_steps, 'total_runtime': total_runtime, 'iter_over': iter_over}
                        df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in parameters.items()]))
                        df.to_csv(sys_params, index = False)

                        #write data to pickle file
                        data = {'inner_data': inner_data}
                        outbox = os.path.join(newdir2, 'SEED_{}.p'.format(str(seed)))
                        pickle.dump(data,open(outbox,'wb'))              

                return
            
        def iterator(value, seeds):
            i = 0
            if value == 'epsilon':
                for epsilon in epsilonvals:
                    if i == 0:
                        data_toexport(epsilon, spacingvals, V_type, seeds, path, iter_over = value, fresh = True)
                        i = 1
                    elif i == 1:
                        data_toexport(epsilon, spacingvals, V_type, seeds, path, iter_over = value, fresh = False)
                        
            
            elif value == 'spacing':
                for spacing in spacingvals:
                    if i == 0:
                        data_toexport(epsilonvals, spacing, V_type, seeds, path, iter_over = value, fresh = True)
                        i = 1
                    elif i == 1:
                        data_toexport(epsilonvals, spacing, V_type, seeds, path, iter_over = value, fresh = False)
                print('COMPLETE')
            return

        iterator(value=iter_over, seeds=seeds)