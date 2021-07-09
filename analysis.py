import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.cbook as cbook
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import math
import os
import pandas as pd

#%matplotlib qt

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

no_seeds = 1
directory = r'D:\Quantum Tunneling Project\Code\data_for_analysis'
# where data files are saved
path = os.path.join(directory, str(no_seeds))
fig_save_path = '\\'.join(path.split('\\')[:3])

'''Analysis, recovers inner product data from data files'''
def analysis(file_dir, foldername, plot = True):
    value_list = []
    folder_path = os.path.join(file_dir, foldername)
    for folder in os.listdir(folder_path):
        if not folder.endswith('.csv'):
            value_list.append(folder)
    
    ## list of seeds
    seed_list = os.listdir(os.path.join(folder_path, value_list[0]))        
    
    ## getting parameters and datasets
    system_param = os.path.join(folder_path, 'sys_params.csv')     
    # parameters
    parameters = pd.read_csv(system_param)
    epsilontemp = parameters['epsilonvals']
    # reducing to only non-NaN values
    epsilonvals = []
    spacingvals = []
    for epsilon in epsilontemp:
        if not math.isnan(epsilon):
            epsilonvals.append(epsilon)
    spacingtemp = parameters['spacingvals']
    for spacing in spacingtemp:
        if not math.isnan(spacing):
            spacingvals.append(spacing)

    dt = parameters['dt'][0]
    no_steps = int(parameters['no_steps'][0])
    total_runtime = parameters['total_runtime'][0]
    iter_over = parameters['iter_over'][0]
    
    #iterating over all data sets
    all_inner_ave = []
    anderson_measure = []
    for folder in os.listdir(folder_path):
        inner_ave = np.zeros(no_steps)
        if not folder.endswith('.csv'):
            for file in os.listdir(folder_path + '\\{}'.format(folder)):
                with open(folder_path + '\\{}\\{}'.format(folder, file), "r") as data:
                    data_str = data.read()

                # get data into usable form and averaging
                data_list = data_str.split(' ')
                data_list.remove('')

                data_arr = np.zeros(no_steps)
                for i in range(no_steps):
                    data_arr[i] = float(data_list[i])

                inner_ave += data_arr

            inner_ave /= len(seed_list)
            all_inner_ave.append(inner_ave)
            
            anderson_measure.append(np.mean(inner_ave))
    
    if plot == True:    
        ## plotting
        fig, ax = plt.subplots(figsize = (8,6))
#         fig.subplots_adjust(bottom=0.2)
#         fig.subplots_adjust(right=0.2)

        x = (np.linspace(0, 1, int(no_steps)))*float(total_runtime)

        if len(epsilonvals) == 1:
            ax.set_title(
                'Inner product of unperturbed and perturbed state\n at each timestep with fixed ε = {}'.format(epsilonvals[0]))
            colors = [cm.jet(x) for x in np.linspace(0,1,len(spacingvals))]
            for i in range(len(spacingvals)):
                ax.plot(
                    x, 
                    all_inner_ave[i],
                    label='{}'.format(spacingvals[i]),
                    color=colors[i]
                    )
            ax.legend(
                bbox_to_anchor=(1.15,0.5),
                loc='center right',
                title='Spacing, s'
                )

        elif len(spacingvals) == 1:
            ax.set_title(
                'Inner product of unperturbed and perturbed state\n at each timestep with fixed spacing, s = {}'.format(spacingvals[0]))   
            colors = [cm.jet(x) for x in np.linspace(0,1,len(epsilonvals[::3]))] 
            for i in range(len(epsilonvals[::3])):
                ax.plot(
                    x, 
                    all_inner_ave[i],                 
                    label='{}'.format(epsilonvals[3*i]),
                    color=colors[i]
                    )
            ax.legend(
                bbox_to_anchor=(1.15,0.5),
                loc='center right',
                title='Disorder, ε'
                )
        else:
            return 'ERROR: INPUT DIMENSIONS ARE INCORRECT'

        ax.set_xlabel('Time (fs)')
        ax.set_ylabel('⟨ Ψ$_0$ | Ψ$_p$ ⟩')
        ax.set_ylim([-0.05,1.05])
        
        ax.grid()
        if iter_over == 'epsilon':
            fig.savefig(os.path.join(fig_save_path, 'Analysis/{}/e_{}.pdf'.format(str(no_seeds), epsilonvals[0])), format = 'pdf')
            print('Figure saved')
        if iter_over == 'spacing':
            fig.savefig(os.path.join(fig_save_path, 'Analysis/{}/s_{}.pdf'.format(str(no_seeds), spacingvals[0])), format = 'pdf')
            print('Figure saved')
        return anderson_measure, iter_over, spacingvals, epsilonvals
    
    else:
        return anderson_measure, iter_over, spacingvals, epsilonvals
    
    
'''Iterates analysis in order to form graphs for all data. Also takes average for heat map'''
def iterator(path):
    tag_list = []
    a_measure_list = []
    spacingvals = []
    epsilonvals = []
    
    # ensuring files are in ascending numerical order
    iter_tag = os.listdir(path)[0].split('_')[3]
    for foldername in os.listdir(path):
        tag_list.append(int(foldername.split('_')[4]))
    
    tag_list.sort()
    foldername_list = []
    for tag in tag_list:
        foldername_list.append('QTP_SYSTEM_RUN_{}_{}'.format(iter_tag, tag))
        
    for foldername in foldername_list:
        anderson_measure, iter_over, d, e = analysis(path, foldername)
        for i in anderson_measure:
            a_measure_list.append(i)

        if iter_over == 'spacing':
            spacingvals.append(int(d[0]))
            epsilonvals = e

        elif iter_over == 'epsilon':
            epsilonvals.append(float(e[0]))
            spacingvals = d

    A = np.reshape(a_measure_list, (len(spacingvals), len(epsilonvals)))
    return  A, spacingvals, epsilonvals

# plotting data
A, d, e = iterator(path)
fig, ax = plt.subplots(figsize = (10,10))
ax.set_xticks(np.arange(len(e)))
plt.xticks(rotation=90)
ax.set_xticklabels(e)
ax.set_yticks(np.arange(len(d)))
ax.set_yticklabels(d)
ax.set_ylabel('Spacing, $s$')
ax.set_xlabel('Disorder, $ε$')
ax.set_title('Average inner product, ⟨ $Ψ_0$ | $Ψ_p$ ⟩, \n between unperturbed and perturbed state')

# colour map
ave_inner_heat = ax.imshow(A, interpolation='bilinear', cmap = 'viridis')

# colour bar
cbar = plt.colorbar(ave_inner_heat)
cbar.set_label('Average ⟨ $Ψ_0$ | $Ψ_p$ ⟩', rotation=270, labelpad = 25)      

ax.grid(True)
fig.savefig(os.path.join(fig_save_path, 'Analysis/{}/inner_cmap.pdf'.format(str(no_seeds))), format = 'pdf')

plt.show()