import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import os
import pandas as pd
import pickle

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})


variable = 'epsilon'
no_seeds = '100'

path = os.path.join(os.getcwd(), *['data_for_analysis_QHO', variable, no_seeds])

if __name__ == '__main__':
    tag_list = []
    for folder in os.listdir(path):
        for subfolder in os.listdir(os.path.join(path, folder)):
            tag = str(subfolder.split('_')[1])
            if tag not in set(tag_list) and not subfolder.endswith('.csv'):
                tag_list.append(tag)
    
    tag_list1 = []
    tag_list2 = []
    for tag in tag_list:
        if len(tag) < 3:
            tag_list1.append(tag)
        elif len(tag) == 3:
            tag_list2.append(tag)
        else:
            pass
    
    tag_list1.sort()
    tag_list2.sort()

    tag_list = tag_list1 + tag_list

    tag_list = tag_list[0:5]

    folder_list = []
    
    for tag in tag_list:
        folder = 'spacing_' + tag
        folder_list.append(folder)

    epsilon_list = []
    spacing_ave_list = []
    for folder in os.listdir(path)[::3]:
        print('Working...')
        epsilon_list.append(folder)
        spacing_ave = []
        for tag in folder_list:
            inner_ave = 0
            for seed in os.listdir(os.path.join(path, *[folder, tag])):
                inner_data = pickle.load(open(os.path.join(path, *[folder, tag, seed]), 'rb'))['inner_data']
                data_arr = np.array(inner_data)
                tval = len(data_arr)
                inner_ave += data_arr
            spacing_ave.append(np.sum(inner_ave/100)/tval)
        spacing_ave_list.append(spacing_ave)
    
    label_list = []
    for epsilon in epsilon_list:
        label_list.append(epsilon.split('_')[4])

    fig, ax = plt.subplots(figsize=(8,8))
    colors = [cm.jet(x) for x in np.linspace(0,1,len(epsilon_list))] 
    for i in range(len(spacing_ave_list)):
        ax.plot(
            tag_list, 
            spacing_ave_list[i], 
            color = colors[i],
            label = label_list[i]
            )
        ax.set_ylabel('Average inner product, ⟨ $Ψ_0$ | $Ψ_p$ ⟩')
        ax.set_xlabel('Spacing, $d$')
        ax.set_ylabel('Average inner product, ⟨ $Ψ_0$ | $Ψ_p$ ⟩')
        ax.set_title('Inner product ⟨ $Ψ_0$ | $Ψ_p$ ⟩ averaged\nover time for variable spacing values')

        ax.legend(
            title='Disorder, $ε$',
            bbox_to_anchor=(1.15,0.5),
            loc='center right'
            )
    plt.show()
