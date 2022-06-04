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


variable = 'bigboi//spacing'
no_seeds = '50'
n = 4

figsave = os.getcwd()
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
        if len(tag) < 2:
            tag_list1.append(tag)
        elif len(tag) == 2:
            tag_list2.append(tag)
        else:
            pass
    
    tag_list1.sort()
    tag_list2.sort()

    tag_list = tag_list1 + tag_list

    folder_list = []
    
    for tag in tag_list:
        folder = 'epsilon_' + tag
        folder_list.append(folder)

    spacing_list = os.listdir(path)
    spacing_list = sorted(spacing_list, key=lambda x: int(x.split('_')[4]))
    epsilon_ave_list = []

    label_list = []

    for folder in spacing_list[::n]:
        label_list.append(folder.split('_')[4])
        print('Working...')
        epsilon_ave = []
        for tag in folder_list:
            inner_ave = 0
            for seed in os.listdir(os.path.join(path, *[folder, tag])):
                inner_data = pickle.load(open(os.path.join(path, *[folder, tag, seed]), 'rb'))['inner_data']
                data_arr = np.array(inner_data)
                tval = len(data_arr)
                inner_ave += data_arr
            epsilon_ave.append(np.sum(inner_ave/100)/tval)
        epsilon_ave_list.append(epsilon_ave)

    fig, ax = plt.subplots(figsize=(8,8))
    colors = [cm.jet(x) for x in np.linspace(0,1,len(spacing_list[::n]))] 
    for i in range(len(epsilon_ave_list)-1):
        ax.plot(
            tag_list, 
            epsilon_ave_list[i+1], 
            color = colors[i+1],
            label = label_list[i+1]
            )
    ax.set_ylabel('Average inner product, ⟨ $Ψ_0$ | $Ψ_p$ ⟩')
    ax.set_xlabel('Dpsilon, $ε$')
    ax.set_ylabel('Average inner product, ⟨ $Ψ_0$ | $Ψ_p$ ⟩')
    ax.set_title('Inner product ⟨ $Ψ_0$ | $Ψ_p$ ⟩ averaged\nover time for variable spacing values')

    labels=ax.xaxis.get_ticklabels()
    i=0
    for i in range(len(labels)):
        if i%3 > 0:
            labels[i].set_visible(False)
            i+=1
        else:
            i+=1

    ax.legend(
        title='Spacing, $d$',
        bbox_to_anchor=(1.15,0.5),
        loc='center right'
        )

    ax.grid()
    fig.savefig(os.path.join(figsave, *['Analysis', 'SOAnalysis.pdf']), format = 'pdf')
    plt.show()
