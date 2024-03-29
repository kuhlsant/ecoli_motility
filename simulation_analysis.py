# %% initialize
import pickle
import os
import numpy as np
import random
from matplotlib import pyplot as plt
import matplotlib.cm as cm

# change figticks params back to before, motile/notmotile labels

timestep, box_size = 1/1800, 100

# %% Simulation Analysis Code
class simulation_data_lite:

    def __init__(self, data_dictionary):
        self.data_dictionary = data_dictionary
        self.mapsize = self.data_dictionary['size']
        self.timelen = self.data_dictionary['timerange']
        self.timestep, self.box_size = timestep, box_size
        self.k_on, self.t_on, self.t_std = self.data_dictionary['k_on'], self.data_dictionary['on_time_mean'], self.data_dictionary['on_time_std']

    def __str__(self):
        ktext, ttext = '$k_{on}$', '$t_{on}$'        
        return f'k_on = {self.k_on}, t_on = {self.t_on} +/- {self.t_std}'
    
    def full_popsize(self):
        return self.data_dictionary['fullpopsize']
    
    def full_expressfrac(self):
        return self.data_dictionary['expressfrac']
    
    def full_statfrac(self):
        return self.data_dictionary['statfrac']
    
    def explored_territory(self):
        return self.data_dictionary['explored_territory']
    
    def deadsize(self):
        return self.data_dictionary['deadsize']

class simulation_data:

    def __init__(self, data_dictionary):
        self.data_dictionary = data_dictionary
        self.mapsize = len(self.data_dictionary['population_sizes'])
        self.timelen = len(self.data_dictionary['population_sizes'][0])
        self.timestep, self.box_size = timestep, box_size
        self.k_on, self.t_on, self.t_std = self.data_dictionary['k_on'], self.data_dictionary['on_time_mean'], self.data_dictionary['on_time_std']


    def __str__(self):
        ktext, ttext = '$k_{on}$', '$t_{on}$'        
        return f'k_on = {self.k_on}, t_on = {self.t_on} +/- {self.t_std}'
        
    def plot_populationmap(self, y_grain=900, method='absolute'):
        x_len = self.mapsize 
        # 1) get timestep for readout
        
        big_t_step = int((self.timelen-1)//y_grain)
        values = []
        for i in range(y_grain):
            values.append([])
            for j in range(x_len):
                values[i].append(self.data_dictionary['population_sizes'][j][i*big_t_step+1])
                
        for i in range(y_grain):
            for j in range(x_len):
                if values[i][j] == 0: values[i][j] = np.nan
                
        if method == 'relative':
            for i, t_list in enumerate(values):
                values[i] = [np.nan if np.nansum(t_list) == 0 else val/np.nansum(t_list) for val in t_list]
        
        values = np.array(values)
        
        # Plotting
        plt.figure(figsize=(3, 9), dpi=240)  
        plt.imshow(values, cmap='viridis', aspect='auto', interpolation='none')
        cbar = plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.08)
        cbar.set_label('Population Size', rotation=0, labelpad=5)

        plt.xlabel('X [mm]')
        plt.ylabel('time [h]')
        plt.xticks([0, 25, 50], [0, 2.5, 5]) #np.arange(0, x_len, x_len/5), np.arange(0, x_len, x_len/5)/5) # 5 boxes = 1mm
        plt.yticks(np.arange(0, y_grain, y_grain/5), [0, 12, 24, 36, 48]) #np.arange(0, y_grain, y_grain/5)*big_t_step*self.timestep) # 900 timesteps = 1h
        
        plt.show()
    
    def plot_foodmap(self, y_grain=900):
        x_len = self.mapsize  # number of tiles in the x direction
        # y_len = 60   # number of timesteps in the y direction
        
        # create a data array
        # 1) get timestep for readout
        
        big_t_step = int((self.timelen-1)//y_grain)
        values_food = []
        for i in range(y_grain):
            values_food.append([])
            for j in range(x_len):
                values_food[i].append(self.data_dictionary['food_distribution'][j][i*big_t_step+1])
        
        values_food = np.array(values_food)
        
        plt.imshow(values_food, cmap='magma', aspect='auto', interpolation='none')
        plt.colorbar(label='available food units')
        
        plt.xlabel('X [mm]')
        plt.ylabel('time [h]')
        plt.xticks(np.arange(0, x_len, x_len/5), np.arange(0, x_len, x_len/5)/5) # 5 boxes = 1mm
        plt.yticks(np.arange(0, y_grain, y_grain/5), np.arange(0, y_grain, y_grain/5)*big_t_step*self.timestep) # 900 timesteps = 1h
        
        plt.show()
       
    def plot_expressionmap(self, y_grain=900):
        x_len = self.mapsize  # number of tiles in the x direction
        
        # create a data array
        # 1) get timestep for readout
        
        big_t_step = int((self.timelen-1)//y_grain)
        values_expr = []
        for i in range(y_grain):
            values_expr.append([])
            for j in range(x_len):
                values_expr[i].append(self.data_dictionary['expressing_fractions'][j][i*big_t_step+1])
        
        values_expr = np.array(values_expr)

        plt.figure(figsize=(3, 9), dpi=240)          
        plt.imshow(values_expr, cmap='cividis', aspect='auto', interpolation='none')
        cbar = plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.08)
        cbar.set_label('fraction of expressing bacteria', rotation=0, labelpad=5)
        
        plt.xlabel('X [mm]')
        plt.ylabel('time [h]')
        plt.xticks([0, 25, 50], [0, 2.5, 5]) #np.arange(0, x_len, x_len/5), np.arange(0, x_len, x_len/5)/5) # 5 boxes = 1mm
        plt.yticks(np.arange(0, y_grain, y_grain/5), [0, 12, 24, 36, 48]) #np.arange(0, y_grain, y_grain/5), np.arange(0, y_grain, y_grain/5)*big_t_step*self.timestep) # 900 timesteps = 1h
        
        plt.show()

    def full_popsize(self):
        fullpopsize = [sum([voxl[t] for voxl in self.data_dictionary['population_sizes']]) for t in range(self.timelen)]
        return fullpopsize

    def full_expressfrac(self):
        expressfrac = [np.nansum([voxl[t] for voxl in self.data_dictionary['expressing_fractions']])/self.mapsize for t in range(self.timelen)]
        return expressfrac

    def full_statfrac(self):
        statfrac = [np.nansum([voxl[t] for voxl in self.data_dictionary['stationary_fractions']])/self.mapsize for t in range(self.timelen)]
        return statfrac

    def explored_territory(self):
        explored_voxels, fraction_explored = [], []
        for t in range(self.timelen):
            explored_voxels.append([])
            if t == 0:
                for i in range(self.mapsize):
                    if self.data_dictionary['population_sizes'][i][t] != 0: explored_voxels[-1].append(1)
                    else: explored_voxels[-1].append(0)
            elif len(explored_voxels)>1:
                for i in range(self.mapsize):
                    if explored_voxels[-2][i] == 1: explored_voxels[-1].append(1)
                    elif explored_voxels[-2][i] == 0 and self.data_dictionary['population_sizes'][i][t] != 0: explored_voxels[-1].append(1)
                    else: explored_voxels[-1].append(0)
                explored_voxels.remove(explored_voxels[0])
            fraction_explored.append(sum(explored_voxels[-1])/self.mapsize)
        return fraction_explored
    
    def motile_intervals(self):
        motilehist = self.data_dictionary['motile_history']
        motileints = [] # timeintervals of motility

        for bact in motilehist:
            counter = 0
            for i, val in enumerate(bact): # val == True/False
                if val: 
                    counter += 1
                elif val == False and counter != 0:
                    motileints.append(counter)
                    counter = 0
        
        # Multiply the values in timesteps by self.timestep
        motileints = [interval * self.timestep for interval in motileints]
        
        return motileints

class concurring_species:

    def __init__(self, condition_name, lite=True, namelist = False):
        #self.simulist1, self.simulist2 = [], [] # list storing all simulation_data
        self.namelist = namelist
        # get number of employed combinations and repeats (with while they will count one too far)       
        datafiles = os.listdir()
        self.combi = 0
        self.fractiles_minmax = False
        rep = 0
        for fil in datafiles:
            if condition_name in fil:
                while 'combi'+str(self.combi) in fil:
                    self.combi += 1
                while 'r'+str(rep) in fil:
                    rep += 1
        print(str(self.combi)+' combinations in this dataset, '+str(rep)+' replicates were done.')

        self.timestep, self.box_size = timestep, box_size
        if lite: 
            self.lite = True
            print('species 1')
            self.spec1 = [condition(condition_name, 'spec1', 'combi'+str(c), lite=True) for c in range(self.combi)]
            print('species 2')
            self.spec2 = [condition(condition_name, 'spec2', 'combi'+str(c), lite=True) for c in range(self.combi)]

        elif lite == False:
            self.lite = False
            self.spec1 = [condition(condition_name, 'spec1', 'combi'+str(c), lite=False) for c in range(self.combi)]
            self.spec2 = [condition(condition_name, 'spec2', 'combi'+str(c), lite=False) for c in range(self.combi)]

        self.spec1_ts, self.spec1_ks = [], [] # on_time_mean, k_on
        for cond in self.spec1: self.spec1_ts += cond.t_list; self.spec1_ks += cond.k_list
        self.spec2_ts, self.spec2_ks = [], []
        for cond in self.spec2: self.spec2_ts += cond.t_list; self.spec2_ks += cond.k_list

        #self.some_combi = combination(condition_name, random.randint(0, self.combi))

    def get_populationdifference(self, combin, plot=True):
        '''gives fraction of species 1 over time'''
        if self.lite:
            # got lost in the structure of this input, the code gets too large now
            pop1 = self.spec1[combin].mean_popsizes()[0][0][0]
            pop2 = self.spec2[combin].mean_popsizes()[0][0][0]
            popdiff = np.array(pop1)/(np.array(pop1)+np.array(pop2))

            if plot:
                fig, axs = plt.subplots(1, 1)
                axs.plot(self.spec1[combin].time, pop1, label='species 1', linewidth=0.5)
                axs.plot(self.spec2[combin].time, pop2, label='species 2', linewidth=0.5)
                axs.plot(1, 1, label='combi '+str(combination), linewidth=3, alpha=0.7, color='olive')

                axs2 = axs.twinx()
                axs2.plot(self.spec1[combin].time, popdiff, label='combi '+str(combination), linewidth=3, alpha=0.7, color='olive')
                axs.set_xlabel('time [h]')
                axs2.set_ylabel('fraction of species 1')
                axs.legend()
                plt.show()

        else: 
            print('nothing to see here')
            pass

        return popdiff
    
    def plot_tiles_with_values(self, value_matrix, value_name):
        '''
        we assume that combinations were carried out in the order that 0 is in the top left and n in the bottom right (so we just fill the full square up)
        '''
        # Define the values for the variables

        edge_length = int(len(value_matrix))
        x_len, y_len = edge_length, edge_length

        value_matrix = np.array(value_matrix)

        if self.namelist:
            variable1_values = self.namelist[0]
            variable2_values = self.namelist[1]

        # k_size, t_size
        else:
            print('define self.namelist manually: [[ylabels], [xlabels]], color boundaries: self.fractiles_minmax = [min, max]')
            variable1_values = ['y' for i in range(y_len)]
            variable2_values = ['x' for i in range(x_len)]
        
        # Create a 3x3 grid of subplots
        fig, axs = plt.subplots(y_len, x_len, figsize=(2*y_len, 1.5*x_len))
        
        # Loop through each subplot
        for i in range(y_len):
            for j in range(x_len):
                # Plot the output value as a color-coded tile

                if self.fractiles_minmax: im = axs[i, j].imshow([[value_matrix[i, j]]], cmap='viridis', vmin=self.fractiles_minmax[0], vmax=self.fractiles_minmax[1])
                else: im = axs[i, j].imshow([[value_matrix[i, j]]], cmap='viridis', vmin=np.nanmin(value_matrix), vmax=np.nanmax(value_matrix))
        
                # Print the output value on the tile
                x_center = 0
                y_center = 0
                axs[i, j].text(x_center, y_center, f'{value_matrix[i, j]:.2f}', color='white', ha='center', va='center', fontsize=12)
        
                # Set the tick labels for the axes
                if i == edge_length-1:  # Bottom row
                    axs[i, j].set_xticks([0])
                    axs[i, j].set_xticklabels([variable2_values[j]])
                else:
                    axs[i, j].set_xticks([])
                if j == 0:  # Leftmost column
                    axs[i, j].set_yticks([0])
                    axs[i, j].set_yticklabels([variable1_values[i]])
                else:
                    axs[i, j].set_yticks([])
        
        # Add text descriptions for variable 1
        #axs[-1, 1].text(-2.1, -1.5, '$k_{on}$', color='black', ha='center', va='center', fontsize=12, rotation=90)
        # Add text descriptions for variable 2
        #axs[1, 0].text(1.2, 2.1, '$t_{on}$', color='black', ha='center', va='center', fontsize=12)
        
        # Add a color bar
        cbar = fig.colorbar(im, ax=axs.ravel().tolist(), orientation='vertical')
        cbar.set_label(value_name)
        
        plt.show()
  
    def plot_popfractiles(self, ignore=False, final_tile_nb = False, average_over=10):
        '''
        average over x last hours
        '''
        tile_nbs = [int(x*x) for x in range(20)]
        if len(self.spec1) not in tile_nbs:
            print('tileplot not possible with '+str(len(self.spec1))+' combinations, check stuff. Filling up to square with np.nans')
            if final_tile_nb == False:
                i = 0
                while tile_nbs[i]<len(self.spec1): i+=1
                final_tile_nb = tile_nbs[i]

        # get values
        value_list = []

        for c in range(self.combi):
            popdiff = self.get_populationdifference(c, plot=False) #, plot=False)
            mean_over = int(average_over/self.timestep) # timesteps, 6h
            popdiffmean, popdiffstd = np.nanmean(popdiff[-mean_over:]), np.nanstd(popdiff[-mean_over:])
            #if popdiffmean/popdiffstd > 0.2 and ignore == False:
            #    print('std is greater than 20%, is this sim equilibrated? In case activate ignore')
            #    self.get_populationdifference(c, plot=True)
            
            value_list.append(popdiffmean)
        
        if final_tile_nb:
            value_list += [np.nan]*(final_tile_nb-len(value_list))

        # distribute in value matrix, which is then handed over to the tile_plotter
        count = 0
        if final_tile_nb: edge_length = int(np.sqrt(final_tile_nb))
        else: edge_length = int(np.sqrt(len(self.spec1)))
        value_matrix = [[] for i in range(edge_length)]
        for i in range(edge_length):
            for j in range(edge_length):
                value_matrix[i].append(value_list[count])
                count += 1

        self.plot_tiles_with_values(value_matrix, 'fraction of species 1')

    def get_mean_and_std_bacfrac(self, combi, t_max):
        if t_max == 'last': t_max = -1
        elif t_max > 0: t_max = int(t_max/timestep)

        bac1sims = self.spec1[combi].simulist
        bac2sims = self.spec2[combi].simulist
        allbact_1 = [np.array(sim.data_dictionary['fullpopsize'])+np.array(sim.data_dictionary['deadsize']) for sim in bac1sims]
        allbact_2 = [np.array(sim.data_dictionary['fullpopsize'])+np.array(sim.data_dictionary['deadsize']) for sim in bac2sims]
   
        bacfracs = [alb1[t_max]/(allbact_2[i][t_max]+alb1[t_max]) for i, alb1 in enumerate(allbact_1)]
        mean_bacfracs = np.mean(bacfracs)
        std_bacfracs = np.std(bacfracs)

        print('combi: ', combi, ', mean = ', mean_bacfracs, ' std = ', std_bacfracs)
        return mean_bacfracs, std_bacfracs

    def plot_fullbacfractiles(self, ignore=False, final_tile_nb = False, timepoint='last'):
        '''
        fraction of all bacteria having lived up to a certain timepoint (last timepoint if not specified) (last timepoint in h)
        '''
        tile_nbs = [int(x*x) for x in range(20)]
        if len(self.spec1) not in tile_nbs:
            print('tileplot not possible with '+str(len(self.spec1))+' combinations, check stuff. Filling up to square with np.nans')
            if final_tile_nb == False:
                i = 0
                while tile_nbs[i]<len(self.spec1): i+=1
                final_tile_nb = tile_nbs[i]

        # get values
        value_list = []
        meanbs = [] # the means of the combinations in htis condition, for easier copy & paste
        for c in range(self.combi):
            bacs1 = self.spec1[c].all_bact_mean(t_max=timepoint)
            bacs2 = self.spec2[c].all_bact_mean(t_max=timepoint)
            bacfrac = bacs1/(bacs1+bacs2)
            
            value_list.append(bacfrac)

            mb, sb = self.get_mean_and_std_bacfrac(c, timepoint)
            meanbs.append(mb)
        
        print(meanbs)
        print(value_list)
        if final_tile_nb:
            value_list += [np.nan]*(final_tile_nb-len(value_list))

        # distribute in value matrix, which is then handed over to the tile_plotter
        count = 0
        if final_tile_nb: edge_length = int(np.sqrt(final_tile_nb))
        else: edge_length = int(np.sqrt(len(self.spec1)))
        value_matrix = [[] for i in range(edge_length)]
        for i in range(edge_length):
            for j in range(edge_length):
                value_matrix[i].append(value_list[count])
                count += 1

        self.plot_tiles_with_values(value_matrix, 'fraction of species 1')
    
    def plot_fullpoptiles(self, ignore=False, final_tile_nb = False):
        print('mean is taken over last 5400 timesteps, change if necessary')
        tile_nbs = [int(x*x) for x in range(20)]
        if len(self.spec1) not in tile_nbs:
            print('tileplot not possible with '+str(len(self.spec1))+' combinations, check stuff. Filling up to square with np.nans')
            if final_tile_nb == False:
                i = 0
                while tile_nbs[i]<len(self.spec1): i+=1
                final_tile_nb = tile_nbs[i]

        # get values
        value_list = []

        for c in range(self.combi):       
            fullpop = self.spec1[c].mean_popsizes()[0][0][0]+self.spec2[c].mean_popsizes()[0][0][0]
            mean_over = 5400 # timesteps, 6h
            fullpopmean, fullpopstd = np.nanmean(fullpop[-5400:]), np.nanstd(fullpop[-5400:])
            if fullpopmean/fullpopstd > 0.2 and ignore == False:
                print('std is greater than 20%, is this sim equilibrated? In case activate ignore')
            
            value_list.append(fullpopmean)
        
        if final_tile_nb:
            value_list += [np.nan]*(final_tile_nb-len(value_list))

        # distribute in value matrix, which is then handed over to the tile_plotter
        count = 0
        if final_tile_nb: edge_length = int(np.sqrt(final_tile_nb))
        else: edge_length = int(np.sqrt(len(self.spec1)))
        value_matrix = [[] for i in range(edge_length)]
        for i in range(edge_length):
            for j in range(edge_length):
                value_matrix[i].append(value_list[count])
                count += 1

        self.plot_tiles_with_values(value_matrix, 'final full population size')

    def plot_population(self, combi):
        sd1 = self.spec1[combi].simulist[0].data_dictionary['fullpopsize']
        sd2 = self.spec2[combi].simulist[0].data_dictionary['fullpopsize']
        time = np.arange(0, len(sd1), 1)*timestep

        plt.plot(time, sd1, label='species 1')
        plt.plot(time, sd2, label='species 2')
        plt.legend()
        plt.xlabel('time [h]')
        plt.ylabel('# bacteria')
        plt.show()
        
class combination:
      
    def __init__(self, *condition_names, combi_nb=0):

        nofiles = True
        datafiles = os.listdir()
        for fil in datafiles:
            stringcheck = all(subname in fil for subname in condition_names)
            if stringcheck and 'lite' not in fil and '_combi'+str(combi_nb)+'_' in fil:
                with open(fil, 'rb') as file: 
                    data_dict = pickle.load(file) 
                    print('opened: ', fil)
                    nofiles = False
                    if 'spec1' in fil: self.spec1 = data_dict
                    elif 'spec2' in fil: self.spec2 = data_dict
        
        if nofiles: raise ValueError('No files found for the given conditions. Remember: combination does not take lite files')

        #self.data_dictionary = data_dictionary
        self.mapsize = len(self.spec1['population_sizes'])
        self.timelen = len(self.spec1['population_sizes'][0])
        self.timestep, self.box_size = timestep, box_size
        self.spec1sim, self.spec2sim = simulation_data(self.spec1), simulation_data(self.spec2)

        # get population sizes. relpopsize: relative population size, showing prevalence of species. fullpopsize: overall count of bacteria
        self.relpopsize, self.fullpopsize = [], []
        for i, voxel in enumerate(self.spec1['population_sizes']):
            self.relpopsize.append([]); self.fullpopsize.append([])
            for j, timep in enumerate(voxel):
                if self.spec1['population_sizes'][i][j]+self.spec2['population_sizes'][i][j] == 0: 
                    relsize, fullsize = np.nan, np.nan
                else: 
                    relsize = self.spec1['population_sizes'][i][j]/(self.spec1['population_sizes'][i][j]+self.spec2['population_sizes'][i][j])
                    fullsize = self.spec1['population_sizes'][i][j]+self.spec2['population_sizes'][i][j]
                self.relpopsize[i].append(relsize); self.fullpopsize[i].append(fullsize)

    def __str__(self):
        k1, t1 = self.spec1['k_on'], self.spec1['on_time_mean']
        k2, t2 = self.spec2['k_on'], self.spec2['on_time_mean']
        return f'species 1: k_on = {k1}, t_on = {t1}, species 2: k_on = {k2}, t_on = {t2}'
        
    def plot_populationmap(self, y_grain=900, method='relative'):
        '''method either relative (prevalence of species) or absolute (overall nb of individuals)'''
        x_len = self.mapsize 
        
        big_t_step = int((self.timelen-1)//y_grain)
        values = []
        if method == 'relative':
            colormap = 'Spectral'
            for i in range(y_grain):
                values.append([])
                for j in range(x_len):
                    values[i].append(self.relpopsize[j][i*big_t_step+1])
        
        elif method == 'absolute':
            colormap = 'viridis'
            for i in range(y_grain):
                values.append([])
                for j in range(x_len):
                    values[i].append(self.fullpopsize[j][i*big_t_step+1])
            
        values = np.array(values)
        
        # Plotting
        plt.figure(figsize=(3, 9), dpi=240)  
        plt.imshow(values, cmap=colormap, aspect='auto', interpolation='none')
        cbar = plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.08)
        cbar.set_label('Population Size', rotation=0, labelpad=5)
        
        plt.xlabel('X [mm]')
        plt.ylabel('time [h]')
        plt.xticks([0, 25, 50], [0, 2.5, 5]) #np.arange(0, x_len, x_len/5), np.arange(0, x_len, x_len/5)/5) # 5 boxes = 1mm
        plt.yticks(np.arange(0, y_grain, y_grain/5), [0, 12, 24, 36, 48]) #np.arange(0, y_grain, y_grain/5), np.arange(0, y_grain, y_grain/5)*big_t_step*timestep) # 900 timesteps = 1h
        
        plt.show()
    
    def plot_newfood(self):
        newfoods = []
        wasfoods = [0]
        for t in range(self.timelen):
            isfood = sum([self.spec1['food_distribution'][vox][t] for vox in range(len(self.spec1['food_distribution']))])
            newfood = isfood-wasfoods[-1]
            newfoods.append(newfood)
            wasfoods.append(isfood)
        plt.plot(newfoods)

    def plot_foodmap(self, y_grain=900):
        x_len = len(self.spec2['food_distribution'])  # number of tiles in the x direction
        # y_len = 60   # number of timesteps in the y direction
        
        # create a data array
        # 1) get timestep for readout
        
        big_t_step = int((self.timelen-1)//y_grain)
        values_food = []
        for i in range(y_grain):
            values_food.append([])
            for j in range(x_len):
                values_food[i].append(self.spec1['food_distribution'][j][i*big_t_step+1])
        
        values_food = np.array(values_food)

        plt.figure(figsize=(3, 9), dpi=240)        

        plt.imshow(values_food, cmap='magma', aspect='auto', interpolation='none')
        cbar = plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.08)
        cbar.set_label('food distribution', rotation=0, labelpad=5)
        
        plt.xlabel('X [mm]')
        plt.ylabel('time [h]')
        plt.xticks([0, 25, 50], [0, 2.5, 5]) #np.arange(0, x_len, x_len/5), np.arange(0, x_len, x_len/5)/5) # 5 boxes = 1mm
        plt.yticks(np.arange(0, y_grain, y_grain/5), [0, 12, 24, 36, 48]) #np.arange(0, y_grain, y_grain/5), np.arange(0, y_grain, y_grain/5)*big_t_step*timestep) # 900 timesteps = 1h
        
        plt.show()

    def plot_foodfluxes(self, y_grain=900):
        x_len = 89 #len(self.spec1['food_flux'][0])  # number of tiles in the x direction
        # y_len = 60   # number of timesteps in the y direction
        
        # create a data array
        # 1) get timestep for readout
        big_t_step = int((self.timelen-1)//y_grain)
        values_food = []
        for i in range(y_grain):
            values_food.append([])
            for j in range(x_len):
                values_food[i].append(self.spec1['food_flux'][i*big_t_step+1][j])
        
        values_food = np.array(values_food)
        
        plt.imshow(values_food, cmap='magma', aspect='auto', interpolation='none')
        plt.colorbar(label='food fluxes')
        
        plt.xlabel('X [mm]')
        plt.ylabel('time [h]')
        plt.xticks(np.arange(0, x_len, x_len/5), np.arange(0, x_len, x_len/5)/5) # 5 boxes = 1mm
        plt.yticks(np.arange(0, y_grain, y_grain/5), np.arange(0, y_grain, y_grain/5)*big_t_step*timestep) # 900 timesteps = 1h
        
        plt.show()

    def plot_popsizes(self):

        fullpop1 = [sum([voxl[t] for voxl in self.spec1['population_sizes']]) for t in range(self.spec1['timerange'])]
        fullpop2 = [sum([voxl[t] for voxl in self.spec2['population_sizes']]) for t in range(self.spec2['timerange'])]

        fig, axs = plt.subplots(1, 1)
        time = [i*self.timestep for i in range(len(fullpop1))]
        axs.plot(time, fullpop1, label='species 1', color='tab:blue')
        axs.plot(time, fullpop2, label='species 2', color='tab:orange')

        axs.legend()
        axs.set_xlabel('time [h]')
        axs.set_ylabel('population size')
        plt.show()
    
    def plot_motile_histograms(self, intime=10, bins=100):
        motileints1 = self.spec1sim.motile_intervals()
        motileints2 = self.spec2sim.motile_intervals()
        bins = np.linspace(0, intime, bins)
        fig, axs = plt.subplots(1, 1)
        axs.hist(motileints1, bins, alpha=0.5, label='species 1', color='tab:blue')
        axs.hist(motileints2, bins, alpha=0.5, label='species 2', color='tab:orange')   
        axs.legend()
        axs.set_xlabel('motile interval [h]')
        axs.set_ylabel('count')
        plt.show()

        tres1 = ss_tres(self.spec1['on_time_mean'], self.spec1['k_on']*self.timestep)
        tres2 = ss_tres(self.spec2['on_time_mean'], self.spec2['k_on']*self.timestep)

        print('species 1: <motile interval> = '+str(np.mean(motileints1))+' h, std = '+str(np.std(motileints1))+' h, expected = '+str(tres1)+' h')
        print('species 2: <motile interval> = '+str(np.mean(motileints2))+' h, std = '+str(np.std(motileints2))+' h, expected = '+str(tres2)+' h')

class condition:

    def __init__(self, *condition_names, lite=True):
        self.simulist = [] # list storing all simulation_data
        self.lite = lite
        # get all data
        datafiles = os.listdir()
        for fil in datafiles:
            stringcheck = all(subname in fil for subname in condition_names)
            if stringcheck:
                if self.lite and 'lite' in fil:                    
                    with open(fil, 'rb') as file: data_dict = pickle.load(file)
                    self.simulist.append(simulation_data_lite(data_dict))
                    print('opened: ', fil)
                elif self.lite == False and 'lite' not in fil:
                    with open(fil, 'rb') as file: data_dict = pickle.load(file)
                    self.simulist.append(simulation_data(data_dict))
                    print('opened: ', fil)
        
        self.timestep, self.box_size = self.simulist[0].data_dictionary['timestep'], self.simulist[0].data_dictionary['box_size']

        #print('condition initialized with '+str(len(self.simulist))+' datasets')
        
        self.datasize = len(self.simulist)
        self.mapsize, self.time = self.simulist[0].mapsize, np.arange(0, self.simulist[0].timelen*self.timestep, self.timestep)[0:self.simulist[0].timelen]
        self.somesim = random.choice(self.simulist)

        # sort data by k and t
        self.k_list, self.t_list = [], []
        for sim in self.simulist:
            if sim.k_on not in self.k_list: self.k_list.append(sim.k_on)
            if sim.t_on not in self.t_list: self.t_list.append(sim.t_on)
        self.k_list.sort(), self.t_list.sort()

        self.simulist_ksort, self.simulist_tsort = [[sim for sim in self.simulist if sim.k_on == k] for k in self.k_list], [[sim for sim in self.simulist if sim.t_on == t] for t in self.t_list]

    def mean_popsizes(self):
        formean = []
        means_and_stds = []            
        for i, k in enumerate(self.k_list):
            formean.append([])
            means_and_stds.append([])
            for j, t in enumerate(self.t_list):
                formean[i].append([])
                means_and_stds[i].append([])
                for sim in self.simulist:
                    if sim in self.simulist_ksort[i] and sim in self.simulist_tsort[j]:
                        if self.lite: formean[i][j].append(sim.data_dictionary['fullpopsize'])
                        else: formean[i][j].append(sim.full_popsize())
                means_and_stds[i][j].append(np.nanmean(formean[i][j], axis=0)); means_and_stds[i][j].append(np.nanstd(formean[i][j], axis=0))        
        return means_and_stds

    def mean_deadsizes(self):
        formean = []
        means_and_stds = []            
        for i, k in enumerate(self.k_list):
            formean.append([])
            means_and_stds.append([])
            for j, t in enumerate(self.t_list):
                formean[i].append([])
                means_and_stds[i].append([])
                for sim in self.simulist:
                    if sim in self.simulist_ksort[i] and sim in self.simulist_tsort[j]:
                        if self.lite: formean[i][j].append(sim.data_dictionary['deadsize'])
                        else: formean[i][j].append(sim.full_popsize())
                means_and_stds[i][j].append(np.nanmean(formean[i][j], axis=0)); means_and_stds[i][j].append(np.nanstd(formean[i][j], axis=0))        
        return means_and_stds
    
    def just_meanpop(self, k=False, t=False):
        if len(self.k_list) == 1 and len(self.t_list) == 1:
            k = self.k_list[0]; t = self.t_list[0]
        elif k == False and t == False:
            print('ks: ', self.k_list)
            print('ts: ', self.t_list)
            raise ValueError('please define k and t, there are more than 1 possibilities')
        
        popdata = self.mean_popsizes()
        for i, ks in enumerate(self.k_list):
            for j, ts in enumerate(self.t_list):
                if ks == k and ts == t: 
                    meanpop = popdata[i][j][0]
        
        return meanpop
    
    def just_meandead(self, k=False, t=False):
        if len(self.k_list) == 1 and len(self.t_list) == 1:
            k = self.k_list[0]; t = self.t_list[0]
        elif k == False and t == False:
            print('ks: ', self.k_list)
            print('ts: ', self.t_list)
            raise ValueError('please define k and t, there are more than 1 possibilities')
        
        popdata = self.mean_deadsizes()
        for i, ks in enumerate(self.k_list):
            for j, ts in enumerate(self.t_list):
                if ks == k and ts == t: meandead = popdata[i][j][0]
        
        return meandead

    def all_bact_mean(self, k=False, t=False, t_max = 'last'):
        '''
        returns all bacteria having lived until a specified timepoint (averaged over repeats)
        '''
        if len(self.k_list) == 1 and len(self.t_list) == 1:
            k = self.k_list[0]; t = self.t_list[0]
        elif k == False and t == False:
            print('ks: ', self.k_list)
            print('ts: ', self.t_list)
            raise ValueError('please define k and t, there are more than 1 possibilities')
        
        meanpop = self.just_meanpop() # allstd is std for alive and dead
        meandead = self.just_meandead()

        if t_max == 'last': timepoint = -1
        else: timepoint = int(t_max/self.timestep)
        all_bact = meanpop[timepoint]+meandead[timepoint]

        bac1sims = self.simulist
        allbact_1 = [np.array(sim.data_dictionary['fullpopsize'])+np.array(sim.data_dictionary['deadsize']) for sim in bac1sims]
        #print('allbacs', [bbb[timepoint] for bbb in allbact_1])

        return all_bact

    def mean_expressfracs(self):
        formean = []
        means_and_stds = []
        for i, k in enumerate(self.k_list):
            formean.append([]); means_and_stds.append([])
            for j, t in enumerate(self.t_list):
                formean[i].append([]); means_and_stds[i].append([])
                for sim in self.simulist:
                    if sim in self.simulist_ksort[i] and sim in self.simulist_tsort[j]:
                        if self.lite: formean[i][j].append(sim.data_dictionary['expressfrac'])
                        else: formean[i][j].append(sim.full_expressfrac())
                means_and_stds[i][j].append(np.nanmean(formean[i][j], axis=0)); means_and_stds[i][j].append(np.nanstd(formean[i][j], axis=0))
        return means_and_stds

    def mean_statfracs(self):
        formean = []
        means_and_stds = []
        for i, k in enumerate(self.k_list):
            formean.append([]); means_and_stds.append([])
            for j, t in enumerate(self.t_list):
                formean[i].append([]); means_and_stds[i].append([])
                for sim in self.simulist:
                    if sim in self.simulist_ksort[i] and sim in self.simulist_tsort[j]:
                        if self.lite: formean[i][j].append(sim.data_dictionary['statfrac'])
                        else: formean[i][j].append(sim.full_statfrac())
                means_and_stds[i][j].append(np.nanmean(formean[i][j], axis=0)); means_and_stds[i][j].append(np.nanstd(formean[i][j], axis=0))
        return means_and_stds

    def mean_exploredter(self):
        formean = []
        means_and_stds = []
        for i, k in enumerate(self.k_list):
            formean.append([]); means_and_stds.append([])
            for j, t in enumerate(self.t_list):
                formean[i].append([]); means_and_stds[i].append([])
                for sim in self.simulist:
                    if sim in self.simulist_ksort[i] and sim in self.simulist_tsort[j]:
                        if self.lite: formean[i][j].append(sim.data_dictionary['explored_territory'])
                        else: formean[i][j].append(sim.explored_territory())
                means_and_stds[i][j].append(np.nanmean(formean[i][j], axis=0)); means_and_stds[i][j].append(np.nanstd(formean[i][j], axis=0))
        return means_and_stds
    
    def plot_popsizes(self, meansim=True, legend=True, xlim=False):
        if meansim:
            popdata = self.mean_popsizes()
            print(type(popdata))
            fig, axs = plt.subplots(1, 1)
            for i, k in enumerate(self.k_list):
                for j, t in enumerate(self.t_list):
                    axs.plot(self.time, popdata[i][j][0], linewidth=2, label=' $k_{on}$ = '+str(k)+', $t_{on}$ = '+str(t))
                    axs.fill_between(self.time, popdata[i][j][0]-popdata[i][j][1], popdata[i][j][0]+popdata[i][j][1], alpha=0.2)
            if legend: axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            if xlim: axs.set_xlim(xlim[0], xlim[1])
            axs.set_xlabel('time [h]')
            axs.set_ylabel('population size')
           
        else:
            popsizes = [sim.full_popsize() for sim in self.simulist]
            k_ons, t_ons = [sim.k_on for sim in self.simulist], [sim.t_on for sim in self.simulist]
            
            fig, axs = plt.subplots(1, 1)
            for i, pop in enumerate(popsizes):
                axs.plot(pop, label=' $k_{on}$ = '+str(k_ons[i])+', $t_{on}$ = '+str(t_ons[i]))
            axs.legend()
            plt.show()

    def plot_expressfracs(self, meansim=True, legend=True, xlim=False):
        if meansim:
            exprdata = self.mean_expressfracs()
            fig, axs = plt.subplots(1, 1)
            for i, k in enumerate(self.k_list):
                for j, t in enumerate(self.t_list):
                    axs.plot(self.time, exprdata[i][j][0], linewidth=2, label=' $k_{on}$ = '+str(k)+', $t_{on}$ = '+str(t))
                    axs.fill_between(self.time, exprdata[i][j][0]-exprdata[i][j][1], exprdata[i][j][0]+exprdata[i][j][1], alpha=0.2)
            if legend: axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            if xlim: axs.set_xlim(xlim[0], xlim[1])
            axs.set_xlabel('time [h]')
            axs.set_ylabel('expressing fraction')
        else:
            exprfracs = [sim.full_expressfrac() for sim in self.simulist]
            k_ons, t_ons = [sim.k_on for sim in self.simulist], [sim.t_on for sim in self.simulist]
            
            fig, axs = plt.subplots(1, 1)
            for i, expr in enumerate(exprfracs):
                axs.plot(expr, label=' $k_{on}$ = '+str(k_ons[i])+', $t_{on}$ = '+str(t_ons[i]))
            axs.legend()
            plt.show()

    def plot_statfracs(self, meansim=True, legend=True, xlim=False):
        if meansim:
            statdata = self.mean_statfracs()
            fig, axs = plt.subplots(1, 1)
            for i, k in enumerate(self.k_list):
                for j, t in enumerate(self.t_list):
                    axs.plot(self.time, statdata[i][j][0], linewidth=2, label=' $k_{on}$ = '+str(k)+', $t_{on}$ = '+str(t))
                    axs.fill_between(self.time, statdata[i][j][0]-statdata[i][j][1], statdata[i][j][0]+statdata[i][j][1], alpha=0.2)
            if legend: axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            if xlim: axs.set_xlim(xlim[0], xlim[1])
            axs.set_xlabel('time [h]')
            axs.set_ylabel('stationary fraction')

    def plot_exploredter(self, meansim=True, legend=True, xlim=False):
        if meansim:
            expldata = self.mean_exploredter()
            fig, axs = plt.subplots(1, 1)
            for i, k in enumerate(self.k_list):
                for j, t in enumerate(self.t_list):
                    axs.plot(self.time, expldata[i][j][0], linewidth=2, label=' $k_{on}$ = '+str(k)+', $t_{on}$ = '+str(t))
                    axs.fill_between(self.time, expldata[i][j][0]-expldata[i][j][1], expldata[i][j][0]+expldata[i][j][1], alpha=0.2)
            if legend: axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            if xlim: axs.set_xlim(xlim[0], xlim[1])
            axs.set_xlabel('time [h]')
            axs.set_ylabel('explored territory')
           
        else:
            expldata = [sim.explored_territory() for sim in self.simulist]
            k_ons, t_ons = [sim.k_on for sim in self.simulist], [sim.t_on for sim in self.simulist]
            
            fig, axs = plt.subplots(1, 1)
            for i, expl in enumerate(expldata):
                axs.plot(expl, label=' $k_{on}$ = '+str(k_ons[i])+', $t_{on}$ = '+str(t_ons[i]))
            axs.legend()
            plt.show()

    def plot_tiles_with_value(self, value_matrix, value_name):
        # Define the values for the variables
        variable1_values = self.k_list
        variable2_values = self.t_list
        ksize, tsize = len(self.k_list), len(self.t_list)
        
        # Create a 3x3 grid of subplots
        fig, axs = plt.subplots(ksize, tsize, figsize=(2*tsize, 1.5*ksize))
        
        # Loop through each subplot
        for i in range(ksize):
            for j in range(tsize):
                # Plot the output value as a color-coded tile
                im = axs[i, j].imshow([[value_matrix[i, j]]], cmap='viridis', vmin=np.min(value_matrix), vmax=np.max(value_matrix))
        
                # Print the output value on the tile
                x_center = 0
                y_center = 0
                axs[i, j].text(x_center, y_center, f'{value_matrix[i, j]:.2f}', color='white', ha='center', va='center', fontsize=12)
        
                # Set the tick labels for the axes
                if i == ksize-1:  # Bottom row
                    axs[i, j].set_xticks([0])
                    axs[i, j].set_xticklabels([variable2_values[j]])
                else:
                    axs[i, j].set_xticks([])
                if j == 0:  # Leftmost column
                    axs[i, j].set_yticks([0])
                    axs[i, j].set_yticklabels([variable1_values[i]])
                else:
                    axs[i, j].set_yticks([])
        
        # Add text descriptions for variable 1
        axs[-1, 1].text(-2.1, -1.5, '$k_{on}$', color='black', ha='center', va='center', fontsize=12, rotation=90)
        # Add text descriptions for variable 2
        axs[1, 0].text(1.2, 2.1, '$t_{on}$', color='black', ha='center', va='center', fontsize=12)
        
        # Add a color bar
        cbar = fig.colorbar(im, ax=axs.ravel().tolist(), orientation='vertical')
        cbar.set_label(value_name)
        
        plt.show()

    def tiles_maxpopsize(self):
        meanpops = self.mean_popsizes()
        valmat = []
        for i, k in enumerate(self.k_list):
            valmat.append([])
            for j, t in enumerate(self.t_list):
                valmat[i].append(max(meanpops[i][j][0]))
        
        valmat = np.array(valmat)   
        self.plot_tiles_with_value(valmat, 'maximum population size')

    def tiles_maxpoptime(self):
        meanpops = self.mean_popsizes()
        valmat = []
        for i, k in enumerate(self.k_list):
            valmat.append([])
            for j, t in enumerate(self.t_list):
                #for sim in self.simulist:
                    #if sim in self.simulist_ksort[i] and sim in self.simulist_tsort[j]:
                        #maxpop = max(sim.full_popsize())
                        #tstep_maxpop = sim.full_popsize().index(maxpop)
                        #valmat[i].append(tstep_maxpop*self.timestep)
                maxpop = max(meanpops[i][j][0])
                tstep_maxpop = np.where(meanpops[i][j][0] == maxpop)[0][0] # meanpops[i][j][0].index(maxpop)
                valmat[i].append(tstep_maxpop*self.timestep)
                        
        for i, ls in enumerate(valmat):
            print(i, len(ls))
            
        valmat = np.array(valmat)
        self.plot_tiles_with_value(valmat, 'time to reach maximum population size [h]')

    def tiles_allexplored(self):
        meanexpl = self.mean_exploredter()
        valmat = []
        for i, k in enumerate(self.k_list):
            valmat.append([])
            for j, t in enumerate(self.t_list):
                expl_index = next((i for i, v in enumerate(meanexpl[i][j][0]) if v == 1), None)
                valmat[i].append(expl_index*self.timestep)

        valmat = np.array(valmat)   
        self.plot_tiles_with_value(valmat, 'time to explore full territory [h]')

    def tiles_v_exploration(self):
        meanexpl = self.mean_exploredter()
        valmat = []
        for i, k in enumerate(self.k_list):
            valmat.append([])
            for j, t in enumerate(self.t_list):
                expl = meanexpl[i][j][0]
                t_0 = list(expl).index([x for x in expl if x > expl[0]][0])
                if expl[-1] == 1:
                    t_end = next((i for i, v in enumerate(meanexpl[i][j][0]) if v == 1), None)
                    vexpl = (1-expl[0])/((t_end-t_0)*self.timestep)
                else:
                    t_end = len(expl)-1
                    vexpl = (expl[-1]-expl[0])/((t_end-t_0)*self.timestep)
                    #print(vexpl)
                valmat[i].append(100*vexpl)

        valmat = np.array(valmat)  
        self.plot_tiles_with_value(valmat, 'exploration velocity [%/h]')
        
# %% plot a single competitive growth simulation

name = 'type the condition_name of your sim here'   
competitive_growth_sim = combination(name, combi_nb=0)
print(competitive_growth_sim)

# adjust resolution of plots with 'y_grain'

# distribution of food
competitive_growth_sim.plot_foodmap()
# distribution of species 1 and 2
competitive_growth_sim.spec1sim.plot_populationmap()
competitive_growth_sim.spec2sim.plot_populationmap()
# prominence of species 1 or 2
competitive_growth_sim.plot_populationmap(method='relative')
# expressing fractions of species 1 and 2
competitive_growth_sim.spec1sim.plot_expressionmap()
competitive_growth_sim.spec2sim.plot_expressionmap()
# the size of the population alive
competitive_growth_sim.plot_popsizes()

# %% visualize the results of a range of simulations
# in this case: environment the same, we look at how different choices of constant \mu and k_on affect the fitness

name = 'type the condition_name of your sim here'   
sims_results = concurring_species(name)

# plot the fraction of species 1 on the cumulative population, as a measure of fitness
sims_results.plot_fullbacfractiles()

# %% get the fitness of a species for a single simulation (with replicates)

name = 'type the condition_name of your sim here'  # multiple entries for names possible as well if necessary 
single_sim_spec1 = condition(name, 'spec1')
single_sim_spec2 = condition(name, 'spec2')

# average number of bacteria having lived until the specified timepoint t_max for each species
cp1 = single_sim_spec1.all_bact_mean(t_max='last')
cp2 = single_sim_spec2.all_bact_mean(t_max='last')

print('fraction of cumulative population: ', cp1/(cp1+cp2))
