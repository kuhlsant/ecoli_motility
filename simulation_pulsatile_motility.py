# %% inititation
import copy
import random
import pickle

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

growth_rate_mean, growth_rate_std = 2, 0.5
growth_rate_reduction_mean = 0.1 # amount by which growth rate is reduced during expression
death_rate, wake_rate = 0.5, 1 # 0.018 # 0.43 d^-1, Schink et al., 2019 'Death Rate of E. coli during starvation...

# %% define classes and functions

class bacterium:
    
    def __init__(self, ID, position, motile=False, t=0, efficient=False, k_variable=False):
        self.ID = ID # in binary
        self.position = position
        self.motile = motile
        self.stationary = False
        self.exprage=-1
        self.expradvance = 0
        self.efficient = efficient # if efficient, we will not consider ID + histories to make the simulation faster (they are only necessary to follow single cell faith)
        self.k_variable = k_variable
        self.position_history = []
        self.motile_history = []
        self.division_tracker = []
        habitat.voxel_list[self.position].bacterialist.append(self)
        self.growth_advance = 0
    
        self.initiate_history()

             
    def __str__(self):        
        return f'Bacterium {self.ID}'

    def initiate_history(self, age='random'):
        if age == 'random': self.age = random.random()
        else: self.age = age

    def division(self, t):
        if self.efficient == False: self.division_tracker.append(t)
        # get the new IDs
        if self.efficient: ID_1, ID_2 = 0, 0
        else: ID_1, ID_2 = 2*self.ID, 2*self.ID + 1
        
        # generate offspring
        offspring_1, offspring_2 = copy.deepcopy(self), copy.deepcopy(self)
        offspring_1.initiate_history(age=0)
        offspring_2.initiate_history(age=0)
        offspring_1.ID, offspring_2.ID = ID_1, ID_2
        
        # habitat is defined in the simulation timeloop
        habitat.voxel_list[self.position].bacterialist.append(offspring_1)
        habitat.voxel_list[self.position].bacterialist.append(offspring_2)        
        habitat.voxel_list[self.position].bacterialist.remove(self)
        
    def now_stationary(self):
        '''
        switches bacterium to stationary phase
        '''
        self.motile = False
        self.expression_switch(t, off=True)
        self.stationary = True
            
    def timestepper(self, t):  
        if self.efficient == False: self.motile_history.append(self.motile)
                          
        if self.stationary == False:
            food_consumption = 0

            # check if we have to switch:
            # exprage > 0: motile, <0: nonmotile

            rr = random.random()
            if self.k_variable:
                # interpolate on self.growth_advance between 0 and max 
                # max_growth_advance defined in habitat

                # linear
                #self.k_variable = (habitat.max_growth_advance - self.growth_advance)/habitat.max_growth_advance
                
                # sigmoidal
                #var_input = (habitat.max_growth_advance - self.growth_advance)/habitat.max_growth_advance                
                #self.k_variable = 1/(1+np.exp(-10*(var_input-0.5))) # sigmoidal response function

                # exponential
                #var_input = self.growth_advance/habitat.max_growth_advance          
                #self.k_variable = np.exp(-5*var_input) # exp response function

                # linear up to 0.5 lambda_max
                var_input = (habitat.max_growth_advance - self.growth_advance)/habitat.max_growth_advance
                if var_input > 0.5:
                    self.k_variable = (habitat.max_growth_advance - 2*self.growth_advance)/habitat.max_growth_advance
                else:
                    self.k_variable = 0

                p_on_var = self.k_variable*timestep

                if rr < p_on_var:
                    self.expression_switch(t, off=False) 
                    self.motile = True
                
            else:
                if rr < p_on: 
                    self.expression_switch(t, off=False) 
                    self.motile = True
            
            if self.exprage*(self.exprage+self.expradvance) < 0 or self.exprage*(self.exprage+self.expradvance) == 0: # crossing 0, we are entering off state
                self.expression_switch(t, off=True)
                self.motile = False
    
            elif self.exprage > 0: self.exprage += self.expradvance
            
            if self.efficient == False: self.position_history.append(self.position)
    
            # check if division is due: from the division time it is calculated how many percent of growth until division the bacterium advanced in one timestep
            # growth only happens when there is space for at least one bacterium more
            self.growth_advance = 0
            
            growth_rate = np.random.normal(growth_rate_mean, growth_rate_std)
            available_food = food_per_vox[self.position]
            growth_rate *= available_food/(F_halfmax+available_food)
                            
            if self.exprage>0: 
                absolute_reduction = growth_rate_reduction_mean*growth_rate_mean
                growth_rate -= absolute_reduction # reduce growth rate by some absolute amount for motility

            calculate_doubling = True # only False if we dont have growth

            if growth_rate < 0: 
                self.growth_advance = 0
                calculate_doubling = False

            if self.stationary == False and calculate_doubling:
                doubling_time = np.log(2)/growth_rate
                self.growth_advance = timestep/doubling_time
                self.age += self.growth_advance

            food_consumption = self.growth_advance  # first term: consumption for motility included
            if self.motile: food_consumption += timestep*growth_rate_reduction_mean*growth_rate_mean/np.log(2)
                
            food_consumption += consumption_per_t # regardless of capacity
            habitat.voxel_list[self.position].food_consumption += food_consumption
    
            if self.age > 1:
                self.division(t)  
        
    def expression_switch(self, t, off):
        if off:
            self.exprage = -1
            
        else:
            # age comes from -1: off -> on
            self.exprage = 1
            self.timer_for_expr = np.random.normal(on_time_mean, on_time_std)
            while self.timer_for_expr < 0: 
                self.timer_for_expr = np.random.normal(on_time_mean, on_time_std)
            self.expradvance = -timestep*1/self.timer_for_expr

class voxel:
    def __init__(self, position, food, initial_pop=0, population_size=0, expressing_fraction=0, stationary_fraction=0):
        
        self.position = position
        self.food = food
        self.food_consumption = 0
        self.bacterialist = []
        self.life = True
        self.nonstat = 0 #len([bac for bac in self.bacterialist if bac.stationary == False])
        
        self.initial_pop = initial_pop
        
        self.population_size, self.expressing_fraction, self.stationary_fraction = initial_pop, expressing_fraction, stationary_fraction
        if self.population_size > 0: self.life = True
        
        world.append(self) # environment.voxel_list.append(self)
        self.initiate_bacteria(initial_pop)
    
    def initiate_bacteria(self, startsize, efficient=False, k_variable=False):
        for i in range(startsize):
            bacterium(i+1, self.position, efficient=efficient, k_variable=k_variable)
        # update info        
        self.population_size, self.expressing_fraction, self.stationary_fraction = self.get_info(t)

        
    def timestepper(self, t): 
        # gather information for direct readout
        self.population_size, self.expressing_fraction, self.stationary_fraction = self.get_info(t)
        self.food_consumption = 0
        self.life = True
        if len(self.bacterialist) == 0:
            self.life == False

        if self.food < 0: 
            self.food = 0
            #print('negative food, t = ', t, ' position: ', self.position)
            
        if self.life:

            # make cells stationary if capacity too small
            to_statprob = habitat.to_stat_probability[self.position]
            if t>np.log(2)/growth_rate_mean and to_statprob > 0: # only after one doubling time
                
                livinglist = [bac for bac in self.bacterialist if bac.stationary == False]
                for bac in livinglist:
                    if random.random() < to_statprob:
                        bac.now_stationary()
             
            # reawake cells / let them die
            stationarylist = [bac for bac in self.bacterialist if bac.stationary == True]

            for bact in stationarylist:
                rr = random.random()
                if rr < death_rate*timestep: 
                    self.bacterialist.remove(bact)
                    habitat.dead_bacteria += 1
                elif rr > death_rate*timestep and rr  < death_rate*timestep + wake_rate*timestep: bact.stationary = False

            # stop if all dead
            if len(self.bacterialist) == 0: 
                #print('all are dead.')
                self.life = False

            # let bacteria timestep
            for bact in self.bacterialist:
                bact.timestepper(t)
                

    def get_info(self, t):
        population_size = len(self.bacterialist)
        
        if population_size < 1: expressing_fraction, stat_fraction = np.nan, np.nan
        else: 
            expressing_fraction = len([bact for bact in self.bacterialist if bact.exprage>0])/population_size
            stat_fraction = len([bac for bac in self.bacterialist if bac.stationary == True])/population_size
        
        return population_size, expressing_fraction, stat_fraction

class environment:
    
    def __init__(self, size, voxel_list=[]):
        '''
        voxel_list (world) has to be created before initiating environment
        '''
        self.size = size
        self.voxel_list = voxel_list
        
        for i in range(size):
            voxel(position=i, food=0)
        
        self.population_sizes, self.expressing_fractions, self.food_distribution, self.stationary_fractions = [[] for i in range(size)], [[] for i in range(size)], [[] for i in range(size)], [[] for i in range(size)]
        self.dead_bacteria = 0
        self.dead_size = []

        self.max_growth_advance = timestep*growth_rate_mean/np.log(2) # needed on single bact level for variable k_on

    def __str__(self):
        return print('test')
        
            
    def timestepper(self, t):
        # timestep in each individual voxel
        for voxel in self.voxel_list:
            voxel.timestepper(t)
            
        # relocate bacteria
        self.relocate_bacteria()
        
        # gather information
        for i, voxel in enumerate(self.voxel_list):
            self.population_sizes[i].append(voxel.population_size)
            self.expressing_fractions[i].append(voxel.expressing_fraction)
            self.stationary_fractions[i].append(voxel.stationary_fraction)
            self.food_distribution[i].append(voxel.food)
        
        self.dead_size.append(self.dead_bacteria)
                   
    def relocate_bacteria(self):
        flags = [0, 3]
        for fl in flags:
            self.relocator(fl)

    def relocator(self, flag_nb):
        '''
        relocate all bacteria with a specific nb of flagella
        '''
        #globe_size = len(self.voxel_list)
        if flag_nb == 0:
            bacterialist_iflagella = [[bact for bact in vxel.bacterialist if bact.exprage == -1] for vxel in self.voxel_list]
        elif flag_nb != 0:
            bacterialist_iflagella = [[bact for bact in vxel.bacterialist if bact.exprage > 0] for vxel in self.voxel_list]

        # initiate list where relocations are stored + count how many bacteria with i flagella there are in each voxel
        boxs_bact = [len(x) for x in bacterialist_iflagella]
        boxs_food = [voxl.food for voxl in self.voxel_list]
        fluxes_bact = [get_bacterial_flux(i, boxs_food, boxs_bact, flag_nb) for i in range(len(boxs_bact))] # from -1->0, 0->1, 1->2, ... reverse direction if < 0
            
        bact_reloc = []
        for j in fluxes_bact:
            to_move = int(j)            
            if j > 0:
                if random.random() < j-int(j): to_move+=1 # e.g. j = 2.3: move 2 bacteria, + 3d one with p=0.3
            elif j < 0:
                if random.random() < int(j)-j: to_move-=1
            bact_reloc.append(to_move)
            

        # relocate
        for i, rel in enumerate(bact_reloc):
            if rel>0:
                for j in range(rel):
                    if bacterialist_iflagella[i-1] != []:
                        bact = random.choice(bacterialist_iflagella[i-1])
                        self.move_bacterium(bact, i)

            elif rel<0:
                for j in range(-rel):
                    if bacterialist_iflagella[i] != []:
                        bact = random.choice(bacterialist_iflagella[i])
                        self.move_bacterium(bact, i-1)
    

                        
    def assign_food(self, position, food):
        voxel = self.voxel_list[position]
        voxel.food += food

    def reset_food(self, position):
        voxel = self.voxel_list[position]
        voxel.food = 0
           
    def assign_bacteria(self, position, bacteria, efficient=False, k_variable=False):
        self.voxel_list[position].initiate_bacteria(bacteria, efficient=efficient, k_variable=k_variable)

    def move_bacterium(self, bacterium, destination):
        '''
        move bacterium from where it currently sits to a new destination
        all written for 1D when position = index of voxel in environment
        '''
        self.voxel_list[bacterium.position].bacterialist.remove(bacterium)
        self.voxel_list[destination].bacterialist.append(bacterium)
        bacterium.position = destination

# independent functions

def get_bacterial_flux(voxel_pos, boxes_food, boxes_bact, flag_nb):
    '''
    voxel pos are shifted by one...but should work out...

    from Ford & Lauffenburger, 1990, all parameters as well
    K_d for sure needs adaptation on the food level chosen
    voxel pos from boxes, calculate flux on the interface between voxels
    '''
    # define constants
    Diff_bact = 0.19 # Ford & Lauffenburger 690 µm^2/s = 6.9*10^-6 cm^2/s against 0.19 µm^2/s, Tavaddod et al., 2011 ????????? Berg & Turner, 1990: 0.01 - 0.07...
    K_d = 0.8 # 0.08 mM
    if flag_nb == 0: v_run = 0
    elif flag_nb > 0: v_run = 22 #17.5+2*np.sqrt(flag_nb) # µm/s, Remy
    sigma_RT = 75 # s
    p_0 = 1 # /s, look up sources for tumbling prob in absence of a gradient
    phi = 0.3 # directional persistence, according to W. Alt, J. Math. Biol., 9, 147 (1980). (check that)
    mu_0 = v_run**2/(p_0*(1-phi)) # 1100 µm^2/s = 1.1*10^-5 cm^2/s
    
    # test
    # get variables
    a, b = boxes_food[voxel_pos-1], boxes_bact[voxel_pos-1]
    da_dx = (boxes_food[voxel_pos]-boxes_food[voxel_pos-1])/box_size
    db_dx = (boxes_bact[voxel_pos]-boxes_bact[voxel_pos-1])/box_size

    # get chemotactic velocity
    V_c = v_run*np.tanh(v_run*sigma_RT*da_dx/(K_d*(1+a/K_d)**2))

    # get flux
    J_bact = -mu_0*db_dx + V_c*b
    J_0 = -Diff_bact*db_dx # just diffusion

    # leads to absolute displacement of that many bacteria:
    J_abs = (J_bact + J_0)*timestep_s/box_size
    

    return J_abs

def get_food_flux(voxel_pos, boxes_food):
    '''
    again: voxels shifted by one...but I accounted for that...I think
    voxel pos from boxes, calculate flux on the interface between voxels
    '''
    # define constants
    Diff_food = 600 # µm^2/s, diffusion coefficient of glucose: https://bionumbers.hms.harvard.edu/files/Diffusion%20coefficients%20of%20various%20substances%20in%20Water.pdf
    
    # get variables
    #a = boxes_food[voxel_pos]
    da_dx = (boxes_food[voxel_pos]-boxes_food[voxel_pos-1])/box_size

    # get flux
    J_food = -Diff_food*da_dx

    # leads to absolute displacement of that many food units:
    J_abs = J_food*timestep_food/box_size

    return J_abs

def food_relocator(habitat):
    '''
    relocate food
    '''
    # initiate list where relocations are stored + count how much food in each voxel
    boxs_food = [voxl.food for voxl in habitat.voxel_list]
    fluxes_food = [get_food_flux(i, boxs_food) for i in range(len(boxs_food))] # from -1->0, 0->1, 1->2, ... reverse direction if < 0
    food_reloc = []
    for j in fluxes_food:
        to_move = j #round(j, 6) #int(j)            
        #if j > 0:
        #    if random.random() < j-int(j): to_move+=1 # e.g. j = 2.3: move 2 food, + 3d one with p=0.3
        #elif j < 0:
        #    if random.random() < int(j)-j: to_move-=1
        food_reloc.append(to_move)
        

    # relocate
    for i, rel in enumerate(food_reloc):
        if rel>0:
            boxs_food[i-1] -= rel
            boxs_food[i] += rel

        elif rel<0:
            rel = np.sqrt(rel*rel) # take absolute value
            boxs_food[i-1] += rel
            boxs_food[i] -= rel

    return boxs_food

# for exporting information later
def exp_full_popsize(habitat):
    fullpopsize = [sum([voxl[t] for voxl in habitat.population_sizes]) for t in range(timerange)]
    return fullpopsize

def exp_full_expressfrac(habitat):
    expressfrac = [np.nansum([voxl[t] for voxl in habitat.expressing_fractions])/habitat.size for t in range(timerange)]
    return expressfrac

def exp_full_statfrac(habitat):
    statfrac = [np.nansum([voxl[t] for voxl in habitat.stationary_fractions])/habitat.size for t in range(timerange)]
    return statfrac

def exp_explored_territory(habitat):
    explored_voxels, fraction_explored = [], []
    for t in range(timerange):
        explored_voxels.append([])
        if t == 0:
            for i in range(habitat.size):
                if habitat.population_sizes[i][t] != 0: explored_voxels[-1].append(1)
                else: explored_voxels[-1].append(0)
        elif len(explored_voxels)>1:
            for i in range(habitat.size):
                if explored_voxels[-2][i] == 1: explored_voxels[-1].append(1)
                elif explored_voxels[-2][i] == 0 and habitat.population_sizes[i][t] != 0: explored_voxels[-1].append(1)
                else: explored_voxels[-1].append(0)
            explored_voxels.remove(explored_voxels[0])
        fraction_explored.append(sum(explored_voxels[-1])/habitat.size)
    return fraction_explored

# %% run simulation (two species)

habitat_size = 50
efficient = True # fast run or slow run with more single cell infos
save_efficiently = True # only save lite, to further increase speed
timestep = 1/1800
box_size = 100 # µm       
timerange = 100000 
food_timescale = 4 # times the food distribution calculation runs faster (has to be int)
growth_rate_reduction_mean = 0.1
F_halfmax = 10 # food concentration at which growth rate is half maximum

condition_name = 'switch02_agmin'
t_list_1, t_list_2 = [1, 1, 1, 2, 2, 2, 3, 3, 3], [0]*9
k_list_1, k_list_2 = [0.05, 0.1, 0.2]*3, [0]*9
p_foodlocswitch = 0.2

repeats = 4

for r in range(repeats):
    print('r = ', r)
    for i in range(len(t_list_1)):   

        combination_nb = i
        timestep_s = timestep*3600
        timestep_food = timestep_s/food_timescale
        
        t = 0
        world_1, world_2 = [], []

        world  = world_1
        habitat_1 = environment(habitat_size, voxel_list = world_1)
        habitat = habitat_1
        #habitat_1.assign_bacteria(45, 4, efficient=efficient)
        for h in range(habitat_size):
            habitat_1.assign_bacteria(h, 1, efficient=efficient)

        world = world_2
        habitat_2 = environment(habitat_size, voxel_list = world_2)
        habitat = habitat_2
        #habitat_2.assign_bacteria(67, 4, efficient=efficient)
        for h in range(habitat_size):
            habitat_2.assign_bacteria(h, 1, efficient=efficient)
        
        # timestep*v_max < boxsize !
        
        food_consumed_per_h_per_bacterium = 2
        consumption_per_t = food_consumed_per_h_per_bacterium*timestep
        
        for k in range(habitat_size):
            habitat_1.assign_food(k, 2)
            habitat_2.assign_food(k, 2)

        #habitat_1.assign_food(10, 1000)
        #habitat_2.assign_food(10, 1000)

        foodconsumption_t = []

        foodloc = int(habitat_size/2)

        for t in tqdm(range(timerange)):

            # get capacities
            # only count nonstationary bacteria

            alive_pop_per_vox = [len([bac for bac in habitat_1.voxel_list[position].bacterialist if bac.stationary == False])+len([bac for bac in habitat_2.voxel_list[position].bacterialist if bac.stationary == False]) for position in range(habitat_size)]           
            food_per_vox = [voxl.food for voxl in habitat_2.voxel_list]
            # bacteria to stationary if amount of food can't account for consumption_per_t for all bacteria
            to_stat_probability = (np.array(alive_pop_per_vox)*consumption_per_t-np.array(food_per_vox))/(consumption_per_t*np.array(alive_pop_per_vox)) # prob of finding bacteria exceeding capacity (if capacity == #food)
                
            habitat_1.food_per_vox = food_per_vox
            habitat_2.food_per_vox = food_per_vox
            habitat_1.to_stat_probability = to_stat_probability
            habitat_2.to_stat_probability = to_stat_probability
            
            food_consumed = np.zeros(habitat_size)

            # timestepper bacteria

            on_time_mean, on_time_std, p_on, habitat = t_list_2[i], t_list_2[i]/2, k_list_2[i]*timestep, habitat_2
            habitat_2.timestepper(t)
            food_consumed += np.array([voxl.food_consumption for voxl in habitat_2.voxel_list])

            on_time_mean, on_time_std, p_on, habitat = t_list_1[i], t_list_1[i]/2, k_list_1[i]*timestep, habitat_1
            habitat_1.timestepper(t)
            food_consumed += np.array([voxl.food_consumption for voxl in habitat_1.voxel_list])

            foodconsumption_t.append(sum(food_consumed))

            # manage food consumption
            for k in range(habitat_size):
                habitat_1.assign_food(k, -food_consumed[k])
                habitat_2.assign_food(k, -food_consumed[k])

            for foodz in range(int(food_timescale)):
                # relocate food (diffusion)
                new_food_distribution = food_relocator(habitat_2)
                for k in range(habitat_size):
                    habitat_1.reset_food(k)
                    habitat_2.reset_food(k)
                    habitat_1.assign_food(k, new_food_distribution[k])
                    habitat_2.assign_food(k, new_food_distribution[k])

            # fixed bead emitting food
            #foodlocs = [12, 36]

            # bead emits new food. with prob of 5% the bead will move
            #if random.random() < 0.1:
            #    rrvar = random.random()
            #    if rrvar < 0.5:
            #        foodloc -= 1
            #    elif rrvar > 0.5:
            #        foodloc += 1

            # 2 beads, emitting bead switches
            #if t == 0: beadchoice = random.choice([0, 1])
            #if random.random() < p_foodlocswitch*timestep:
            #    if beadchoice == 0: beadchoice = 1
            #    elif beadchoice == 1: beadchoice = 0

            # bead moves to random place
            if random.random() < p_foodlocswitch*timestep:
                foodloc = random.randrange(0, habitat_size)

            #if foodloc < 0: foodloc = habitat_size+foodloc
            #if foodloc >= habitat_size: foodloc = foodloc-habitat_size

            habitat_1.assign_food(foodloc, 1) #0.02
            habitat_2.assign_food(foodloc, 1)
                
            # with prob of x% assign food to random spot, assign ~1800 food /h
            #if random.random() < 0.002: # 0.006 ~ 5 times per h
            #    foodloc = random.randint(0, habitat_size-1)
            #    habitat_2.assign_food(foodloc, 1500)

            #for h in range(habitat_size):
            #    habitat_1.assign_food(h, 0.01)
            #    habitat_2.assign_food(h, 0.01)

        on_time_mean, on_time_std, k_on = t_list_1[i], t_list_1[i]/2, k_list_1[i]
        data_export_spec1, data_export_lite_spec1 = {}, {}

        data_export_lite_spec1.update({'fullpopsize': exp_full_popsize(habitat_1), 'expressfrac': exp_full_expressfrac(habitat_1), 'statfrac': exp_full_statfrac(habitat_1), 'explored_territory': exp_explored_territory(habitat_1)})
        data_export_lite_spec1.update({'growth_rate_mean': growth_rate_mean, 'growth_rate_std': growth_rate_std, 'growth_rate_reduction_mean': growth_rate_reduction_mean})
        data_export_lite_spec1.update({'on_time_mean': on_time_mean, 'on_time_std': on_time_std, 'k_on': k_on, 'size': habitat_1.size, 'timerange': timerange})           

        name_1 = condition_name+'_spec1_combi'+str(combination_nb)+'_full_t'+str(on_time_mean)+'_k'+str(k_on)+'_r'+str(r)+'.pkl'
        name_lite_1 = condition_name+'_spec1_combi'+str(combination_nb)+'_lite_t'+str(on_time_mean)+'_k'+str(k_on)+'_r'+str(r)+'.pkl'

        if save_efficiently == False:
            data_export_spec1.update({'population_sizes': habitat_1.population_sizes, 'expressing_fractions': habitat_1.expressing_fractions, 'food_distribution': habitat_1.food_distribution, 'stationary_fractions': habitat_1.stationary_fractions}) #, 'food_flux': habitat_1.food_fluxes})
            data_export_spec1.update({'growth_rate_mean': growth_rate_mean, 'growth_rate_std': growth_rate_std, 'growth_rate_reduction_mean': growth_rate_reduction_mean})
            data_export_spec1.update({'on_time_mean': on_time_mean, 'on_time_std': on_time_std, 'k_on': k_on, 'size': habitat_1.size, 'timerange': timerange})
            with open(name_1, 'wb') as file:
                pickle.dump(data_export_spec1, file)

        if efficient == False:
            motilehist = [bact.motile_history for i in range(habitat_size) for bact in habitat_1.voxel_list[i].bacterialist]
            data_export_spec1.update({'motile_history': motilehist})


        with open(name_lite_1, 'wb') as file:
            pickle.dump(data_export_lite_spec1, file)
            
        on_time_mean, on_time_std, k_on = t_list_2[i], t_list_2[i]/2, k_list_2[i]
        data_export_spec2, data_export_lite_spec2 = {}, {}

        data_export_lite_spec2.update({'fullpopsize': exp_full_popsize(habitat_2), 'expressfrac': exp_full_expressfrac(habitat_2), 'statfrac': exp_full_statfrac(habitat_2), 'explored_territory': exp_explored_territory(habitat_2)})
        data_export_lite_spec2.update({'growth_rate_mean': growth_rate_mean, 'growth_rate_std': growth_rate_std, 'growth_rate_reduction_mean': growth_rate_reduction_mean})
        data_export_lite_spec2.update({'on_time_mean': on_time_mean, 'on_time_std': on_time_std, 'k_on': k_on, 'size': habitat_2.size, 'timerange': timerange})           

        name_2 = condition_name+'_spec2_combi'+str(combination_nb)+'_full_t'+str(on_time_mean)+'_k'+str(k_on)+'_r'+str(r)+'.pkl'
        name_lite_2 = condition_name+'_spec2_combi'+str(combination_nb)+'_lite_t'+str(on_time_mean)+'_k'+str(k_on)+'_r'+str(r)+'.pkl'

        if save_efficiently == False:
            data_export_spec2.update({'population_sizes': habitat_2.population_sizes, 'expressing_fractions': habitat_2.expressing_fractions, 'food_distribution': habitat_2.food_distribution, 'stationary_fractions': habitat_2.stationary_fractions}) #, 'food_flux': habitat_2.food_fluxes})
            data_export_spec2.update({'growth_rate_mean': growth_rate_mean, 'growth_rate_std': growth_rate_std, 'growth_rate_reduction_mean': growth_rate_reduction_mean})
            data_export_spec2.update({'on_time_mean': on_time_mean, 'on_time_std': on_time_std, 'k_on': k_on, 'size': habitat_2.size, 'timerange': timerange})
            with open(name_2, 'wb') as file:
                pickle.dump(data_export_spec2, file)
        if efficient == False:
            motilehist = [bact.motile_history for i in range(habitat_size) for bact in habitat_2.voxel_list[i].bacterialist]
            data_export_spec2.update({'motile_history': motilehist})


        with open(name_lite_2, 'wb') as file:
            pickle.dump(data_export_lite_spec2, file)
            
        print('saved ', name_2)

condition_name = 'switch02_agmax'
t_list_1, t_list_2 = [1, 1, 1, 2, 2, 2, 3, 3, 3], [10]*9
k_list_1, k_list_2 = [0.05, 0.1, 0.2]*3, [10]*9

for r in range(repeats):
    print('r = ', r)
    for i in range(len(t_list_1)):   

        combination_nb = i
        timestep_s = timestep*3600
        timestep_food = timestep_s/food_timescale
        
        t = 0
        world_1, world_2 = [], []

        world  = world_1
        habitat_1 = environment(habitat_size, voxel_list = world_1)
        habitat = habitat_1
        #habitat_1.assign_bacteria(45, 4, efficient=efficient)
        for h in range(habitat_size):
            habitat_1.assign_bacteria(h, 1, efficient=efficient)

        world = world_2
        habitat_2 = environment(habitat_size, voxel_list = world_2)
        habitat = habitat_2
        #habitat_2.assign_bacteria(67, 4, efficient=efficient)
        for h in range(habitat_size):
            habitat_2.assign_bacteria(h, 1, efficient=efficient)
        
        # timestep*v_max < boxsize !
        
        food_consumed_per_h_per_bacterium = 2
        consumption_per_t = food_consumed_per_h_per_bacterium*timestep
        
        for k in range(habitat_size):
            habitat_1.assign_food(k, 2)
            habitat_2.assign_food(k, 2)

        #habitat_1.assign_food(10, 1000)
        #habitat_2.assign_food(10, 1000)

        foodconsumption_t = []

        foodloc = int(habitat_size/2)

        for t in tqdm(range(timerange)):

            # get capacities
            # only count nonstationary bacteria

            alive_pop_per_vox = [len([bac for bac in habitat_1.voxel_list[position].bacterialist if bac.stationary == False])+len([bac for bac in habitat_2.voxel_list[position].bacterialist if bac.stationary == False]) for position in range(habitat_size)]           
            food_per_vox = [voxl.food for voxl in habitat_2.voxel_list]
            # bacteria to stationary if amount of food can't account for consumption_per_t for all bacteria
            to_stat_probability = (np.array(alive_pop_per_vox)*consumption_per_t-np.array(food_per_vox))/(consumption_per_t*np.array(alive_pop_per_vox)) # prob of finding bacteria exceeding capacity (if capacity == #food)
                
            habitat_1.food_per_vox = food_per_vox
            habitat_2.food_per_vox = food_per_vox
            habitat_1.to_stat_probability = to_stat_probability
            habitat_2.to_stat_probability = to_stat_probability
            
            food_consumed = np.zeros(habitat_size)

            # timestepper bacteria

            on_time_mean, on_time_std, p_on, habitat = t_list_2[i], t_list_2[i]/2, k_list_2[i]*timestep, habitat_2
            habitat_2.timestepper(t)
            food_consumed += np.array([voxl.food_consumption for voxl in habitat_2.voxel_list])

            on_time_mean, on_time_std, p_on, habitat = t_list_1[i], t_list_1[i]/2, k_list_1[i]*timestep, habitat_1
            habitat_1.timestepper(t)
            food_consumed += np.array([voxl.food_consumption for voxl in habitat_1.voxel_list])

            foodconsumption_t.append(sum(food_consumed))

            # manage food consumption
            for k in range(habitat_size):
                habitat_1.assign_food(k, -food_consumed[k])
                habitat_2.assign_food(k, -food_consumed[k])

            for foodz in range(int(food_timescale)):
                # relocate food (diffusion)
                new_food_distribution = food_relocator(habitat_2)
                for k in range(habitat_size):
                    habitat_1.reset_food(k)
                    habitat_2.reset_food(k)
                    habitat_1.assign_food(k, new_food_distribution[k])
                    habitat_2.assign_food(k, new_food_distribution[k])

            # fixed bead emitting food
            #foodlocs = [12, 36]

            # bead emits new food. with prob of 5% the bead will move
            #if random.random() < 0.1:
            #    rrvar = random.random()
            #    if rrvar < 0.5:
            #        foodloc -= 1
            #    elif rrvar > 0.5:
            #        foodloc += 1

            # 2 beads, emitting bead switches
            #if t == 0: beadchoice = random.choice([0, 1])
            #if random.random() < p_foodlocswitch*timestep:
            #    if beadchoice == 0: beadchoice = 1
            #    elif beadchoice == 1: beadchoice = 0

            # bead moves to random place
            if random.random() < p_foodlocswitch*timestep:
                foodloc = random.randrange(0, habitat_size)

            #if foodloc < 0: foodloc = habitat_size+foodloc
            #if foodloc >= habitat_size: foodloc = foodloc-habitat_size

            habitat_1.assign_food(foodloc, 1) #0.02
            habitat_2.assign_food(foodloc, 1)
                
            # with prob of x% assign food to random spot, assign ~1800 food /h
            #if random.random() < 0.002: # 0.006 ~ 5 times per h
            #    foodloc = random.randint(0, habitat_size-1)
            #    habitat_2.assign_food(foodloc, 1500)

            #for h in range(habitat_size):
            #    habitat_1.assign_food(h, 0.01)
            #    habitat_2.assign_food(h, 0.01)

        on_time_mean, on_time_std, k_on = t_list_1[i], t_list_1[i]/2, k_list_1[i]
        data_export_spec1, data_export_lite_spec1 = {}, {}

        data_export_lite_spec1.update({'fullpopsize': exp_full_popsize(habitat_1), 'expressfrac': exp_full_expressfrac(habitat_1), 'statfrac': exp_full_statfrac(habitat_1), 'explored_territory': exp_explored_territory(habitat_1)})
        data_export_lite_spec1.update({'growth_rate_mean': growth_rate_mean, 'growth_rate_std': growth_rate_std, 'growth_rate_reduction_mean': growth_rate_reduction_mean})
        data_export_lite_spec1.update({'on_time_mean': on_time_mean, 'on_time_std': on_time_std, 'k_on': k_on, 'size': habitat_1.size, 'timerange': timerange})           

        name_1 = condition_name+'_spec1_combi'+str(combination_nb)+'_full_t'+str(on_time_mean)+'_k'+str(k_on)+'_r'+str(r)+'.pkl'
        name_lite_1 = condition_name+'_spec1_combi'+str(combination_nb)+'_lite_t'+str(on_time_mean)+'_k'+str(k_on)+'_r'+str(r)+'.pkl'

        if save_efficiently == False:
            data_export_spec1.update({'population_sizes': habitat_1.population_sizes, 'expressing_fractions': habitat_1.expressing_fractions, 'food_distribution': habitat_1.food_distribution, 'stationary_fractions': habitat_1.stationary_fractions}) #, 'food_flux': habitat_1.food_fluxes})
            data_export_spec1.update({'growth_rate_mean': growth_rate_mean, 'growth_rate_std': growth_rate_std, 'growth_rate_reduction_mean': growth_rate_reduction_mean})
            data_export_spec1.update({'on_time_mean': on_time_mean, 'on_time_std': on_time_std, 'k_on': k_on, 'size': habitat_1.size, 'timerange': timerange})
            with open(name_1, 'wb') as file:
                pickle.dump(data_export_spec1, file)

        if efficient == False:
            motilehist = [bact.motile_history for i in range(habitat_size) for bact in habitat_1.voxel_list[i].bacterialist]
            data_export_spec1.update({'motile_history': motilehist})


        with open(name_lite_1, 'wb') as file:
            pickle.dump(data_export_lite_spec1, file)
            
        on_time_mean, on_time_std, k_on = t_list_2[i], t_list_2[i]/2, k_list_2[i]
        data_export_spec2, data_export_lite_spec2 = {}, {}

        data_export_lite_spec2.update({'fullpopsize': exp_full_popsize(habitat_2), 'expressfrac': exp_full_expressfrac(habitat_2), 'statfrac': exp_full_statfrac(habitat_2), 'explored_territory': exp_explored_territory(habitat_2)})
        data_export_lite_spec2.update({'growth_rate_mean': growth_rate_mean, 'growth_rate_std': growth_rate_std, 'growth_rate_reduction_mean': growth_rate_reduction_mean})
        data_export_lite_spec2.update({'on_time_mean': on_time_mean, 'on_time_std': on_time_std, 'k_on': k_on, 'size': habitat_2.size, 'timerange': timerange})           

        name_2 = condition_name+'_spec2_combi'+str(combination_nb)+'_full_t'+str(on_time_mean)+'_k'+str(k_on)+'_r'+str(r)+'.pkl'
        name_lite_2 = condition_name+'_spec2_combi'+str(combination_nb)+'_lite_t'+str(on_time_mean)+'_k'+str(k_on)+'_r'+str(r)+'.pkl'

        if save_efficiently == False:
            data_export_spec2.update({'population_sizes': habitat_2.population_sizes, 'expressing_fractions': habitat_2.expressing_fractions, 'food_distribution': habitat_2.food_distribution, 'stationary_fractions': habitat_2.stationary_fractions}) #, 'food_flux': habitat_2.food_fluxes})
            data_export_spec2.update({'growth_rate_mean': growth_rate_mean, 'growth_rate_std': growth_rate_std, 'growth_rate_reduction_mean': growth_rate_reduction_mean})
            data_export_spec2.update({'on_time_mean': on_time_mean, 'on_time_std': on_time_std, 'k_on': k_on, 'size': habitat_2.size, 'timerange': timerange})
            with open(name_2, 'wb') as file:
                pickle.dump(data_export_spec2, file)
        if efficient == False:
            motilehist = [bact.motile_history for i in range(habitat_size) for bact in habitat_2.voxel_list[i].bacterialist]
            data_export_spec2.update({'motile_history': motilehist})


        with open(name_lite_2, 'wb') as file:
            pickle.dump(data_export_lite_spec2, file)
            
        print('saved ', name_2)






# %% code for plotting

def plot_population_size_thr(y_grain, method='absolute', threshold=None, above_threshold_color='red'):
    x_len = habitat.size  # number of tiles in the x direction
    # y_len = 60   # number of timesteps in the y direction
    # create a data array
    # 1) get timestep for readout
    
    big_t_step = int((timerange-1)//y_grain)
    values = []
    for i in range(y_grain):
        values.append([])
        for j in range(x_len):
            values[i].append(habitat.population_sizes[j][i*big_t_step+1])
            
    for i in range(y_grain):
        for j in range(x_len):
            if values[i][j] == 0: values[i][j] = np.nan
    
    values = np.array(values)
    
    # Applying threshold and color mapping
    if threshold is not None:
        values[values > threshold] = threshold
    
    cmap = 'viridis'
    if above_threshold_color is not None:
        cmap = plt.cm.get_cmap(cmap)
        cmap.set_over(above_threshold_color)
    
    # Plotting
    plt.imshow(values, cmap=cmap, aspect='auto', interpolation='none', vmax=threshold)
    if threshold is not None:
        plt.colorbar(label='Population Size (cutoff at {})'.format(threshold))
    else:
        plt.colorbar(label='Population Size')
    
    plt.xlabel('X [mm]')
    plt.ylabel('h')
    plt.xticks(np.arange(0, x_len, step=5), [np.round(x, 2) for x in np.arange(0, x_len*box_size/1000, step=5*box_size/1000)])
    plt.yticks(np.arange(0, y_grain, step=9), [np.round(x,2) for x in np.arange(0, y_grain*big_t_step*timestep, step=9*big_t_step*timestep)])
    
    plt.show()


def plot_population_size(y_grain, method='absolute'):
    x_len = habitat.size  # number of tiles in the x direction
    # y_len = 60   # number of timesteps in the y direction
    # create a data array
    # 1) get timestep for readout
    
    big_t_step = int((timerange-1)//y_grain)
    values = []
    for i in range(y_grain):
        values.append([])
        for j in range(x_len):
            values[i].append(habitat.population_sizes[j][i*big_t_step+1])
            
    for i in range(y_grain):
        for j in range(x_len):
            if values[i][j] == 0: values[i][j] = np.nan
            
    if method == 'relative':
        for i, t_list in enumerate(values):
            values[i] = [np.nan if np.nansum(t_list) == 0 else val/np.nansum(t_list) for val in t_list]
    
    values = np.array(values)
    
    # Plotting
    plt.imshow(values, cmap='viridis', aspect='auto', interpolation='none')
    plt.colorbar(label='Population Size')
    
    plt.xlabel('X [mm]')
    plt.ylabel('h')
    plt.xticks(np.arange(0, x_len, step=5), [np.round(x, 2) for x in np.arange(0, x_len*box_size/1000, step=5*box_size/1000)])
    plt.yticks(np.arange(0, y_grain, step=9), [np.round(x,2) for x in np.arange(0, y_grain*big_t_step*timestep, step=9*big_t_step*timestep)])
    
    plt.show()

def plot_food(y_grain):
    x_len = habitat.size  # number of tiles in the x direction
    # y_len = 60   # number of timesteps in the y direction
    
    # create a data array
    # 1) get timestep for readout
    
    big_t_step = int((timerange-1)//y_grain)
    values_food = []
    for i in range(y_grain):
        values_food.append([])
        for j in range(x_len):
            values_food[i].append(habitat.food_distribution[j][i*big_t_step+1])
    
    values_food = np.array(values_food)
    
    plt.imshow(values_food, cmap='magma', aspect='auto', interpolation='none')
    plt.colorbar(label='available food units')
    
    plt.xlabel('X [mm]')
    plt.ylabel('h')
    plt.xticks(np.arange(0, x_len, step=5), [np.round(x, 2) for x in np.arange(0, x_len*box_size/1000, step=5*box_size/1000)])
    plt.yticks(np.arange(0, y_grain, step=9), [np.round(x,2) for x in np.arange(0, y_grain*big_t_step*timestep, step=9*big_t_step*timestep)])
    
    plt.show()

def plot_expressing_fractions(y_grain):
    x_len = habitat.size  # number of tiles in the x direction
    
    # create a data array
    # 1) get timestep for readout
    
    big_t_step = int((timerange-1)//y_grain)
    values_expr = []
    for i in range(y_grain):
        values_expr.append([])
        for j in range(x_len):
            values_expr[i].append(habitat.expressing_fractions[j][i*big_t_step+1])
    
    values_expr = np.array(values_expr)
    
    plt.imshow(values_expr, cmap='viridis', aspect='auto', interpolation='none')
    plt.colorbar(label='fraction of expressing bacteria')
    
    plt.xlabel('X [mm]')
    plt.ylabel('h')
    plt.xticks(np.arange(0, x_len, step=5), [np.round(x, 2) for x in np.arange(0, x_len*box_size/1000, step=5*box_size/1000)])
    plt.yticks(np.arange(0, y_grain, step=9), [np.round(x,2) for x in np.arange(0, y_grain*big_t_step*timestep, step=9*big_t_step*timestep)])
    
    plt.show()
    
def plot_stationary_fractions(y_grain):
    x_len = habitat.size  # number of tiles in the x direction
    
    # create a data array
    # 1) get timestep for readout
    
    big_t_step = int((timerange-1)//y_grain)
    values_expr = []
    for i in range(y_grain):
        values_expr.append([])
        for j in range(x_len):
            values_expr[i].append(habitat.stationary_fractions[j][i*big_t_step+1])
    
    values_expr = np.array(values_expr)
    
    plt.imshow(values_expr, cmap='viridis', aspect='auto', interpolation='none')
    plt.colorbar(label='fraction of stationary bacteria')
    
    plt.xlabel('X [mm]')
    plt.ylabel('h')
    plt.xticks(np.arange(0, x_len, step=5), [np.round(x, 2) for x in np.arange(0, x_len*box_size/1000, step=5*box_size/1000)])
    plt.yticks(np.arange(0, y_grain, step=9), [np.round(x,2) for x in np.arange(0, y_grain*big_t_step*timestep, step=9*big_t_step*timestep)])
    
    plt.show()
    
# %% plotting
#habitat = habitat_1
print('population size')
#plot_population_size(90, method='relative')
print('food consumption')
plot_food(90)
print('expressing fraction')
plot_expressing_fractions(90)
print('stationary fraction')
plot_stationary_fractions(90)
plot_population_size_thr(90, threshold=20)
# %%% more plotting

fullpopsize = [sum([voxl[t] for voxl in habitat.population_sizes]) for t in range(len(habitat.population_sizes[0]))]
expressfrac = [np.nansum([voxl[t] for voxl in habitat.expressing_fractions])/habitat.size for t in range(len(habitat.population_sizes[0]))]
statfrac = [np.nansum([voxl[t] for voxl in habitat.stationary_fractions])/habitat.size for t in range(len(habitat.population_sizes[0]))]

explored_voxels, fraction_explored = [], []
for t in range(timerange):
    explored_voxels.append([])
    if t == 0:
        for i in range(habitat.size):
            if habitat.population_sizes[i][t] != 0: explored_voxels[-1].append(1)
            else: explored_voxels[-1].append(0)
    elif len(explored_voxels)>1:
        for i in range(habitat.size):
            if explored_voxels[-2][i] == 1: explored_voxels[-1].append(1)
            elif explored_voxels[-2][i] == 0 and habitat.population_sizes[i][t] != 0: explored_voxels[-1].append(1)
            else: explored_voxels[-1].append(0)
        explored_voxels.remove(explored_voxels[0])
    fraction_explored.append(sum(explored_voxels[-1])/habitat.size)
        

fig, axs = plt.subplots(1, 1)
axs.plot(fullpopsize)
axs2 = axs.twinx()
axs2.plot(expressfrac, color='red')
axs2.plot(fraction_explored, color='green')
axs2.plot(statfrac, color='orange')
plt.show()

# %% plot single bact

## bact = random.choice(graveyard) 
bact = random.choice([bbb for vxx in habitat_1.voxel_list for bbb in vxx.bacterialist])
#bact = random.choice([bbb for vxx in habitat_1.voxel_list for bbb in vxx.bacterialist if bbb.ID == 133])
#bact = random.choice([bb for bb in graveyard if bb.flagella>0]
print(bact)
exprhist, divhist = bact.motile_history, bact.division_tracker
posnhist = bact.position_history
time_axis = np.arange(0, len(posnhist)*timestep, timestep)

fig, axs = plt.subplots(1, 1, figsize=(10, 4))
axs.scatter(time_axis, posnhist)
axs.set_yticks(np.arange(min(posnhist), max(posnhist), 5))
axs.set_yticklabels(np.arange(min(posnhist), max(posnhist), 5)/5) # 5 boxes = 1mm

#print(exprhist)
axs02 = axs.twinx()
axs02.plot(np.arange(0, len(exprhist), )*timestep, exprhist, color='red')
axs02.set_yticks([0, 1])
axs02.set_yticklabels(['no', 'yes']) # 5 boxes = 1mm
axs.set_xlabel('time [h]')
axs.set_ylabel('box position')
axs02.set_ylabel('motility')
plt.show()



                 

# %%
