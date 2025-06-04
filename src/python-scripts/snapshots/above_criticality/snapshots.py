import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import glob
from matplotlib import rcParams
import matplotlib.ticker as ticker
from matplotlib.colors import hsv_to_rgb


# The function below takes the spreading rate of prey (p), the spontaneous conversion rate (lambda_rate),
# and the lattice size L for an L x L square lattice as input and simulates a single realization of the 
# Chase-Escape process on the lattice. As an output, it returns the variable "escaped"
# which is True if the prey manage to escape and False if all the prey die before escaping.
# It also saves a snapshot of the final lattice state as an image.

def single_realization_simulation(p, lambda_rate, L):
    # Define constants for site colors
    # 0 = white (Empty site)
    # 1 = blue (Predators)
    # 2 = red (Prey)
    
    # Initialize the lattice with all empty sites
    sites = np.zeros((L, L), dtype=int)
    
    # Initialize lists to store BR and RE bonds
    BR_bonds = [] # List to store BR bonds
    RE_bonds = [] # List to store RE bonds
    red_sites = [] # List to store all red site coordinates
  

    # NEW INITIAL CONDITIONS: Single red site in the middle
    middle_i = L // 2
    middle_j = L // 2
    sites[middle_i, middle_j] = 2  # Place one red site in the middle
    red_sites.append((middle_i, middle_j))  # Add to red sites list

    # Create initial RE bonds from the middle red site to its empty neighbors
    up = ((middle_i - 1) % L, middle_j)
    down = ((middle_i + 1) % L, middle_j)
    left = (middle_i, (middle_j - 1) % L)
    right = (middle_i, (middle_j + 1) % L)

    # Add RE bonds to all neighboring empty sites
    RE_bonds.append([(middle_i, middle_j), up])
    RE_bonds.append([(middle_i, middle_j), down])
    RE_bonds.append([(middle_i, middle_j), left])
    RE_bonds.append([(middle_i, middle_j), right])

    # Main simulation loop. The program exits this loop when all the prey die
    # or when any of the prey particles reaches the rightmost column of the lattice (escape!)
    
    iteration = 0
    escaped = 'unknown'
    converted_red_x = []
    converted_red_y = []
    while escaped == 'unknown':
        rnd = np.random.uniform()
        num_BR_bonds = len(BR_bonds)
        num_RE_bonds = len(RE_bonds)
        num_red_sites = len(red_sites)
        
        # Check termination condition - no red sites left
        if num_red_sites == 0:
            escaped = False
            save_lattice_snapshot(sites, iteration, p, lambda_rate, L, escaped, converted_red_x, converted_red_y)
            continue
            
        # Calculate total rate and probabilities
        total_rate = num_BR_bonds + p * num_RE_bonds + lambda_rate * num_red_sites
        
        # Calculate cumulative probabilities for event selection
        prob_BR = num_BR_bonds / total_rate
        prob_RE = (p * num_RE_bonds) / total_rate
        prob_spontaneous = (lambda_rate * num_red_sites) / total_rate
        
        # Select event type based on random number
        if rnd < prob_BR:
            # Choose a BR bond
            index = np.random.randint(0, num_BR_bonds)
            bond = BR_bonds.pop(index)
            i, j = bond[1]
            sites[i, j] = 1 # the red site in the BR bond is now blue
            
            # Remove this site from red_sites list
            red_sites.remove((i, j))
            
            # identify the neighbors of the newly turned blue site (i,j)
            up = ((i - 1) % L, j)
            down = ((i + 1) % L, j)
            left = (i, (j - 1))
            right = (i, (j + 1))
            # make a list of all the neighbors
            neighbors = [up, down, left, right]

            # loop over all the neighbors to update the list of active bonds
            for nbr_index in neighbors:

                # If the neighbor is red, we need to add a new BR bond.
                if sites[nbr_index] == 2:
                    BR_bond = [bond[1], nbr_index]
                    BR_bonds.append(BR_bond)

                # If the neighbor is blue, then previously there was an active BR bond 
                # which now has to be removed from the list of active bonds.
                elif sites[nbr_index] == 1:
                    for BR_bond in BR_bonds:
                        if BR_bond[1] == (i, j) and BR_bond[0] == nbr_index:
                            BR_bonds.remove(BR_bond)
                            break
        
                # If the neighbor is empty, then previously there was an active RE bond
                # which now has to be removed from the list of active bonds.
                elif sites[nbr_index] == 0:
                    for RE_bond in RE_bonds:
                        if RE_bond[0] == (i, j) and RE_bond[1] == nbr_index:
                            RE_bonds.remove(RE_bond)
                            break

        elif rnd < prob_BR + prob_RE:
            # Choose a RE bond
            index = np.random.randint(0, num_RE_bonds)
            bond = RE_bonds.pop(index)
            i, j = bond[1]
            sites[i, j] = 2 # the empty site in the RE bond is now red
            
            # Add this new red site to red_sites list
            red_sites.append((i, j))

            # check escape condition -- if the red site has reached the rightmost column of the
            # square lattice
            if j == L - 1 or j == 0 or i == L - 1 or i == 0:
                escaped = True
                save_lattice_snapshot(sites, iteration, p, lambda_rate, L, escaped, converted_red_x, converted_red_y)
                continue

            # identify the neighbors of the newly turned red site (i,j)    
            up = ((i - 1) % L, j)
            down = ((i + 1) % L, j)
            left = (i, (j - 1))
            right = (i, (j + 1))
            neighbors = [up, down, left, right]

            # loop over all the neighbors to update the list of active bonds
            for nbr_index in neighbors:

                # If the neighbor is red, then previously there was an active RE bond 
                # which now has to be removed.
                if sites[nbr_index] == 2:
                    for RE_bond in RE_bonds:
                        if RE_bond[1] == (i, j) and RE_bond[0] == nbr_index:
                            RE_bonds.remove(RE_bond)
                            break

                # If the neighbor is blue, then we need to add a new BR bond.
                elif sites[nbr_index] == 1:
                    BR_bond = [nbr_index, bond[1]]
                    BR_bonds.append(BR_bond)

                # If the neighbor is empty, then we need to add a new RE bond.
                elif sites[nbr_index] == 0:
                    RE_bond = [bond[1], nbr_index]
                    RE_bonds.append(RE_bond)
                    
        else:
            # Spontaneous conversion: choose a random red site to convert to blue
            index = np.random.randint(0, num_red_sites)
            red_site = red_sites.pop(index)
            i, j = red_site
            sites[i, j] = 1 # convert red site to blue
            converted_red_y.append(i)
            converted_red_x.append(j)
            
            # identify the neighbors of the newly turned blue site (i,j)
            up = ((i - 1) % L, j)
            down = ((i + 1) % L, j)
            left = (i, (j - 1))
            right = (i, (j + 1))
            neighbors = [up, down, left, right]

            # loop over all the neighbors to update the list of active bonds
            for nbr_index in neighbors:

                # If the neighbor is red, we need to add a new BR bond.
                if sites[nbr_index] == 2:
                    BR_bond = [(i, j), nbr_index]
                    BR_bonds.append(BR_bond)

                # If the neighbor is blue, then previously there was an active BR bond 
                # which now has to be removed from the list of active bonds.
                elif sites[nbr_index] == 1:
                    for BR_bond in BR_bonds:
                        if BR_bond[1] == (i, j) and BR_bond[0] == nbr_index:
                            BR_bonds.remove(BR_bond)
                            break
        
                # If the neighbor is empty, then previously there was an active RE bond
                # which now has to be removed from the list of active bonds.
                elif sites[nbr_index] == 0:
                    for RE_bond in RE_bonds:
                        if RE_bond[0] == (i, j) and RE_bond[1] == nbr_index:
                            RE_bonds.remove(RE_bond)
                            break

        iteration += 1
        if iteration % 1000 == 0 or escaped == True:
            save_lattice_snapshot(sites, iteration, p, lambda_rate, L, escaped, converted_red_x, converted_red_y)

    return escaped

def save_lattice_snapshot(sites, iteration, p, lambda_rate, L, escaped, converted_red_x, converted_red_y):

    os.environ['PATH'] += ':/Library/TeX/texbin'
    rcParams['text.usetex'] = True
    rcParams['font.family'] = "Times New Roman"

    rcParams['xtick.labelsize'] = 15
    rcParams['ytick.labelsize'] = 15


    rcParams['xtick.major.width'] = 1
    rcParams['xtick.minor.width'] = 1
    rcParams['xtick.major.size'] = 4
    rcParams['xtick.minor.size'] = 0

    rcParams['ytick.major.width'] = 1
    rcParams['ytick.minor.width'] = 1
    rcParams['ytick.major.size'] = 4
    rcParams['ytick.minor.size'] = 0

    rcParams['ytick.right'] = True
    rcParams['ytick.labelright'] = False
    rcParams['xtick.top'] = True
    rcParams['xtick.labeltop'] = False

    rcParams['axes.labelsize'] = 15
    rcParams['axes.edgecolor'] = 'black'
    rcParams['axes.linewidth'] = 1.1
    rcParams['axes.xmargin'] = 0.01
    rcParams['axes.spines.top'] = True
    rcParams["axes.axisbelow"]  = False
    rcParams["boxplot.meanline"] = True 
    rcParams['legend.frameon'] = False
    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.direction'] = 'in'
    # Create figure with better styling
    plt.figure(figsize=(12, 12))
    
    # Define size of the markers
    predator_size = 2
    prey_size = 2
    
    # Define colors
    predator_color = 'blue'  # Blue
    prey_color = 'red'      # Red
    
    # Get coordinates for each site type
    empty_sites = np.where(sites == 0)
    predator_sites = np.where(sites == 1)
    prey_sites = np.where(sites == 2)
    
    if len(predator_sites[0]) > 0:
        plt.plot(predator_sites[1], predator_sites[0], 'o', color=predator_color, markersize=predator_size)
    
    if len(prey_sites[0]) > 0:
        plt.plot(prey_sites[1], prey_sites[0], 'o', color=prey_color, markersize=prey_size)

    if len(converted_red_x) > 0:
        plt.plot(converted_red_x, converted_red_y, 'o', color='yellow', markersize=prey_size)
    
    # Set up the plot
    plt.xlim(0, L)
    plt.ylim(0, L)
    plt.gca().set_aspect('equal')
    

    # Labels and title
    plt.xlabel('Column', fontsize=12)
    plt.ylabel('Row', fontsize=12)
    plt.title(f'Chase-Escape Lattice ($p$={p}, $\\lambda$={lambda_rate}, $L$={L})\n'
              f'Iteration: {iteration}, Escaped: {escaped}', fontsize=14, pad=20)
    
    
    plt.tight_layout()
    plt.savefig(f'./lattice_snapshot_l_{L}_p_{p}_lambda_{lambda_rate}_iter_{iteration:06d}.png', 
                dpi=150, bbox_inches='tight')
    plt.close()

# Run single realization with image output
p = 2.5
lambda_rate = 1.0
L = 250

print(f"Running single realization with p = {p}, lambda_rate = {lambda_rate}, L = {L}")
escaped = single_realization_simulation(p, lambda_rate, L)
print(f"Result: Escaped = {escaped}")