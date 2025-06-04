import numpy as np
import matplotlib.pyplot as plt

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
    
    # ORIGINAL INITIAL CONDITIONS (COMMENTED OUT)
    # Set initial conditions (first column blue, second column red)
    # sites[:, 0] = 1
    # sites[:, 1] = 2

    # Initialize lists to store BR and RE bonds
    BR_bonds = [] # List to store BR bonds
    RE_bonds = [] # List to store RE bonds
    red_sites = [] # List to store all red site coordinates
        
    # ORIGINAL BOND INITIALIZATION (COMMENTED OUT)    
    # Create initial BR and RE bonds along the lattice
    # Also populate the initial red sites list
    # for i in range(L):
    #     BR_bond = [(i, 0), (i, 1)] # BR bond connects first and second column
    #     RE_bond = [(i, 1), (i, 2)] # RE bond connects second and third column
    #     BR_bonds.append(BR_bond)
    #     RE_bonds.append(RE_bond)
    #     red_sites.append((i, 1)) # All sites in column 1 are initially red

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
    while escaped == 'unknown':
        rnd = np.random.uniform()
        num_BR_bonds = len(BR_bonds)
        num_RE_bonds = len(RE_bonds)
        num_red_sites = len(red_sites)
        
        # Check termination condition - no red sites left
        if num_red_sites == 0:
            escaped = False
            save_lattice_snapshot(sites, iteration, p, lambda_rate, L, escaped)
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
                save_lattice_snapshot(sites, iteration, p, lambda_rate, L, escaped)
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
            save_lattice_snapshot(sites, iteration, p, lambda_rate, L, escaped)

    return escaped

def save_lattice_snapshot(sites, iteration, p, lambda_rate, L, escaped):
    # Create figure with better styling
    plt.figure(figsize=(12, 12))
    plt.style.use('default')  # Clean style
    
    # Define visualization parameters
    empty_size = 8      # Small size for empty sites
    active_size = 25    # Larger size for active sites
    
    # Define colors - using colorblind-friendly palette
    empty_color = '#f0f0f0'    # Light gray
    predator_color = '#1f77b4'  # Blue
    prey_color = '#d62728'      # Red
    
    # Get coordinates for each site type
    empty_sites = np.where(sites == 0)
    predator_sites = np.where(sites == 1)
    prey_sites = np.where(sites == 2)
    
    # Plot empty sites (small circles, light gray)
    if len(empty_sites[0]) > 0:
        plt.scatter(empty_sites[1], empty_sites[0], 
                   s=empty_size, c=empty_color, marker='o', 
                   alpha=0.3, edgecolors='none', label='Empty')
    
    # Plot predator sites (blue squares, larger)
    if len(predator_sites[0]) > 0:
        plt.scatter(predator_sites[1], predator_sites[0], 
                   s=active_size, c=predator_color, marker='s', 
                   alpha=0.8, edgecolors='darkblue', linewidth=0.5, label='Predators')
    
    # Plot prey sites (red circles, larger)
    if len(prey_sites[0]) > 0:
        plt.scatter(prey_sites[1], prey_sites[0], 
                   s=active_size, c=prey_color, marker='o', 
                   alpha=0.8, edgecolors='darkred', linewidth=0.5, label='Prey')
    
    # Set up the plot
    plt.xlim(-0.5, L-0.5)
    plt.ylim(-0.5, L-0.5)
    plt.gca().invert_yaxis()  # Invert y-axis to match array indexing
    plt.gca().set_aspect('equal')
    
    # Add grid for better visualization
    plt.grid(True, alpha=0.2, linestyle='-', linewidth=0.5)
    
    # Labels and title
    plt.xlabel('Column', fontsize=12)
    plt.ylabel('Row', fontsize=12)
    plt.title(f'Chase-Escape Lattice (p={p}, λ={lambda_rate}, L={L})\n'
              f'Iteration: {iteration}, Escaped: {escaped}', fontsize=14, pad=20)
    
    # Add legend
    plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1), fontsize=10)
    
    # Add border around the lattice
    border = plt.matplotlib.patches.Rectangle((-0.5, -0.5), L, L, 
                                            fill=False, edgecolor='black', linewidth=2)
    plt.gca().add_patch(border)
    
    # Highlight escape boundaries (edges where prey can escape)
    escape_color = 'gold'
    escape_width = 3
    
    # Left boundary (x=0)
    plt.axvline(x=-0.5, color=escape_color, linewidth=escape_width, alpha=0.7)
    # Right boundary (x=L-1)  
    plt.axvline(x=L-0.5, color=escape_color, linewidth=escape_width, alpha=0.7)
    # Top boundary (y=0)
    plt.axhline(y=-0.5, color=escape_color, linewidth=escape_width, alpha=0.7)
    # Bottom boundary (y=L-1)
    plt.axhline(y=L-0.5, color=escape_color, linewidth=escape_width, alpha=0.7)
    
    # Add statistics text
    num_empty = np.sum(sites == 0)
    num_predators = np.sum(sites == 1)
    num_prey = np.sum(sites == 2)
    
    stats_text = f'Empty: {num_empty}\nPredators: {num_predators}\nPrey: {num_prey}'
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
             fontsize=10, verticalalignment='top', 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'./snapshots/near_criticality/lattice_snapshot_l_{L}_p_{p}_lambda_{lambda_rate}_iter_{iteration:06d}.png', 
                dpi=150, bbox_inches='tight')
    plt.close()

# Run single realization with image output
p = 2.5
lambda_rate = 1.0
L = 250

print(f"Running single realization with p = {p}, lambda_rate = {lambda_rate}, L = {L}")
escaped = single_realization_simulation(p, lambda_rate, L)
print(f"Result: Escaped = {escaped}")