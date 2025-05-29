import numpy as np
import matplotlib.pyplot as plt

# The function below takes the spreading rate of prey (p), the spontaneous conversion rate (lambda_rate),
# and the lattice size L for an L x L square lattice as input and simulates a single realization of the 
# Chase-Escape process on the lattice. As an output, it returns the variable "escaped"
# which is True if the prey manage to escape and False if all the prey die before escaping.

def single_realization_simulation(p, lambda_rate, L):
    # Define constants for site colors
    # 0 = white (Empty site)
    # 1 = blue (Predators)
    # 2 = red (Prey)
    
    # Initialize the lattice with all empty sites
    sites = np.zeros((L, L), dtype=int)
    
    # Set initial conditions (first column blue, second column red)
    sites[:, 0] = 1
    sites[:, 1] = 2

    # Initialize lists to store BR and RE bonds
    BR_bonds = [] # List to store BR bonds
    RE_bonds = [] # List to store RE bonds
    red_sites = [] # List to store all red site coordinates
        
    # Create initial BR and RE bonds along the lattice
    # Also populate the initial red sites list
    for i in range(L):
        BR_bond = [(i, 0), (i, 1)] # BR bond connects first and second column
        RE_bond = [(i, 1), (i, 2)] # RE bond connects second and third column
        BR_bonds.append(BR_bond)
        RE_bonds.append(RE_bond)
        red_sites.append((i, 1)) # All sites in column 1 are initially red

    # Main simulation loop. The program exits this loop when all the prey die
    # or when any of the prey particles reaches the rightmost column of the lattice (escape!)
    
    escaped = 'unknown'
    while escaped == 'unknown':
        rnd = np.random.uniform()
        num_BR_bonds = len(BR_bonds)
        num_RE_bonds = len(RE_bonds)
        num_red_sites = len(red_sites)
        
        # Check termination condition - no red sites left
        if num_red_sites == 0:
            escaped = False
            break
            
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
            if j == L - 1:
                escaped = True
                break

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

    return escaped

# The function below takes as input p, lambda_rate, L, and the number of times you
# want to run the Chase-Escape process. It returns as output the escape
# probability, which denotes the fraction of runs in which the prey have
# managed to escape. 

def simulation(p, lambda_rate, L, num_realizations):
    num_escaped = 0
    for i in range(num_realizations):
        escaped = single_realization_simulation(p, lambda_rate, L)
        if escaped:
            num_escaped += 1

    escape_probability = num_escaped / num_realizations
    return escape_probability

# Example usage with spontaneous conversion
p = 0.49
lambda_rate = 0.01  # Spontaneous red-to-blue conversion rate
L = 100
num_realizations = 100
print(f"p = {p}, lambda_rate = {lambda_rate}, L = {L}, num_realizations = {num_realizations}")

escape_probability = simulation(p, lambda_rate, L, num_realizations)

print(f"Escape probability = {escape_probability}")

# Compare with original model (no spontaneous conversion)
lambda_rate_zero = 0.0
escape_probability_original = simulation(p, lambda_rate_zero, L, num_realizations)
print(f"\nComparison:")
print(f"Escape probability (with spontaneous conversion, λ={lambda_rate}): {escape_probability}")
print(f"Escape probability (original model, λ=0): {escape_probability_original}")