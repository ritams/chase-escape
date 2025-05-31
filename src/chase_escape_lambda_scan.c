#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

// Structure to represent a bond between two lattice sites
typedef struct {
    int i1, j1;  // First site coordinates
    int i2, j2;  // Second site coordinates
} Bond;

// Structure to represent a coordinate pair
typedef struct {
    int i, j;
} Coordinate;

// Dynamic array to store bonds
typedef struct {
    Bond* bonds;
    int count;
    int capacity;
} BondList;

// Dynamic array to store red site coordinates
typedef struct {
    Coordinate* coords;
    int count;
    int capacity;
} RedSiteList;

// Initialize a bond list
void init_bond_list(BondList* list, int initial_capacity) {
    list->bonds = (Bond*)malloc(initial_capacity * sizeof(Bond));
    list->count = 0;
    list->capacity = initial_capacity;
}

// Initialize a red site list
void init_red_site_list(RedSiteList* list, int initial_capacity) {
    list->coords = (Coordinate*)malloc(initial_capacity * sizeof(Coordinate));
    list->count = 0;
    list->capacity = initial_capacity;
}

// Add a bond to the list
void add_bond(BondList* list, int i1, int j1, int i2, int j2) {
    if (list->count >= list->capacity) {
        list->capacity *= 2;
        list->bonds = (Bond*)realloc(list->bonds, list->capacity * sizeof(Bond));
    }
    list->bonds[list->count].i1 = i1;
    list->bonds[list->count].j1 = j1;
    list->bonds[list->count].i2 = i2;
    list->bonds[list->count].j2 = j2;
    list->count++;
}

// Add a red site to the list
void add_red_site(RedSiteList* list, int i, int j) {
    if (list->count >= list->capacity) {
        list->capacity *= 2;
        list->coords = (Coordinate*)realloc(list->coords, list->capacity * sizeof(Coordinate));
    }
    list->coords[list->count].i = i;
    list->coords[list->count].j = j;
    list->count++;
}

// Remove a bond at specific index
void remove_bond_at_index(BondList* list, int index) {
    if (index < 0 || index >= list->count) return;
    
    // Move last element to the removed position
    if (index != list->count - 1) {
        list->bonds[index] = list->bonds[list->count - 1];
    }
    list->count--;
}

// Remove a red site at specific index
void remove_red_site_at_index(RedSiteList* list, int index) {
    if (index < 0 || index >= list->count) return;
    
    // Move last element to the removed position
    if (index != list->count - 1) {
        list->coords[index] = list->coords[list->count - 1];
    }
    list->count--;
}

// Remove a specific bond from the list
void remove_bond(BondList* list, int i1, int j1, int i2, int j2) {
    for (int i = 0; i < list->count; i++) {
        if (list->bonds[i].i1 == i1 && list->bonds[i].j1 == j1 &&
            list->bonds[i].i2 == i2 && list->bonds[i].j2 == j2) {
            remove_bond_at_index(list, i);
            return;
        }
    }
}

// Remove a specific red site from the list
void remove_red_site(RedSiteList* list, int i, int j) {
    for (int idx = 0; idx < list->count; idx++) {
        if (list->coords[idx].i == i && list->coords[idx].j == j) {
            remove_red_site_at_index(list, idx);
            return;
        }
    }
}

// Free bond list memory
void free_bond_list(BondList* list) {
    free(list->bonds);
    list->bonds = NULL;
    list->count = 0;
    list->capacity = 0;
}

// Free red site list memory
void free_red_site_list(RedSiteList* list) {
    free(list->coords);
    list->coords = NULL;
    list->count = 0;
    list->capacity = 0;
}

// Generate random double between 0 and 1
double uniform_random() {
    return (double)rand() / RAND_MAX;
}

// Generate random integer between 0 and max_val-1
int random_int(int max_val) {
    return rand() % max_val;
}

// Single realization simulation with spontaneous conversion
bool single_realization_simulation(double p, double lambda_rate, int L) {
    // Define constants for site colors
    // 0 = white (Empty site)
    // 1 = blue (Predators)
    // 2 = red (Prey)
    
    // Initialize the lattice with all empty sites
    int** sites = (int**)malloc(L * sizeof(int*));
    for (int i = 0; i < L; i++) {
        sites[i] = (int*)calloc(L, sizeof(int));
    }
    
    // Set initial conditions (first column blue, second column red)
    for (int i = 0; i < L; i++) {
        sites[i][0] = 1;  // First column: predators
        sites[i][1] = 2;  // Second column: prey
    }
    
    // Initialize bond lists and red site list
    BondList BR_bonds, RE_bonds;
    RedSiteList red_sites;
    init_bond_list(&BR_bonds, L * 4);
    init_bond_list(&RE_bonds, L * 4);
    init_red_site_list(&red_sites, L * 4);
    
    // Create initial BR and RE bonds along the lattice
    // Also populate the initial red sites list
    for (int i = 0; i < L; i++) {
        // BR bond connects first and second column
        add_bond(&BR_bonds, i, 0, i, 1);
        // RE bond connects second and third column (if L > 2)
        if (L > 2) {
            add_bond(&RE_bonds, i, 1, i, 2);
        }
        // All sites in column 1 are initially red
        add_red_site(&red_sites, i, 1);
    }
    
    // Main simulation loop
    bool escaped = false;
    
    while (true) {
        double rnd = uniform_random();
        int num_BR_bonds = BR_bonds.count;
        int num_RE_bonds = RE_bonds.count;
        int num_red_sites = red_sites.count;
        
        // Check termination condition - no red sites left
        if (num_red_sites == 0) {
            escaped = false;
            break;
        }
        
        // Calculate total rate and probabilities
        double total_rate = num_BR_bonds + p * num_RE_bonds + lambda_rate * num_red_sites;
        double prob_BR = (double)num_BR_bonds / total_rate;
        double prob_RE = (double)(p * num_RE_bonds) / total_rate;
        
        if (rnd < prob_BR) {
            // Choose a BR bond
            int index = random_int(num_BR_bonds);
            Bond bond = BR_bonds.bonds[index];
            remove_bond_at_index(&BR_bonds, index);
            
            int i = bond.i2;
            int j = bond.j2;
            sites[i][j] = 1;  // The red site in the BR bond is now blue
            
            // Remove this site from red_sites list
            remove_red_site(&red_sites, i, j);
            
            // Identify the neighbors of the newly turned blue site (i,j)
            int up_i = (i - 1 + L) % L;
            int up_j = j;
            int down_i = (i + 1) % L;
            int down_j = j;
            int left_i = i;
            int left_j = j - 1;
            int right_i = i;
            int right_j = j + 1;
            
            // Array of neighbors
            int neighbors[4][2] = {{up_i, up_j}, {down_i, down_j}, {left_i, left_j}, {right_i, right_j}};
            
            // Loop over all neighbors to update the list of active bonds
            for (int n = 0; n < 4; n++) {
                int nbr_i = neighbors[n][0];
                int nbr_j = neighbors[n][1];
                
                // Check bounds (only for left/right neighbors)
                if (nbr_j < 0 || nbr_j >= L) continue;
                
                int nbr_site = sites[nbr_i][nbr_j];
                
                // If the neighbor is red, we need to add a new BR bond
                if (nbr_site == 2) {
                    add_bond(&BR_bonds, i, j, nbr_i, nbr_j);
                }
                // If the neighbor is blue, remove existing BR bond
                else if (nbr_site == 1) {
                    remove_bond(&BR_bonds, nbr_i, nbr_j, i, j);
                }
                // If the neighbor is empty, remove existing RE bond
                else if (nbr_site == 0) {
                    remove_bond(&RE_bonds, i, j, nbr_i, nbr_j);
                }
            }
        }
        else if (rnd < prob_BR + prob_RE) {
            // Choose a RE bond
            int index = random_int(num_RE_bonds);
            Bond bond = RE_bonds.bonds[index];
            remove_bond_at_index(&RE_bonds, index);
            
            int i = bond.i2;
            int j = bond.j2;
            sites[i][j] = 2;  // The empty site in the RE bond is now red
            
            // Add this new red site to red_sites list
            add_red_site(&red_sites, i, j);
            
            // Check escape condition
            if (j == L - 1) {
                escaped = true;
                break;
            }
            
            // Identify the neighbors of the newly turned red site (i,j)
            int up_i = (i - 1 + L) % L;
            int up_j = j;
            int down_i = (i + 1) % L;
            int down_j = j;
            int left_i = i;
            int left_j = j - 1;
            int right_i = i;
            int right_j = j + 1;
            
            // Array of neighbors
            int neighbors[4][2] = {{up_i, up_j}, {down_i, down_j}, {left_i, left_j}, {right_i, right_j}};
            
            // Loop over all neighbors to update the list of active bonds
            for (int n = 0; n < 4; n++) {
                int nbr_i = neighbors[n][0];
                int nbr_j = neighbors[n][1];
                
                // Check bounds (only for left/right neighbors)
                if (nbr_j < 0 || nbr_j >= L) continue;
                
                int nbr_site = sites[nbr_i][nbr_j];
                
                // If the neighbor is red, remove existing RE bond
                if (nbr_site == 2) {
                    remove_bond(&RE_bonds, nbr_i, nbr_j, i, j);
                }
                // If the neighbor is blue, add new BR bond
                else if (nbr_site == 1) {
                    add_bond(&BR_bonds, nbr_i, nbr_j, i, j);
                }
                // If the neighbor is empty, add new RE bond
                else if (nbr_site == 0) {
                    add_bond(&RE_bonds, i, j, nbr_i, nbr_j);
                }
            }
        }
        else {
            // Spontaneous conversion: choose a random red site to convert to blue
            int index = random_int(num_red_sites);
            Coordinate red_site = red_sites.coords[index];
            remove_red_site_at_index(&red_sites, index);
            
            int i = red_site.i;
            int j = red_site.j;
            sites[i][j] = 1;  // Convert red site to blue
            
            // Identify the neighbors of the newly turned blue site (i,j)
            int up_i = (i - 1 + L) % L;
            int up_j = j;
            int down_i = (i + 1) % L;
            int down_j = j;
            int left_i = i;
            int left_j = j - 1;
            int right_i = i;
            int right_j = j + 1;
            
            // Array of neighbors
            int neighbors[4][2] = {{up_i, up_j}, {down_i, down_j}, {left_i, left_j}, {right_i, right_j}};
            
            // Loop over all neighbors to update the list of active bonds
            for (int n = 0; n < 4; n++) {
                int nbr_i = neighbors[n][0];
                int nbr_j = neighbors[n][1];
                
                // Check bounds (only for left/right neighbors)
                if (nbr_j < 0 || nbr_j >= L) continue;
                
                int nbr_site = sites[nbr_i][nbr_j];
                
                // If the neighbor is red, we need to add a new BR bond
                if (nbr_site == 2) {
                    add_bond(&BR_bonds, i, j, nbr_i, nbr_j);
                }
                // If the neighbor is blue, remove existing BR bond
                else if (nbr_site == 1) {
                    remove_bond(&BR_bonds, nbr_i, nbr_j, i, j);
                }
                // If the neighbor is empty, remove existing RE bond
                else if (nbr_site == 0) {
                    remove_bond(&RE_bonds, i, j, nbr_i, nbr_j);
                }
            }
        }
    }
    
    // Clean up memory
    for (int i = 0; i < L; i++) {
        free(sites[i]);
    }
    free(sites);
    free_bond_list(&BR_bonds);
    free_bond_list(&RE_bonds);
    free_red_site_list(&red_sites);
    
    return escaped;
}

// Run multiple simulations and calculate escape probability
double simulation(double p, double lambda_rate, int L, int num_realizations) {
    int num_escaped = 0;
    
    for (int i = 0; i < num_realizations; i++) {
        bool escaped = single_realization_simulation(p, lambda_rate, L);
        if (escaped) {
            num_escaped++;
        }
    }
    
    return (double)num_escaped / num_realizations;
}

int main(int argc, char *argv[]) {
    // Check command line arguments
    if (argc < 4) {
        printf("Usage: %s <L> <start_index> <end_index>\n", argv[0]);
        printf("  L: lattice size\n");
        printf("  start_index: starting realization index (inclusive)\n");
        printf("  end_index: ending realization index (exclusive)\n");
        printf("\nExample for array job parallelization:\n");
        printf("  Job 0: %s 50 0 10\n", argv[0]);
        printf("  Job 1: %s 50 10 20\n", argv[0]);
        printf("  Job 999: %s 50 9990 10000\n", argv[0]);
        printf("\nThis version scans lambda from 0.26 to 0.29 with fixed p = 1.0\n");
        return 1;
    }
    
    // Parse command line arguments
    int L = atoi(argv[1]);
    int start_index = atoi(argv[2]);
    int end_index = atoi(argv[3]);
    
    if (L <= 0) {
        printf("Error: L must be a positive integer\n");
        return 1;
    }
    
    if (start_index < 0 || end_index <= start_index) {
        printf("Error: Invalid index range. start_index must be >= 0 and end_index > start_index\n");
        return 1;
    }
    
    // Calculate number of realizations for this job
    int num_realizations = end_index - start_index;
    
    // Seed random number generator with unique seed for each job
    srand(time(NULL) + start_index);
    
    // Parameters - fixed p value and lambda scan range
    double p = 1.0;  // Fixed p value
    int num_lambda_values = 11;
    double lambda_start = 0.26;
    double lambda_end = 0.29;
    
    // Create filename based on L and indices
    char filename[150];
    sprintf(filename, "escape_prob_lambda_scan_L%d_p%.2f_%d_%d.txt", L, p, start_index, end_index);
    
    // Open file for writing
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        return 1;
    }
    
    // Write header to file
    fprintf(file, "# lambda\tescape_probability\n");
    fprintf(file, "# Fixed p = %.2f, Realizations: %d (indices %d to %d)\n", p, num_realizations, start_index, end_index - 1);
    
    printf("Running lambda scan simulations for L = %d, p = %.2f\n", L, p);
    printf("Lambda range: %.3f to %.3f\n", lambda_start, lambda_end);
    printf("Realizations: %d (indices %d to %d)\n", num_realizations, start_index, end_index - 1);
    printf("Saving results to %s\n", filename);
    printf("Progress: ");
    
    // Run simulations for lambda values from 0.26 to 0.29
    for (int i = 0; i < num_lambda_values; i++) {
        double lambda_rate = lambda_start + (lambda_end - lambda_start) * i / (num_lambda_values - 1.0);
        
        // Run simulation
        double escape_probability = simulation(p, lambda_rate, L, num_realizations);
        
        // Save result to file
        fprintf(file, "%.4f\t%.6f\n", lambda_rate, escape_probability);
        
        // Print progress
        printf("%.1f%% ", 100.0 * (i + 1) / num_lambda_values);
        fflush(stdout);
    }
    
    printf("\nSimulations complete. Results saved to %s\n", filename);
    
    // Close the file
    fclose(file);
    
    return 0;
} 