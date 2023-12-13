import numpy as np

# Parameters
num_atoms_per_side = 10
lattice_constant = 3.4  # Angstroms, roughly for Argon

# Derived quantities
num_atoms = num_atoms_per_side ** 3
box_size = num_atoms_per_side * lattice_constant

# Create positions (simple cubic lattice)
positions = np.array([[x, y, z] for x in np.linspace(0, box_size, num_atoms_per_side, endpoint=False)
                                 for y in np.linspace(0, box_size, num_atoms_per_side, endpoint=False)
                                 for z in np.linspace(0, box_size, num_atoms_per_side, endpoint=False)])

# Write LAMMPS data file
with open('lammps_initial_data.txt', 'w') as file:
    file.write('LAMMPS Description\n\n')
    file.write(f'{num_atoms} atoms\n\n')
    file.write('1 atom types\n\n')
    file.write(f'0.0 {box_size:.4f} xlo xhi\n')
    file.write(f'0.0 {box_size:.4f} ylo yhi\n')
    file.write(f'0.0 {box_size:.4f} zlo zhi\n\n')
    file.write('Masses\n\n')
    file.write('1 39.948 # Argon\n\n')
    file.write('Atoms # atomic\n\n')
    
    for i, pos in enumerate(positions, start=1):
        # Assuming atom type 1 for all atoms
        file.write(f'{i} 1 {pos[0]:.4f} {pos[1]:.4f} {pos[2]:.4f}\n')
