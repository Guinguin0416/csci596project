from lammps import lammps

# Create a LAMMPS object
lmp = lammps()

# Path to your LAMMPS input file
input_file = '/Users/yisu/Desktop/csci596project/Task1_ColorKE/a.in'

# Load the LAMMPS input script
lmp.file(input_file)

# The simulation runs according to the commands in your input file
