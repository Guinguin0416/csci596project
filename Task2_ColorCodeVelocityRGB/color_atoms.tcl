# Example Tcl script for VMD
# This script changes the color of all atoms based on their x-coordinate

# Select all atoms
set all [atomselect top "all"]

# Get the x-coordinates of all atoms
set x_coords [$all get x]

# Initialize a list to store beta values
set beta_values {}

# Loop through each atom to set its beta value
for {set i 0} {$i < [$all num]} {incr i} {
    # Get the x-coordinate of the ith atom
    set x [lindex $x_coords $i]

    # Determine beta value based on x-coordinate
    if {$x > 3} {
        lappend beta_values 1
    } else {
        lappend beta_values 0
    }
}

# Assign the beta values to the atoms
$all set beta $beta_values

# Apply the changes made to the atom selection
$all update

# Dispose of the atom selection to free memory
$all delete

# Change the representation to color by beta value
mol modcolor 0 top Beta
mol modstyle 0 top VDW 1.0 12.0
