# Preprocessing

This takes ligandBonds.coords and ligandCoordinates.coords and generates P22CD40L.in and lammps_script.in

### ligandCoordinates.coords

This is a user generated file of all xyz coordinates of a ligand complex around P22 centered at (0,0,0)

### ligandBonds.coords

This is a user generated file of bonds between the ligands and vlp

### P22CD40L.in

This is the input file containing all atom and bond information for the lammps script

### lammps_script.in

This is the lammps script that will be run

## Running the postprocessor

* First, git clone the project:
'''git clone https://github.com/softmaterialslab/nanobind.git '''

* Install the project:
'''make'''

* Run the project:
''' ./generate_lammps_script'''

* If you want to clean everything and create a new build, use:
'''make clean'''

