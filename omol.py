# %%
"""
Run molecular dynamics (MD) simulations using ASE and a FAIRChem machine
learning interatomic potential.

The script loads a FAIRChem prediction unit from the dataset/model specified
on the command line and uses it through the FAIRChemCalculator with ASE.
An input atomic structure is read, its charge and spin are assigned, and
molecular dynamics is performed using the Velocity Verlet integrator.

The MD trajectory is stored in ``md.traj`` during the simulation and is
subsequently converted to the user-specified output format.

Command-line arguments
----------------------
1. timestep : float
    MD integration timestep in femtoseconds.
2. n_steps : int
    Number of molecular dynamics steps to perform.
3. charge : int
    Total charge of the molecular system.
4. spin : int
    Spin of the molecular system.
5. input_file : str
    Path to the input structure file readable by ASE.
6. output_file : str
    Path to the output trajectory file. The format is inferred by ASE
    from the filename extension.
7. interval : int
    Interval, in MD steps, between trajectory writes/logging.
8. dataset : str
    Dataset/model identifier passed to FAIRChem's ``load_predict_unit``.

Notes
-----
Periodic boundary conditions are disabled. The charge and spin information
is stored in ``atoms.info`` for use with the OMol task.
"""

# Import command-line argument handling.
import sys

# ASE utilities for physical units and atomic structure/trajectory I/O.
from ase import units
from ase.io import Trajectory
from ase.io import read, write

# Velocity Verlet molecular dynamics integrator.
from ase.md.verlet import VelocityVerlet

# FAIRChem calculator providing energies and forces to ASE.
from fairchem.core import FAIRChemCalculator
from fairchem.core.units.mlip_unit import load_predict_unit


# %%
# Read molecular dynamics parameters and file/model information from the
# command-line arguments.

# MD timestep in femtoseconds.
timestep = float(sys.argv[1])

# Total number of MD integration steps.
n_steps = int(sys.argv[2])

# Total molecular charge and spin used by the OMol calculator.
charge = int(sys.argv[3])
spin = int(sys.argv[4])

# Input structure and output trajectory filenames.
input_file = str(sys.argv[5])
output_file = str(sys.argv[6])

# Number of MD steps between trajectory/logging operations.
interval = int(sys.argv[7])

# Dataset/model identifier used to load the FAIRChem prediction unit.
dataset = str(sys.argv[8])


# %%
# Import FAIRChem inference configuration.
from fairchem.core.units.mlip_unit.api.inference import InferenceSettings

# Load the prediction unit from the requested dataset/model.
predictor = load_predict_unit(dataset)


# %%
# Wrap the FAIRChem predictor in an ASE-compatible calculator.
# The "omol" task is used for molecular calculations.
calc = FAIRChemCalculator(predictor, task_name="omol")

# Read the initial atomic structure using ASE.
atoms = read(input_file)

# Store the total charge and spin in the Atoms metadata for the OMol task.
atoms.info.update(charge=charge, spin=spin)  # For omol

# Attach the FAIRChem calculator so ASE can evaluate energies and forces.
atoms.calc = calc


# Disable periodic boundary conditions, treating the system as non-periodic.
atoms.set_pbc(False)


# Create the Velocity Verlet molecular dynamics integrator.
# The user-provided timestep is converted from femtoseconds to ASE units.
# ASE is also instructed to use "md.traj" as its trajectory file.
dyn = VelocityVerlet(atoms=atoms,
       timestep=timestep * units.fs,
       trajectory="md.traj",
       loginterval=interval)

# Open a trajectory file for explicitly writing MD configurations.
trajectory = Trajectory("md.traj", "w", atoms)

# Write the current atomic configuration every `interval` MD steps.
dyn.attach(trajectory.write, interval=interval)

# Run the molecular dynamics simulation for the requested number of steps.
dyn.run(steps=n_steps)

# Close the trajectory writer after the simulation finishes.
trajectory.close()


# %%
# Reopen the generated ASE trajectory for reading.
trajectory_images = Trajectory("md.traj")

# Convert/write the complete trajectory to the requested output file.
# ASE determines the output format from `output_file`.
write(output_file, trajectory_images)

# Close the trajectory reader.
trajectory_images.close()
