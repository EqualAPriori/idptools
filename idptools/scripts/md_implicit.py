# simple MD script
# this one both solvates AND runs the explicit solvent simulation

from __future__ import division, print_function

import sys

# OpenMM Imports
import openmm as mm
import openmm.app as app

# ParmEd Imports
from parmed import load_file, unit as u
#from parmed.openmm import StateDataReporter, NetCDFReporter
import mdtraj

# Other imports
import argparse
parser = argparse.ArgumentParser("run implicit solvent")
parser.add_argument("prefix",type=str,help="prefix for .parm7 and .rst7 files")
parser.add_argument("-ps",default=1000,type=int,help="number of picoseconds to run for")
parser.add_argument("-solv",default=8,type=int,choices=[0,1,2,5,7,8],help="igb version")
parser.add_argument("-T",default=298.15,type=float,help="Temp (K)")
parser.add_argument("-salt",default=0.1,type=float,help="salt conc (M)")
parser.add_argument("-stride",default=100,type=int,help="stride for trajectory output (ps)")
parser.add_argument("-init",default=None,help="pdb file for initial coordinates")
parser.add_argument("-device",default=0,type=int,help="device")
args = parser.parse_args()

impSolvDict = {0:None, 1:app.HCT, 2:app.OBC1, 5:app.OBC2, 7:app.GBn, 8:app.GBn2}
impSolv = impSolvDict[args.solv]
dt = 2.0 #fs
steps_per_ps = int(1000/dt)
stride = args.stride * steps_per_ps

# Load the Amber files
print('Loading AMBER files...')
amber_def = load_file(f'{args.prefix}.parm7', f'{args.prefix}.rst7')


# Create the OpenMM system
print('Creating OpenMM System')
system = amber_def.createSystem(nonbondedMethod=app.NoCutoff,
                               constraints=app.HBonds, implicitSolvent=impSolv,
                               implicitSolventSaltConc=args.salt*u.moles/u.liter,
)

# Create the integrator to do Langevin dynamics
integrator = mm.LangevinIntegrator(
                        args.T*u.kelvin,       # Temperature of heat bath
                        1.0/u.picoseconds,  # Friction coefficient
                        dt*u.femtoseconds, # Time step
)

# Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
# the platform to use the default (fastest) platform
platform = mm.Platform.getPlatformByName('CUDA')
prop = dict(CudaPrecision='mixed',CudaDeviceIndex=str(args.device)) # Use mixed single/double precision

# Create the Simulation object
sim = app.Simulation(amber_def.topology, system, integrator, platform, prop)

# Set the particle positions
if args.init is None:
    sim.context.setPositions(amber_def.positions) #CHANGE TO USE PDB AND MDTRAJ
else:
    #struc = mdtraj.load(args.init) #for now, should be a pdb file
    pdb = app.PDBFile(args.init) 
    sim.context.setPositions(pdb.positions)

# Minimize the energy
print('Minimizing energy')
sim.minimizeEnergy(maxIterations=500)


# ===== Reporters and Run =====
totalSteps = steps_per_ps * args.ps

# Set up the reporters to report energies and coordinates every 100 steps
sim.reporters.append(
        app.StateDataReporter(sys.stdout, int(args.stride/10)*steps_per_ps, step=True, potentialEnergy=True,
                                kineticEnergy=True, temperature=True, speed=True, time=True,
                                progress=True, remainingTime=True, totalSteps=totalSteps)
)
sim.reporters.append(
        #NetCDFReporter(f'{args.prefix}.nc', 1000, crds=True) # CHANGE TO USE NETCDF
        mdtraj.reporters.DCDReporter(f'{args.prefix}.dcd', stride)
)

# Run dynamics
positions = sim.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(sim.topology, positions, open("init.pdb", 'w'))

print('Running dynamics')
sim.step(totalSteps)

positions = sim.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(sim.topology, positions, open("final.pdb", 'w'))
