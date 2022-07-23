# simple MD script
# this works for both implicit AND explicit solvent simulation
# TODO:
#   currently hardcodes amber99sb/tip3p .xml files!
#       need to add in settings.json parsing, i.e. to read out what .xml forcefield to use
#       potentially, write a settings parser, or a wrap this script in a function...
#   to reduce redundancy in feeding in arguments...
# 

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
parser = argparse.ArgumentParser("run simulation with openMM library .xml files")
parser.add_argument("prefix",type=str,help="prefix for .parm7 and .rst7 files")
parser.add_argument("-ps",default=1000,type=int,help="number of picoseconds to run for")
parser.add_argument("-stride",default=100,type=int,help="stride for trajectory output (ps)")
parser.add_argument("-solv",default=None,type=int,choices=[None,0,1,2,5,7,8],help="igb version. 0 = vacuum, None = use forcefield as is (e.g. explicit solvent)")
parser.add_argument("-T",default=298.15,type=float,help="Temperature (K)")
parser.add_argument("-p",default=0.0,type=float,help="Pressuer (bar)")
parser.add_argument("-salt",default=0.1,type=float,help="salt conc (M)")
parser.add_argument("-cutoff",default=1.0,type=float,help="cutoff (nm)")
parser.add_argument("-init",default=None,help="pdb file for initial coordinates")
parser.add_argument("-device",default=0,type=int,help="device")
parser.add_argument("-dt",default=2.0,type=float,help="device")
parser.add_argument("-dispcorr",action="store_true",help="turn on dispersion correction")
args = parser.parse_args()

impSolvDict = {None:None, 0:None, 1:app.HCT, 2:app.OBC1, 5:app.OBC2, 7:app.GBn, 8:app.GBn2}
impSolv = impSolvDict[args.solv]
dt = 2.0 #fs
steps_per_ps = int(1000/dt)
stride = args.stride * steps_per_ps

#===== Create System =====
pdb = app.PDBFile(args.init) 
forcefield = app.ForceField('amber99sb.xml', 'amber14/tip3p.xml')
if args.solv is None:
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
        nonbondedCutoff=args.cutoff*u.nanometer, constraints=app.HBonds)
else:
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff,
        nonbondedCutoff=args.cutoff*nanometer, constraints=app.HBonds,
        implicitSolvent=impSolv, implicitSolventSaltConc=args.salt*u.moles/u.liter)


"""
# Load the Amber files
print('Loading AMBER files...')
amber_def = load_file(f'{args.prefix}.parm7', f'{args.prefix}.rst7')


# Create the OpenMM system
print('Creating OpenMM System')
if args.solv is None:
    print("... no implicit solvent method specified, just using input forcefield")
    print("... using cutoff of {}nm".format(args.cutoff))
    system = amber_def.createSystem(nonbondedMethod=app.PME,
                                nonbondedCutoff = args.cutoff*u.nanometer,
                                constraints=app.HBonds,
    )
else:
    system = amber_def.createSystem(nonbondedMethod=app.NoCutoff,
                               constraints=app.HBonds, implicitSolvent=impSolv,
                               implicitSolventSaltConc=args.salt*u.moles/u.liter,
    )
"""

# Add Barostat if applicable
if args.p > 0.0:
    print("... adding barostat ({} bar)".format(args.p))
    if args.solv is not None:
        raise ValueError("... barostat can not be used with implicit solvent")
    barostat_freq = 25
    #barostat = mm.MonteCarloBarostat(args.p * u.bar, args.T * u.kelvin, barostat_freq)
    barostat = mm.MonteCarloAnisotropicBarostat(mm.vec3.Vec3(args.p*u.bar, args.p*u.bar, args.p*u.bar), args.T*u.kelvin, True, True, True, barostat_freq)
    system.addForce(barostat)

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

# Dispersion correction
for i in system.getForces():
    if isinstance(i, mm.NonbondedForce):
        print("... Default dispersionCorrection: {}".format(i.getUseDispersionCorrection())) #Note OpenMM Default is Yes Dispersion Correction!
        i.setUseDispersionCorrection(args.dispcorr)
        print("... New dispersionCorrection: {}".format(i.getUseDispersionCorrection()))

# Create the Simulation object
sim = app.Simulation(pdb.topology, system, integrator, platform, prop)
#sim = app.Simulation(amber_def.topology, system, integrator, platform, prop)

# Set the particle positions
if args.init is None:
    sim.context.setPositions(amber_def.positions) #CHANGE TO USE PDB AND MDTRAJ
else:
    #struc = mdtraj.load(args.init) #for now, should be a pdb file
    pdb = app.PDBFile(args.init) 
    sim.context.setPositions(pdb.positions)

# Minimize the energy
print('Minimizing energy')
print("...Minimization start, the energy is: {}".format(sim.context.getState(getEnergy=True).getPotentialEnergy()))
#sim.minimizeEnergy(maxIterations=500)
sim.minimizeEnergy()
print("...Minimization done, the energy is {}".format(sim.context.getState(getEnergy=True).getPotentialEnergy()))
positions = sim.context.getState(getPositions=True).getPositions()
print("...Minimized geometry is written to 'minimized.pdb'")
app.PDBFile.writeModel(sim.topology, positions, open('minimized.pdb','w'))


# ===== Reporters and Run =====
totalSteps = steps_per_ps * args.ps

# Set up the reporters to report energies and coordinates every 100 steps
logfile = "thermo.log"
if args.solv is not None: #implicit solvent
    sim.reporters.append(
            app.StateDataReporter(logfile, int(args.stride/10)*steps_per_ps, step=True, potentialEnergy=True,
                                   kineticEnergy=True, temperature=True, volume=True,
                                   speed=True, time=True)
    )
    sim.reporters.append(
            app.StateDataReporter(sys.stdout, int(args.stride/10)*steps_per_ps, step=True, potentialEnergy=True,
                                   kineticEnergy=True, temperature=True, volume=True,
                                   speed=True, time=True,
                                   remainingTime=True, progress=True, totalSteps=totalSteps)
    )
else: #probably explicit solvent
    sim.reporters.append(
            app.StateDataReporter(logfile, int(args.stride/10)*steps_per_ps, step=True, potentialEnergy=True,
                                   kineticEnergy=True, temperature=True, density=True, volume=True,
                                   speed=True, time=True)
    )
    sim.reporters.append(
            app.StateDataReporter(sys.stdout, int(args.stride/10)*steps_per_ps, step=True, potentialEnergy=True,
                                   kineticEnergy=True, temperature=True, density=True, volume=True,
                                   speed=True, time=True,
                                   remainingTime=True, progress=True, totalSteps=totalSteps)
    )

sim.reporters.append(
        #NetCDFReporter(f'{args.prefix}.nc', 1000, crds=True) # CHANGE TO USE NETCDF
        mdtraj.reporters.DCDReporter(f'{args.prefix}.dcd', stride)
)

# Run dynamics and save
positions = sim.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(sim.topology, positions, open("init.pdb", 'w'))

print('Running dynamics')
sim.step(totalSteps)

positions = sim.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(sim.topology, positions, open("final.pdb", 'w'))
