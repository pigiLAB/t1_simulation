#Molecule definition file
name = 'Acetate'

#Preset variables
muH = 1.4106E-26 # JA/T
muD = muH * 6.536 / 42.577
muC = muH * 10.705 / 42.577
muN = muH * 3.077 / 42.577
muNa = muH * 11.262 / 42.577
muP = muH * 17.235 / 42.577

nucList = ['C','O','O','C','H','H','H']

#Boxsize of the simulation in angstroms
boxsize = 23

#Timestep of the simulation
dt = 1E-12

#Number of cations (assume Na if positive, Cl if negative))
numIon = 1

#Nucleus of interest
nuc = 0

#Magnetic moment of atoms in the molecule. Change this, along with nuc, to
#alter nucleus of interest. Also change this to deuterate the molecule
molGamma = [muC, 0, 0, 0, muH, muH, muH] #1-13C
#molGamma = [muC, 0, 0, 0, muD, muD, muD] #1-13C Deuterated
#molGamma = [0, 0, 0, muC, muH, muH, muH] #3-13C
#molGamma = [0, 0, 0, muC, muD, muD, muD] #3-13C Deuterated
