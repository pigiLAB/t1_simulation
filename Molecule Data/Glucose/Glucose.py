#Molecule definition file
name = 'Glucose'

#Preset variables
muH = 1.4106E-26 # JA/T
muD = muH * 6.536 / 42.577
muC = muH * 10.705 / 42.577
muN = muH * 3.077 / 42.577
muNa = muH * 11.262 / 42.577
muP = muH * 17.235 / 42.577

#Carbon numbering order in nucList: 4, 3, 2, 1, 5, 6
nucList = ['O','H','C','H','C','H','O','H','C',
	'H','O','H','C','H','O','H','O','C','H','C','H','H','O','H']

#Boxsize of the simulation in nm
boxsize = 23

#Timestep of the simulation
dt = 1E-12

#Number of cations (assume Na)
numIon = 0

#Nucleus of interest
nuc = 2

#Magnetic moment of atoms in the molecule. Change this, along with nuc, to
#alter nucleus of interest. Also change this to deuterate the molecule
molGamma = [0,muH,muC,muH,0,muH,0,muH,0,
	muH,0,muH,0,muH,0,muH,0,0,muH,0,muH,muH,0,muH]


