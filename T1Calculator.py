#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python script to calcualte estiamted T1 values from a molecular
dynamics simulation.

    - Assume trajectory files are numpy arrays of the form:
        - [[[1x 1y 1z],[2x 2y 2z], ...
    - Need to create T1Calculator object using molecule definition file
    
    Functions:
        - writeGauss/buildSubSystem - helper functions to generate
            Gaussian input files for each timestep. These are generally not
            called directly, but instead through calcQM
        - calcQM - in parallel (if desired), run through the entire
            trajectory using Gaussian to calculate the QM values
        - simT1 - once all of the QM data is generated, estimate 1/T1
            for the given trajectory and molecule
        - plotT1 - helper function to plot the estiamted T1 with appropriate
        	labels
        - main - if run from the command line, this script will read in
            a raw trajectory and remove the dummy atoms from TIP4P water.


@author: Joe Wildenberg

Adapted from pseudocode developed by by Steve Kadlecek
"""
import multiprocessing as mult
from pathos.multiprocessing import ProcessingPool as Pool
import os
import sys
import subprocess
import numpy as np
from math import pi
from math import sqrt
import matplotlib.pyplot as plt
import importlib
import MDAnalysis as mda

#T1Calculator class contains all the methods to perform the pipeline
# Usage for piecemeal method usage:
#       - import T1Calculator
#       - import "molecular definition file" as molecule
#       - t1 = T1Calculator(mol)
#     Can change default nucleus of interest and external B field if desired
#
# The whole pipeline can be performed automatically through the __main__ function
class T1Calculator:
    def __init__(self, molecule, B0=9.4):
        #Store all of the module information in this object
        self.name = molecule.name
        self.nucList = molecule.nucList
        self.molGamma = molecule.molGamma
        self.numIon = molecule.numIon
        self.boxsize = molecule.boxsize
        
        #Values that can be modified through object creation
        self.nuc = molecule.nuc
        self.B0 = B0
        
        
        print('T1Calculator created with molecule: ' + self.name)
        
        #----A bunch of constants 
        #Directory to write Gaussian logfiles, used by runGauss
        #Need to modify this to specific machine
        self.gauDir = 'Calculations/'
        
        #Use 0.5 nm from molecule center as cutoff for Gaussian subsystem
        self.gauDist = 5
        
        #Magnetic moments
        self.muH = 1.4106E-26 # JA/T for a proton
        self.muD = self.muH * 6.536 / 42.577
        self.muNa = self.muH * 11.262 / 42.577
        self.muW = self.muH #Set hydrogen isotope for waters
        self.mu0o4pi = 1E-7
        
        #Values for the output of the molecular simulation
        self.dt = molecule.dt #timestep is 1 ps
        self.dx = 1E-10  #units of Gromacs coordinates (Angstroms)

    #Generate the Gaussian input file 
    #	atomLabels - vector of character labels (e.g. 'H', 'O', 'C')
    #	coords - N x 3 matrix in same order as atomLabels
    #	fileName - input filename to use
    def writeGauss(self, atomLabels, coords, fileName):
        #Gaussian uses Angstroms
        coords = coords * (self.dx/1E-10)
        
        #These first lines are options for Gaussian. Edit for your system
        outFile = open(fileName, 'w')
        outFile.write('%mem=1GB\n') #How much memory to give each job
        #outFile.write('%CPU=0-3\n') #If want specific CPUs
        #outFile.write('%GPUCPU=0=0\n') #If want to use GPUs
        outFile.write('%chk=temp'+fileName+'.chk\n')
        outFile.write('# B3LYP/6-31 NMR Pop=Minimal Prop NoSymmetry\n')
        outFile.write('\n')
        outFile.write(self.name + '\n')
        outFile.write('\n')
        outFile.write(str(-self.numIon) + ' 1\n')
        
        for i in range(len(atomLabels)):
            outFile.write(' ' + atomLabels[i] + '\t' +
                          '\t'.join(map(str,coords[i,:])) + '\n')
        outFile.write('\n')
        outFile.write('--Link1--\n')
        outFile.write('%mem=1GB\n')
        outFile.write('%chk=temp'+fileName+'.chk\n')
        outFile.write('# B3LYP/6-31 Geom=AllCheck Pop=MK Density=Checkpoint NoSymmetry\n')
        outFile.write('\n')
        outFile.close()
        
    #Create the subsystem that will be submitted to Gaussian. Calls writeGauss to
    #write the input file
    #  data - data[0] is timepoint (number) for Gaussian file naming
    #		  - data[1] are coordinates in N x 3 matrix
    def buildSubSystem(self, data):
        i = data[0]
        curPoint = data[1]
        numNuc = len(self.molGamma)
        totalNuc = curPoint.shape[0]
        waters = int((totalNuc - numNuc - abs(self.numIon)) / 3)
        
        #Start by centering the molecule in the box and correcting for jumps
        refLoc = np.repeat(curPoint[np.newaxis,self.nuc,:], totalNuc, axis=0)
        cent = curPoint - refLoc
        jumps = np.where(abs(cent) > self.boxsize / 2)
        #Correct for jumps by wrapping
        cent[jumps] = -1*np.sign(cent[jumps]) * (self.boxsize - abs(cent[jumps]))
        
        #This will contain the atoms we want to include in the Gaussian calculation
        atomLocs = []
                
        #Start by adding the molecule itself
        for j in range(numNuc):
            atomLocs.append(cent[j,:])
        
        for w in range(waters):
            oxyLoc = cent[numNuc + 3*w,:] #Location of this oxygen
            if np.linalg.norm(oxyLoc) < self.gauDist: #COI is at origin
                atomLocs.append(oxyLoc)#O
                atomLocs.append(cent[numNuc + 3*w + 1,:])#H1
                atomLocs.append(cent[numNuc + 3*w + 2,:])#H2
        
        atomLabels = self.nucList + (['O', 'H', 'H'] *
                                    int((len(atomLocs) - numNuc) / 3))
        
        self.writeGauss(atomLabels, np.array(atomLocs),
                        self.name + '-' + str(i) + 'NMR.com')
        
        #Run g16
        logFile = self.name + '-' + str(i) + 'NMR.log'
        if not os.path.isfile(logFile): #Check if the calculation has alreay run
            inFile = open(self.name + '-' + str(i) + 'NMR.com', 'r')
            outFile = open(logFile, 'w')
            subprocess.run('g16', stdin=inFile, stdout=outFile)
            inFile.close()
            outFile.close()
        
    #After using a molecular dynamics program to generate the trajectory
    #we need to use Gaussian to perform the QM calculations. This method
    #and its helper methods (buildSubSystem and writeGauss) will generate
    #a Gaussian input file, run the calculation, and extract the output
    #   - data - the npz file with the trajectory in 'traj'
    #   - processes - the number of CPUs to use in parallelization
    def calcQM(self, data, parallel=mult.cpu_count()):
        #Change our working directory to gauDir+name and then back at the end
        curPath = os.getcwd()
        os.chdir(self.gauDir + self.name)
        
        traj = data['traj']
        totT = np.shape(traj)[0]
        numNuc = len(self.molGamma)
        
        #runList keeps track of the timepoints that still need to be calculated
        #only add if the log file does not yet exist or is incomplete
        runList = []
    
        #Create the variables that we will read in from Gaussian
        shield = np.zeros([totT,numNuc,3,3], np.float32)#Shield matrix
        chargesMul = np.zeros([totT,numNuc], np.float32)
        chargesESP = np.zeros([totT,numNuc], np.float32)
        potential = np.zeros([totT,numNuc], np.float32)
        field = np.zeros([totT,numNuc,3], np.float32) #X, Y, Z
        gradient = np.zeros([totT, numNuc, 6], np.float32)# XX, YY, ZZ, XY, XZ, YZ
        
        #Read in any log files that already exist
        for t in range(totT):
            logFile = self.name + '-' + str(t) + 'NMR.log'
            #If the file doesn't exist, add to runlist and move to next timepoint
            if not os.path.isfile(logFile):
                runList.append([t, traj[t,:,:]])
                continue
            
            #Look at log file and read in all variables
            inFile = open(logFile, 'r')
            lines = inFile.readlines()
            inFile.close()
            
            ##NMR properties - Find the start of the NMR calculation
            NMRLoc = -1
            for j, l in enumerate(lines):
                if l.find("SCF GIAO Magnetic shielding tensor (ppm):") > -1:
                    NMRLoc = j
                    break
                
            if NMRLoc > 0:
                NMRdata = lines[NMRLoc+1:]
                for nuc in range(numNuc):
                    #Skip to the position of the nucleus of interest
                    nucData = NMRdata[5*nuc:5*(nuc+1)]
                
                    #Get each 3 separately - nucData[0] is total chemical shift
                    #Need to split on both whitespace and = as sometimes big
                    # values don't have a whitespace.
                    line = nucData[1].replace('=',' ').split() #XX YX ZX
                    shield[t,nuc,0,:] = [float(line[1]), float(line[3]), float(line[5])]
                    line = nucData[2].replace('=',' ').split() #XY YY ZY
                    shield[t,nuc,1,:] = [float(line[1]), float(line[3]), float(line[5])]
                    line = nucData[3].replace('=',' ').split() #XZ YZ ZZ
                    shield[t,nuc,2,:] = [float(line[1]), float(line[3]), float(line[5])]
    
            ##Mulliken charges
            MulLoc = -1
            for j, l in enumerate(lines):
                if l.find("Mulliken charges:") > -1:
                    MulLoc = j
                    break
                
            if MulLoc > 0:
                ChargeData = lines[MulLoc+2:]
                for nuc in range(numNuc):
                    chargesMul[t,nuc] = float(ChargeData[nuc].split()[2])
                    
            ##Electrostatic properties - Find the start of the Prop calculation
            PropLoc = -1
            for j,l in enumerate(lines):
                if l.find("Electrostatic Properties (Atomic Units)") > -1:
                    PropLoc = j
                    break
            
            if PropLoc > 0:
                PotField = lines[PropLoc+6:]
                for nuc in range(numNuc):
                    line = PotField[nuc].split()
                    potential[t,nuc] = float(line[2])
                    field[t,nuc,:] = [float(line[3]), float(line[4]), float(line[5])]
                for j,l in enumerate(PotField):
                    if l.find("Electric Field Gradient") > -1:
                        GradLoc = j
                        break
                Grad = PotField[GradLoc+3:]
                for nuc in range(numNuc):
                    line = Grad[nuc].split()
                    gradient[t,nuc,:3] = [float(line[2]), float(line[3]), float(line[4])]
                for j,l in enumerate(Grad):
                    if l.find("Electric Field Gradient") > -1:
                        GradLoc = j
                        break
                Grad2 = Grad[GradLoc+3:]
                for nuc in range(numNuc):
                    line = Grad2[nuc].split()
                    gradient[t,nuc,3:] = [float(line[2]), float(line[3]), float(line[4])]
    
            ##ESP Charges:
            ESPLoc = -1
            for j,l in enumerate(lines):
                if l.find(" ESP charges:") > -1:
                    ESPLoc = j
                    break
            
            if ESPLoc > 0:
                ChargeData = lines[ESPLoc+2:]
                for nuc in range(numNuc):
                    chargesESP[t,nuc] = float(ChargeData[nuc].split()[2])
                    
            #If the data is not found, rename the log file and recalculate
            else:
                print("Error! NMR/Electrostatic Calculation not found in file " + logFile)
                os.rename(logFile,logFile + '-fail.log')
                runList.append([t, traj[t,:,:]])
        
        #There are datapoints to calculate. Do this in parallel
        if len(runList) > 0:
            #Run the remaining calcualtions in parallel
            print("Calculating data for " + str(len(runList)) + " timepoints")
            p = Pool(parallel)
            p.map(self.buildSubSystem, runList)
            p.clear()
            print("Must run calcShield again to get results")
            os.chdir(curPath)
            return
               
        os.chdir(curPath)
        return {"traj":traj, "shield":shield, "chargesMul":chargesMul, "chargesESP":chargesESP,
                "potential":potential, "field":field, "gradient":gradient}

    
    #This function performs the main calculation of T1 for a single trajectory
    #Inputs:
        # data - list containg 
        #   1) trajectory matrix of shape t x n x 3 where n is number of atoms
        #   2) shield matrix t x 3 x 3 where t is the number of timepoints
        # numSamp - the number of samples to calculate
        # chargesESP - default is to use ESP charges. If not true, use Mullikin
    def simT1(self, data, numSamp=100, fullResults=True, chargesESP=True):
        traj = self.dx * data['traj'].astype(np.float64)
        shield = 1E-6 * data['shield'].astype(np.float64)
        field = data['field'].astype(np.float64) * 5.14E+11 #convert to mks
        
        if chargesESP:
            charges = data['chargesESP'].astype(np.float64)
        else:
            charges = data['chargesMul'].astype(np.float64)
         
        #Pull in variables from mol module
        boxsize = 1E-10 * self.boxsize #to m
        B0 = self.B0
        numNuc = len(self.molGamma)
        totT = traj.shape[0]
        totalNuc = traj.shape[1]
        otherNucsAll = list(range(self.nuc)) + list(range(self.nuc+1,totalNuc))
        waters = int((totalNuc - numNuc - abs(self.numIon)) / 3) #Number of water molecules
        
        #-- Charge information for the system
        waterCharges = [-0.779513, 0.403529, 0.322110] #average water charges
        allCharges = np.zeros([totT, totalNuc], np.float64)
        allCharges[:,:numNuc] = charges #gamma is an array of gyromagnetic ratios
        if self.numIon == 0:
            allCharges[:,numNuc:] = np.tile(waterCharges, [totT, waters])
        else:
            allCharges[:,numNuc:-abs(self.numIon)] = np.tile(waterCharges, [totT, waters])
            allCharges[:,-abs(self.numIon):] = np.tile([int(self.numIon/abs(self.numIon))], [totT,1])
        
        #Convert to mks
        allCharges *= self.mu0o4pi * 1.602E-19
        
        #-- Magnetic moment information for the system
        gamma = np.zeros(totalNuc, np.float64)
        gamma[:numNuc] = self.molGamma #gamma is an array of gyromagnetic ratios
        if self.numIon == 0:
            gamma[numNuc:] = np.repeat([[0, self.muW, self.muW]],
                 waters, axis=0).reshape(-1)
        else:
            gamma[numNuc:-abs(self.numIon)] = np.repeat([[0, self.muW, self.muW]],
                  waters, axis=0).reshape(-1)
            if self.numIon > 0:#Only need to do if Na as Cl gamma=0
                gamma[-abs(self.numIon):] = np.repeat([self.muNa],self.numIon)
        omega = 2.675E+8 * gamma / self.muH #precession rate per Tesla
        
        
        ## Most of the calculations can be done using numpy's array logic
        ## and not inside the loop over the trajectory. Do this to save lots of time
        ##------------------
        
        #--Start by centering the system on the nucleus of interest
        #  and adjusting for "jumps" across the boundry.
        #  We are going to need to calculate velocities where there are no jumps
        #  in between timesteps, so do for current timestep and also 
        #  last timestep but using this timestep's location as reference.
        refLoc = np.repeat(traj[:,self.nuc,np.newaxis,:], totalNuc, axis=1)
        cent = traj - refLoc
        jumps = np.where(abs(cent) > boxsize / 2)
        #Correct for jumps by wrapping
        cent[jumps] = -1*np.sign(cent[jumps]) * (boxsize - abs(cent[jumps]))
       
        #Now do same for previous timestep, centered on current timestep
        centLast = traj[:-1,:,:] - refLoc[1:,:,:]
        jumpsLast = np.where(abs(centLast) > boxsize / 2)
        centLast[jumpsLast] = -1*np.sign(centLast[jumpsLast]) * (boxsize - 
                abs(centLast[jumpsLast]))
        
        deltar = cent[1:,:,:] - centLast
        
        #Last, we have to correct for any jumps between centLast and cent 
        #otherwise the velocities will be way off
        jumpVel = np.where(abs(deltar) > boxsize / 2)
        deltar[jumpVel] = -1*np.sign(deltar[jumpVel]) * (boxsize - abs(deltar[jumpVel]))
        
        #Calculate velocity vector
        vt = deltar / self.dt
        
        #Perform distance calculations and powers (NOI is Nucleus Of Interest)
        rNOI = cent #with NOI at origin, r is just the coordinates
        dNOI = np.linalg.norm(rNOI, axis=2)
        dNOI3 = np.power(dNOI, 3)
        dNOI5 = np.power(dNOI, 5)
                
        
        ##--Spin-rotation. These calculations can be done once and added
        #  to the B-field in the sampling loop
        #  Iterate through each timepoint, though start at 1 as lose
        #  first timepoint to get velocities   
        if fullResults == True or any("SR" in s for s in fullResults):
            #Current contribution - movment of charged nuclei
            #At the element-wise level, this is:
            #    allCharges[t,n]*np.cross(vt, rRel)/magVel3[t,n]
            #       but need to expand allCharges and magVel3 to do element-wise
            #Then sum over all nuclei
            #Also, do not want to include nuclei of interest
            #tempBSR1 = np.sum(np.repeat(allCharges[1:,otherNucsMol,np.newaxis],3, axis=2) *
            #           np.cross(vt, rRel)[:,otherNucsMol,:] / 
            #           np.repeat(magVel3[:,otherNucsMol,np.newaxis],3, axis=2), axis=1)
            tempBSR1 = np.sum(np.repeat(allCharges[1:,otherNucsAll,np.newaxis],3, axis=2) *
                              np.cross(vt,rNOI[1:,:,:])[:,otherNucsAll,:] /
                              np.repeat(dNOI3[1:,otherNucsAll,np.newaxis],3,axis=2), axis=1)
            #Moving E field contribution
            tempBSR2 = np.cross(field[1:,self.nuc,:], 
                                vt[:,self.nuc,:]) / np.power(2.998E+8,2)
        
        #----------------------------Now we can use these variables in each numSamp
        results = np.zeros([numSamp, totT-1],np.float64)
        for samp in range(numSamp):
            #At each timepoint, there is a static B0 field in the z-direction
            B = np.repeat([[0, 0, B0]], totT, axis=0)
            
            #--Dipole-dipole interactions
            if fullResults == True or any("Dipole" in s for s in fullResults):
                #Randomize the vectors. Need to set the generator seed as
                #np.random and multiprocessing do not play well together
                #(e.g. each process gets the same seed by default)
                np.random.seed(int.from_bytes(os.urandom(4), sys.byteorder))
                theta = pi * np.random.random(totalNuc)
                phi = 2 * pi * np.random.random(totalNuc)
            
                #Calculate the mu vector
                mp = self.mu0o4pi * np.multiply(gamma[:,np.newaxis], 
                                   np.array([np.multiply(np.sin(theta), np.cos(phi)),
                                             np.multiply(np.sin(theta), np.sin(phi)),
                                             np.cos(theta)]).transpose())
        
                #Now iterate through each nucleus
                for n in otherNucsAll:
                    if gamma[n] != 0:#Don't need to calculate if no gamma
                        dot = np.inner(mp[n,:],rNOI[:,n,:])[:,np.newaxis]
                        B += (3 * dot * rNOI[:,n,:] / (dNOI5[:,n])[:,np.newaxis] 
                            - mp[n,:] / (dNOI3[:,n])[:,np.newaxis])
                    
            ##--Spin-rotation
            #  Iterate through each timepoint, though start at 1 as lose
            #  first timepoint to get velocities
            if fullResults == True or any("SR" in s for s in fullResults):
                B[1:,:] += tempBSR1 + tempBSR2
                    
            ##--Chemical Shift Aniosotropy
            #  For each timepoint, use the shielding tensor's interaction
            #  with the local field
            if fullResults == True or any("CSA" in s for s in fullResults):
                for t in range(totT):
                    B[t,:] += shield[t,self.nuc,:,:].dot(B[t,:])
        
            #Now evolve mu as d mu / dt =  omega mu x B,
            muclass = np.zeros([totT, 3], np.float64)
            muclass[0,:] = [0, 0, 1]#Initialize the first timepoint
            for t in range(totT-1):
                muclass[t+1,:] =  muclass[t,:] + (omega[self.nuc] * 
                        np.cross(muclass[t,:], B[t,:]) * self.dt)
        
            #relaxation rate is 1 - z projection of mu / t
            #1/T1 = (mu[x]^2 + mu[y]^2) / (2*t*dt) - All timepoints at once
            results[samp,:] = (np.square(muclass[1:,0]) + np.square(muclass[1:,1])) / (
                    (2 * (np.array(range(totT-1))+1) * self.dt))
        
        
        return results.transpose()

    #Calculate the T1 estimation. This is done separately for each simulation
    #in results. After averaging, the returned value is the estimate over
    #the final 80% of the trajectory to allow for equibrilization
    def plotT1(self, results, plot=True, sd=False,
               extrap=False, ymax=-1, label=None, linestyle=None):
           
        x = self.dt * (np.array(range(results.shape[0]))+1) / 1E-9 #in ns
                    
        #Each column is a separate simulation of the same nucleus, plot
        #the average and standard deviation
        mean = 1./np.mean(results, axis=1)
        
        #Calc estimated T1 from last 80%
        t1Val = np.mean(mean[int(1*len(mean)/5):])
        
        ci = 1.96 * (1./np.std(results, axis=1)) / sqrt(results.shape[1])
        if plot:
            plt.plot(x, mean, label=label,linestyle=linestyle)
            if sd:
                ci = 1.96 * (1./np.std(results, axis=1)) / sqrt(results.shape[1])
                plt.fill_between(x, mean-ci, mean+ci, alpha=0.5, label='95% CI')
                
                if extrap:
                    plt.plot((x[int(1*len(x)/5)],np.max(x)), (t1Val,t1Val), 'k--',
                             label=r'$T_1$' + ' Value')
        
            if ymax > 0:
                plt.ylim([0, ymax])
            else:
                plt.autoscale(axis='y')
            plt.xlabel('Simulation Time (ns)',size='x-large')
            plt.ylabel('Estimated ' + r'$T_1$' + ' (s)',size='x-large')
     
        return t1Val, ci



#Function so this script can be run from the command line as part of
#a processing pipeline 
#
# 1) Load the molecule definition from mol.py
# 2) Load molecular simulation trajectory and process
#       - Remove dummy atoms (e.g. TIP4P water dummy)
#       - Save in numpy format
# 3) Option: QM - Calculate QM files and save results
#       - Followed by number of processors. If nothing, assume all available
#
#  Usage: python T1Calculator.py mol.py (QM)
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('No molecule file provided')
    else:
        mol = importlib.import_module(str(sys.argv[1]))
    
    numNuc = len(mol.molGamma)
    
    #Read in the data file generated by the molecular dynamics program
    u = mda.Universe('md.gro', 'md-nojump.trr',in_memory=True)
    atoms = u.select_atoms('not type DUMMY') #Remove dummy atoms from TIP4P water
    
    #The rest of this script has time as first axis, so swap
    #timeseries is in units of Angstroms
    traj = u.trajectory.timeseries(atoms).swapaxes(0,1)
    
    #Save new trajectory
    with mda.Writer(mol.name+'-nodummy.trr', atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(atoms)
    atoms.write(mol.name+'-nodummy.gro')
    
    np.savez_compressed(mol.name + '.npz', traj=traj.astype(np.float16), timestep=mol.dt)
    
    sys.exit(0)
    
