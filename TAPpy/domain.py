import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import sys
from . import input
from . import MSD

class Domain:
    def __init__(self,molDict=None):
        # Objects to hold per atom, per mol, domain size etc.
        self.molDict = molDict
        self.MSDDict = None
        self.MSDHDict = None
        self.diffusivityDict = None

        self.boundaries = np.array([[None,None],
                                    [None,None],
                                    [None,None]])
        
        self.simulationUnits = None
        self.timestepLength = None
        self.timesteps = None

        self.atomIDs = None
        self.atomTypes = None
        self.atomCharges = None
        self.atomMolIDs = None

        self.atomTypeMasses = None
        self.atomTypeEpsilons = None
        self.atomTypeSigmas = None

        self.molIDs = None
        self.molAtomIDs = None

        self.molTimeSeries = None
        self.comTimeSeries = None

        self.systemTimeSeries = None
        

    def readInputs(self,dataFile,logFile,trajectoryFile,trim=False):
        input.readLAMMPSdata(self,dataFile)
        input.readLAMMPSlog(self,logFile)
        input.readLAMMPStrajectory(self,trajectoryFile,trim)
        if not self.molDict:
            print('No molecular information provided, defaulting to identical analysis for all molecules.')
            self.molDict = {'mol1':len(self.molIDs)}
        else:
            stringInformation = ['{} {}'.format(self.molDict[key],key) for key in self.molDict]
            stringInformation = ", ".join(stringInformation)
            print('Performing analysis for: {}.'.format(stringInformation))

    def calcMSDWindow(self,fractionalWindow=1):
        print('Calculating MSD.')
        self.MSDDict = {}

        for i in range(int(len(self.timesteps))-int(len(self.timesteps)*fractionalWindow)):
            sys.stdout.write('\rMSD calculation {:.2f}% complete'.format((i+1)*100.0/(int(len(self.timesteps))-int(len(self.timesteps)*fractionalWindow))))
            j = i + int(len(self.timesteps)*fractionalWindow)
            molSlicedArray = self.molTimeSeries[i:j].copy()
            comSlicedArray = self.comTimeSeries[i:j].copy()
            molSlicedArray = molSlicedArray-molSlicedArray[0]
            comSlicedArray = comSlicedArray-comSlicedArray[0]

            molDisplacementArray = molSlicedArray-comSlicedArray
            molDisplacementArray = molDisplacementArray ** 2
            molDisplacementArray = molDisplacementArray.sum(axis=2)

            prev_nMols = 0
            for key in self.molDict:
                nMols = self.molDict[key]+prev_nMols
                MSDArray = molDisplacementArray[:,prev_nMols:nMols]
                MSDArray = MSDArray.mean(axis=1)
                if i == 0:
                    self.MSDDict[key]=MSDArray
                else:
                    self.MSDDict[key]+=MSDArray
                prev_nMols = nMols
        for key in self.molDict:
            self.MSDDict[key] = self.MSDDict[key]/(int(len(self.timesteps))-int(len(self.timesteps)*fractionalWindow)+1)
        sys.stdout.write('\n')

    def cy_calcMSDWindow(self,fractionalWindow=1):
        print('Calculating MSD.')
        self.MSDDict = {}
        for key in self.molDict:
            self.MSDDict[key] = np.asarray(MSD.calculate(self.molTimeSeries,self.comTimeSeries,fractionalWindow))

    def cy_calcMSDHWindow(self,fractionalWindow=1):
        print('Calculating MSD and diffusive energies of activation.')
        self.MSDDict = {}
        self.MSDHDict = {}
        for key in self.molDict:
            combinedMSD = np.asarray(MSD.calculateMSDH(self.molTimeSeries,self.comTimeSeries,self.systemTimeSeries,fractionalWindow))
            self.MSDDict[key] = combinedMSD[:,0]
            self.MSDHDict[key] = combinedMSD[:,1:]

    def calcDiffusivity(self):
        self.diffusivityDict = {}
        for key in self.MSDDict:
            regression = LinearRegression(fit_intercept=False).fit(x,y)
            self.diffusivityDict[key]=regression.coef_[0][0]*1e-8*1e-8/self.timestepLength/1e-15/6
