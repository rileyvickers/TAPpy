import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import sys
from . import input

class Domain:
    def __init__(self,molDict=None):
        # Objects to hold per atom, per mol, domain size etc.
        self.molDict = molDict
        self.MSDDFDict = None
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

        self.molTimeSeriesDF = None
        self.comTimeSeriesDF = None
        

    def readInputs(self,dataFile,logFile,trajectoryFile):
        input.readLAMMPSdata(self,dataFile)
        input.readLAMMPSlog(self,logFile)
        input.readLAMMPStrajectory(self,trajectoryFile)
        if not self.molDict:
            print('No molecular information provided, defaulting to identical analysis for all molecules.')
            self.molDict = {'mol1':len(self.molIDs)}
        else:
            stringInformation = ['{} {}'.format(self.molDict[key],key) for key in self.molDict]
            stringInformation = ", ".join(stringInformation)
            print('Performing analysis for: {}.'.format(stringInformation))
    
    def calcDisplacement(self):
        self.molTimeSeriesDF[['xdiff','ydiff','zdiff']] = self.molTimeSeriesDF[['x','y','z']]-self.molTimeSeriesDF.loc[0][['x','y','z']]
        self.comTimeSeriesDF[['xdiff','ydiff','zdiff']] = self.comTimeSeriesDF[['x','y','z']]-self.comTimeSeriesDF.iloc[0][['x','y','z']]
        self.molTimeSeriesDF = self.molTimeSeriesDF.reset_index(level=1)
        self.molTimeSeriesDF[['xdiff','ydiff','zdiff']] = self.molTimeSeriesDF[['xdiff','ydiff','zdiff']]-self.comTimeSeriesDF[['xdiff','ydiff','zdiff']]
        self.molTimeSeriesDF['SD'] = self.molTimeSeriesDF['xdiff']**2 + self.molTimeSeriesDF['ydiff']**2 + self.molTimeSeriesDF['zdiff']**2
        self.molTimeSeriesDF.set_index('mol',append=True,inplace=True)
        self.molTimeSeriesDF.reset_index(level=0,inplace=True)
        self.molTimeSeriesDF.rename(columns={'level_0':'timestep'},inplace=True)
    
    def calcMSD(self):
        prev_nMols = 0
        self.MSDDFDict = {}
        for key in self.molDict:
            nMols = self.molDict[key]+prev_nMols
            slicedDF = self.molTimeSeriesDF.loc[(self.molTimeSeriesDF.index>prev_nMols)&(self.molTimeSeriesDF.index<=nMols)]
            MSDDF = slicedDF['SD'].groupby(slicedDF['timestep']).mean()
            MSDDF.rename('MSD',inplace=True)
            self.MSDDFDict[key]=MSDDF
            prev_nMols = nMols
    

    def calcMSDWindow(self,fractionalWindow=1):
        print('Calculating MSD.')
        ts = np.sort(np.array(list(set(self.molTimeSeriesDF.index.droplevel(1).to_numpy()))))
        self.MSDDFDict = {}

        for i in range(int(len(ts))-int(len(ts)*fractionalWindow)+1):
            sys.stdout.write('\rMSD calculation {:.2f}% complete'.format((i+1)*100.0/(int(len(ts))-int(len(ts)*fractionalWindow)+1)))
            j = i + int(len(ts)*fractionalWindow)-1
            molSlicedDF = self.molTimeSeriesDF[['x','y','z']].loc[ts[i]:ts[j]]
            comSlicedDF = self.comTimeSeriesDF[['x','y','z']].iloc[i:j+1]
            molSlicedDF[['xdiff','ydiff','zdiff']] = molSlicedDF[['x','y','z']]-molSlicedDF[['x','y','z']].loc[ts[i]]
            comSlicedDF[['xdiff','ydiff','zdiff']] = comSlicedDF[['x','y','z']]-comSlicedDF[['x','y','z']].iloc[0]
            print(comSlicedDF)
            molSlicedDF = molSlicedDF.reset_index(level=1)
            molSlicedDF[['xdiff','ydiff','zdiff']] = molSlicedDF[['xdiff','ydiff','zdiff']]-comSlicedDF[['xdiff','ydiff','zdiff']]
            molSlicedDF['SD'] = molSlicedDF['xdiff']**2 + molSlicedDF['ydiff']**2 + molSlicedDF['zdiff']**2
            molSlicedDF.set_index('mol',append=True,inplace=True)
            molSlicedDF.reset_index(level=0,inplace=True)
            molSlicedDF.rename(columns={'level_0':'timestep'},inplace=True)
            

            prev_nMols = 0
            for key in self.molDict:
                nMols = self.molDict[key]+prev_nMols
                slicedDF = molSlicedDF.loc[(molSlicedDF.index>prev_nMols)&(molSlicedDF.index<=nMols)]
                MSDDF = slicedDF['SD'].groupby(slicedDF['timestep']).mean()
                MSDDF.rename('MSD',inplace=True)
                MSDDF = pd.DataFrame(MSDDF)
                MSDDF['dt'] = MSDDF.index.to_numpy()-MSDDF.index[0]
                MSDDF.set_index('dt',inplace=True)
                if i == 0:
                    self.MSDDFDict[key]=MSDDF
                else:
                    self.MSDDFDict[key]['MSD']+=MSDDF['MSD']
                prev_nMols = nMols
        
        for key in self.molDict:
            self.MSDDFDict[key]['MSD'] = self.MSDDFDict[key]['MSD']/(int(len(ts))-int(len(ts)*fractionalWindow)+1)
        sys.stdout.write('\n')

    def calcDiffusivity(self):
        self.diffusivityDict = {}
        for key in self.MSDDFDict:
            MSDDF = self.MSDDFDict[key]
            x = MSDDF.index.to_numpy().reshape(-1,1)
            y = MSDDF.to_numpy()
            regression = LinearRegression(fit_intercept=False).fit(x,y)
            self.diffusivityDict[key]=regression.coef_[0][0]*1e-8*1e-8/self.timestepLength/1e-15/6
