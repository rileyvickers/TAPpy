import numpy as np
import pandas as pd
import sys
import os
import re 

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def helloworld():
    print('Sup fuckers')

def readLAMMPSdata(domain,file):
    print('Reading LAMMPS data file.')
    with open(file) as f:
        for i in range(2):
            f.readline()
        nAtoms = int(f.readline().split()[0])
        nAtomTypes = int(f.readline().split()[0])
        nBonds = int(f.readline().split()[0])
        nBondTypes = int(f.readline().split()[0])
        nAngles = int(f.readline().split()[0])
        nAtngleTypes = int(f.readline().split()[0])
        nDihedrals = int(f.readline().split()[0])
        nDihedralTypes = int(f.readline().split()[0])
        nImpropers = int(f.readline().split()[0])
        nImproperTypes = int(f.readline().split()[0])

        domain.atomIDs = np.zeros(nAtoms) 
        domain.atomTypes = np.zeros(nAtoms)
        domain.atomCharges = np.zeros(nAtoms)
        domain.atomMolIDs = np.zeros(nAtoms)

        domain.atomTypeMasses = np.zeros(nAtomTypes)
        domain.atomTypeEpsilons = np.zeros(nAtomTypes)
        domain.atomTypeSigmas = np.zeros(nAtomTypes)

        f.readline()
        xBounds = np.array([float(i) for i in f.readline().split()[0:2]])
        yBounds = np.array([float(i) for i in f.readline().split()[0:2]])
        zBounds = np.array([float(i) for i in f.readline().split()[0:2]])
        
        domain.boundaries[0] = xBounds
        domain.boundaries[1] = yBounds
        domain.boundaries[2] = zBounds

        for i in range(3):
            f.readline()
        
        for i in range(nAtomTypes):
            split_line = f.readline().split()
            domain.atomTypeMasses[i] = float(split_line[1])

        for i in range(3):
            f.readline()

        for i in range(nAtomTypes):
            split_line = f.readline().split()
            domain.atomTypeEpsilons[i] = float(split_line[1])
            domain.atomTypeSigmas[i] = float(split_line[2])
        
        # TODO: Potentially not necessary, 
        # read bond, angle, dihedral, improper information into data

        # for line in f:
        #     if line.startswith('Atoms'):
        #         break
        
        # f.readline()

        # for i in range(nAtoms):
        #     split_line = f.readline().split()
        #     domain.atomIDs[i] = int(split_line[0])
        #     domain.atomMolIDs[i] = int(split_line[1])
        #     domain.atomTypes[i] = int(split_line[2])
        #     domain.atomCharges[i] = float(split_line[3])
        
        # domain.molIDs = list(set(domain.atomMolIDs))
        
        # domain.molAtomIDs = []
        # for i in range(len(domain.molIDs)):
        #     domain.molAtomIDs.append([])

        # for atom, atomMol in zip(domain.atomIDs,domain.atomMolIDs):
        #     indexMolID = domain.molIDs.index(atomMol)
        #     domain.molAtomIDs[indexMolID].append(atom)




def readLAMMPSlog(domain,file):
    print('Reading LAMMPS log file.')
    domain.systemTimeSeriesArray = []
    with open(file) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith('units'):
                domain.simulationUnits = line.split()[1]
            elif line.startswith('timestep'):
                domain.timestepLength = float(line.split()[1])
            elif line.startswith('run'):
                linesplit = line.split()
                numTimesteps = int(linesplit[1])
                if numTimesteps > 0:
                    while True:
                        line = f.readline().split()                  
                        if not line:
                            break
                        elif line[0] == 'Step':
                            domain.systemTimeSeriesArrayColumns = line
                            while True:
                                line = f.readline().split()
                                if not line:
                                    break

                                elif line[0] == 'Loop':
                                    break
                                
                                elif len(line) != len(domain.systemTimeSeriesArrayColumns):
                                    break    

                                domain.systemTimeSeriesArray.append([float(i) for i in line])
                        elif line[0] == 'Loop':
                            break
    domain.systemTimeSeries = np.array(domain.systemTimeSeriesArray)
    if not domain.timestepLength:
        print('Timestep length not defined in log file.')
    if not domain.simulationUnits:
        print('Units style not defined in log file.')
    

def readLAMMPStrajectory(domain,file,trim=None):
    print('Reading LAMMPS trajectory file.')
    # TODO: add variant to read a series of files instead of 1
    nPreambleLines = 8
    atomIDcounter=0
    domain.timesteps = []
    if os.path.isfile(file):
        with open(file) as f:
            i=0
            for line in f:
                if line == 'ITEM: TIMESTEP\n':
                    domain.timesteps.append(int(f.readline()))
                elif line == 'ITEM: NUMBER OF ATOMS\n':
                    nAtoms = int(f.readline())
                elif line.startswith('ITEM: ATOMS'):
                    if i == 0:
                        header = line.split()[2:]
                        
                        idIndex = header.index('id')
                        molIndex = header.index('mol')
                        typeIndex = header.index('type')
                        xIndex = header.index('xu')
                        yIndex = header.index('yu')
                        zIndex = header.index('zu')
                        i+=1
                    else:
                        headerCheck = line.split()[2:]
                        if headerCheck == header:
                            continue
                        else:
                            print('Trajectory header changes mid-simulation: unexpected.')
    elif os.path.isdir(file):
        trajFiles = [i for i in sorted(os.listdir(file),key=natural_key) if i.endswith('.trj')]
        domain.timesteps = [int(re.sub(r'[^0-9]','',i)) for i in trajFiles]
        readFile = file+'/'+trajFiles[0]
        with open(readFile) as f:
            i=0
            for line in f:
                if line == 'ITEM: NUMBER OF ATOMS\n':
                    nAtoms = int(f.readline())
                elif line.startswith('ITEM: ATOMS'):
                    if i == 0:
                        header = line.split()[2:]
                        
                        idIndex = header.index('id')
                        molIndex = header.index('mol')
                        typeIndex = header.index('type')
                        xIndex = header.index('xu')
                        yIndex = header.index('yu')
                        zIndex = header.index('zu')
                        i+=1

    domain.timesteps = np.array(domain.timesteps)
    nTimesteps = len(domain.timesteps)
    # TODO: Add trajectory information to a molecule data structure
    molTimeSeriesArray = []
    comTimeSeriesArray = []
    if os.path.isfile(file):
        with open(file) as f:
            for i in range(nTimesteps):
                sys.stdout.write('\rTrajectory read {:.2f}% complete'.format((i+1)*100.0/nTimesteps))
                atomTimestepArray = []
                for j in range(nPreambleLines+1):
                    f.readline()
                atomIDcounter = 0
                for k in range(nAtoms):
                    
                    split_line = f.readline().split()
                    if int(split_line[molIndex]) == 0:
                        molID = int(split_line[idIndex])
                        atomIDcounter += 1 
                    else:
                        molID = int(split_line[molIndex])
                    atomType = int(split_line[typeIndex])
                    atomX = float(split_line[xIndex])
                    atomY = float(split_line[yIndex])
                    atomZ = float(split_line[zIndex])

                    atomTimestepArray.append([molID,domain.atomTypeMasses[atomType-1],atomX,atomY,atomZ])
                    
                atomTimestepDF = pd.DataFrame(atomTimestepArray,columns=['mol','mass','x','y','z'])

                for dimension in ['x','y','z']:
                    atomTimestepDF[dimension+'mass'] = atomTimestepDF[dimension] * atomTimestepDF['mass']
                molTimestepDF = atomTimestepDF.groupby('mol').sum()
                domain.molIDs = molTimestepDF.index.values
                comTimestepDF = atomTimestepDF.sum()
                for dimension in ['x','y','z']:
                    molTimestepDF[dimension] = molTimestepDF[dimension+'mass']/molTimestepDF['mass']
                    comTimestepDF[dimension] = comTimestepDF[dimension+'mass']/comTimestepDF['mass']
                molTimeSeriesArray.append(molTimestepDF[['x','y','z']])
                comTimeSeriesArray.append(comTimestepDF[['x','y','z']])
    elif os.path.isdir(file):
        for i in range(nTimesteps):
            readFile = file+'/'+trajFiles[i]
            with open(readFile,'rb') as f:
                sys.stdout.write('\rTrajectory read {:.2f}% complete'.format((i+1)*100.0/nTimesteps))
                atomTimestepArray = []
                for j in range(nPreambleLines+1):
                    f.readline()
                atomIDcounter = 0
                for k in range(nAtoms):
                    split_line = f.readline().split()
                    if int(split_line[molIndex]) == 0:
                        molID = int(split_line[idIndex])
                        atomIDcounter += 1 
                    else:
                        molID = int(split_line[molIndex])
                    atomType = int(split_line[typeIndex])
                    atomX = float(split_line[xIndex])
                    atomY = float(split_line[yIndex])
                    atomZ = float(split_line[zIndex])

                    atomTimestepArray.append([molID,domain.atomTypeMasses[atomType-1],atomX,atomY,atomZ])
                  
                atomTimestepDF = pd.DataFrame(atomTimestepArray,columns=['mol','mass','x','y','z'])

                for dimension in ['x','y','z']:
                    atomTimestepDF[dimension+'mass'] = atomTimestepDF[dimension] * atomTimestepDF['mass']
                molTimestepDF = atomTimestepDF.groupby('mol').sum()
                domain.molIDs = molTimestepDF.index.values
                comTimestepDF = atomTimestepDF.sum()
                for dimension in ['x','y','z']:
                    molTimestepDF[dimension] = molTimestepDF[dimension+'mass']/molTimestepDF['mass']
                    comTimestepDF[dimension] = comTimestepDF[dimension+'mass']/comTimestepDF['mass']
                molTimeSeriesArray.append(molTimestepDF[['x','y','z']])
                comTimeSeriesArray.append(comTimestepDF[['x','y','z']])
    sys.stdout.write('\n')
    if atomIDcounter > 0:
        print('Used atom ID to define mol ID {0} times.'.format(atomIDcounter))
    domain.molTimeSeries = pd.concat(molTimeSeriesArray,keys=domain.timesteps)
    domain.molTimeSeries = domain.molTimeSeries.values.reshape(len(domain.timesteps),len(domain.molIDs),3)
    domain.comTimeSeries = pd.concat(comTimeSeriesArray,axis=1,keys=domain.timesteps).T.values.reshape(len(domain.timesteps),1,3)

    if trim:
        # disp = ((domain.molTimeSeries[-1]-domain.molTimeSeries[0])**2).sum(axis=1)
        # indices = (((domain.molTimeSeries[-1]-domain.molTimeSeries[0])**2).sum(axis=1) < meanDisp).nonzero()
        # domain.molIDs = (((domain.molTimeSeries[-1]-domain.molTimeSeries[0])**2).sum(axis=1) >= meanDisp).nonzero()
        # indices = np.argpartition(disp,int(0.9 * len(disp))[:int(0.9 * len(disp))])
        disp = ((domain.molTimeSeries[-1]-domain.molTimeSeries[0])**2).sum(axis=1)
        indices = np.concatenate(((disp).argsort()[:int(np.ceil(len(disp)*trim)/2)], (disp).argsort()[int(np.ceil(len(disp)-len(disp)*trim/2)):]))
        # indices = (disp).argsort()[int(np.ceil(len(disp)*trim)/2):int(np.ceil(len(disp)-len(disp)*trim/2))]
        domain.molIDs= (disp).argsort()[int(np.ceil(len(disp)*trim)/2):int(np.ceil(len(disp)-len(disp)*trim/2))]
        print(len(indices))
        print(len(domain.molIDs))
        domain.molTimeSeries = np.delete(domain.molTimeSeries,indices,axis=1)
    print('Analyzing {} molecular trajectories.'.format(len(domain.molIDs)))
        

