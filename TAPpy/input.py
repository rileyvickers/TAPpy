import numpy as np
import pandas as pd

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
        xBounds = np.array([int(i) for i in f.readline().split()[0:2]])
        yBounds = np.array([int(i) for i in f.readline().split()[0:2]])
        zBounds = np.array([int(i) for i in f.readline().split()[0:2]])
        
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

        for line in f:
            if line.startswith('Atoms'):
                break
        
        f.readline()

        for i in range(nAtoms):
            split_line = f.readline().split()
            domain.atomIDs[i] = int(split_line[0])
            domain.atomMolIDs[i] = int(split_line[1])
            domain.atomTypes[i] = int(split_line[2])
            domain.atomCharges[i] = float(split_line[3])
        
        domain.molIDs = list(set(domain.atomMolIDs))
        
        domain.molAtomIDs = []
        for i in range(len(domain.molIDs)):
            domain.molAtomIDs.append([])

        for atom, atomMol in zip(domain.atomIDs,domain.atomMolIDs):
            indexMolID = domain.molIDs.index(atomMol)
            domain.molAtomIDs[indexMolID].append(atom)




def readLAMMPSlog(domain,file):
    print('Reading LAMMPS log file.')
    with open(file) as f:
        for line in f:
            if line.startswith('units'):
                domain.simulationUnits = line.split()[1]
            elif line.startswith('timestep'):
                domain.timestepLength = float(line.split()[1])

    if not domain.timestepLength:
        print('Timestep length not defined in log file.')
    if not domain.simulationUnits:
        print('Units style not defined in log file.')
    

def readLAMMPStrajectory(domain,file):
    print('Reading LAMMPS trajectory file.')
    # TODO: add variant to read a series of files instead of 1
    nPreambleLines = 8
    domain.timesteps = []
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

    nTimesteps = len(domain.timesteps)
    # TODO: Add trajectory information to a molecule data structure
    molTimeSeriesArray = []
    comTimeSeriesArray = []
    with open(file) as f:
        for i in range(nTimesteps):
            atomTimestepArray = []
            for j in range(nPreambleLines+1):
                f.readline()
            for k in range(nAtoms):
                split_line = f.readline().split()
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
            comTimestepDF = atomTimestepDF.sum()
            for dimension in ['x','y','z']:
                molTimestepDF[dimension] = molTimestepDF[dimension+'mass']/molTimestepDF['mass']
                comTimestepDF[dimension] = comTimestepDF[dimension+'mass']/comTimestepDF['mass']
            molTimeSeriesArray.append(molTimestepDF[['x','y','z']])
            comTimeSeriesArray.append(comTimestepDF[['x','y','z']])
    
    domain.molTimeSeriesDF = pd.concat(molTimeSeriesArray,keys=domain.timesteps)
    domain.comTimeSeriesDF = pd.concat(comTimeSeriesArray,axis=1,keys=domain.timesteps).T

