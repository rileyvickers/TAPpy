import TAPpy
import numpy


testDirectory = './KCl/'
logFile = testDirectory + 'lammps.log'
trajectoryFile = testDirectory + 'dump.trj'
dataFile = testDirectory + 'data.lmp'
plotFile = testDirectory + 'plot.jpg'
molDict = {r'H$_2$O':2133,r'K$^+$':106,r'Cl$^-$':106}

testDomain = TAPpy.Domain(molDict=molDict)
testDomain.readInputs(dataFile,logFile,trajectoryFile)

# testDomain.calcDisplacement()
# testDomain.calcMSD()

testDomain.calcMSDWindow(fractionalWindow=0.9)

testDomain.calcDiffusivity()

TAPpy.output.plotDiffusivity(testDomain,plotFile)