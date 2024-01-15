import numpy as np
import sys

cpdef double[:] calculate(double[:,:,:] molTimeSeries, double[:,:,:] comTimeSeries, double fractionalWindow):
    cdef Py_ssize_t nTimesteps = molTimeSeries.shape[0]
    cdef Py_ssize_t nMolecules = molTimeSeries.shape[1]
    cdef Py_ssize_t nDimensions = molTimeSeries.shape[2]
    cdef Py_ssize_t nWindow = np.intp(nTimesteps*fractionalWindow)
    cdef Py_ssize_t nMSD = nTimesteps-nWindow+1

    molSlicedArray = np.zeros((nWindow,nMolecules,nDimensions),dtype=np.double)
    comSlicedArray = np.zeros((nWindow, 1, nDimensions),dtype=np.double)
    molSumArray = np.zeros((nWindow,nMSD),dtype=np.double)
    MSD = np.zeros(nWindow,dtype=np.double)

    cdef double[:,:,:] molSlicedView = molSlicedArray
    cdef double[:,:,:] comSlicedView = comSlicedArray
    cdef double[:,:] molSumView = molSumArray
    cdef double[:] MSDview = MSD

    cdef Py_ssize_t i, j, ii, jj, kk, ll, prev_nMols, nMols
    cdef double disp, com_disp, val, SDperMolperTimestepSum, SDperWindowSum

    for i in range(nMSD):
        sys.stdout.write('\rMSD calculation {:.2f}% complete'.format((i+1)*100.0/nMSD))
        j = i + nWindow
        molSlicedView[:,:,:] = molTimeSeries[i:j,:,:]
        comSlicedView[:,:,:] = comTimeSeries[i:j,:,:]
        for ii in range(nWindow):
            SDperMolperTimestepSum = 0
            for kk in range(nDimensions):
                com_disp = comSlicedView[ii,0,kk] - comSlicedView[0,0,kk]
                for jj in range(nMolecules):
                    disp = molSlicedView[ii,jj,kk] - molSlicedView[0,jj,kk]
                    val = disp - com_disp
                    SDperMolperTimestepSum += val * val / nMolecules
            
            molSumView[ii,i] = SDperMolperTimestepSum
    
    for ii in range(nWindow):
        SDperWindowSum = 0 
        for i in range(nMSD):
            SDperWindowSum += molSumView[ii,i]
        MSDview[ii] = SDperWindowSum / nMSD
    sys.stdout.write('\n')
    return MSD

cpdef double[:,:] calculateMSDH(double[:,:,:] molTimeSeries, double[:,:,:] comTimeSeries, double[:,:] systemTimeSeries, double fractionalWindow):
    cdef Py_ssize_t nTimesteps = molTimeSeries.shape[0]
    cdef Py_ssize_t nMolecules = molTimeSeries.shape[1]
    cdef Py_ssize_t nDimensions = molTimeSeries.shape[2]
    cdef Py_ssize_t nSystemParams = 8
    cdef Py_ssize_t nWindow = np.intp(nTimesteps*fractionalWindow)
    cdef Py_ssize_t nMSD = nTimesteps-nWindow+1

    molSlicedArray = np.zeros((nWindow,nMolecules,nDimensions),dtype=np.double)
    comSlicedArray = np.zeros((nWindow, 1, nDimensions),dtype=np.double)
    systemSlicedArray = np.zeros((nWindow, nSystemParams),dtype=np.double)
    molSumArray = np.zeros((nMSD,nWindow),dtype=np.double)
    
    systemDeviationArray = np.zeros((nMSD,nWindow,nSystemParams),dtype=np.double)
    weightedMolSumArray  = np.zeros((nMSD,nWindow,nSystemParams),dtype=np.double)
    weightedSDperWindowSumArray = np.zeros(nSystemParams,dtype=np.double)

    MSD = np.zeros(nWindow,dtype=np.double)
    weightedMSD = np.zeros((nWindow,nSystemParams),dtype=np.double)

    cdef double[:,:,:] molSlicedView = molSlicedArray
    cdef double[:,:,:] comSlicedView = comSlicedArray
    cdef double[:,:] systemSlicedView = systemSlicedArray

    cdef double[:,:,:] systemDeviationView = systemDeviationArray
    cdef double[:,:] molSumView = molSumArray
    cdef double[:,:,:] weightedMolSumView = weightedMolSumArray
    cdef double[:] MSDview = MSD
    cdef double[:,:] weightedMSDview = weightedMSD

    cdef Py_ssize_t i, j, ii, jj, kk, ll, prev_nMols, nMols
    cdef double disp, com_disp, val, SDperMolperTimestepSum, SDperWindowSum
    cdef double[:] weightedSDperWindowSumView = weightedSDperWindowSumArray

    for i in range(nMSD):
        sys.stdout.write('\rMSD calculation {:.2f}% complete'.format((i+1)*100.0/nMSD))
        j = i + nWindow
        molSlicedView[:,:,:] = molTimeSeries[i:j,:,:]
        comSlicedView[:,:,:] = comTimeSeries[i:j,:,:]
        systemSlicedView[:,:] = systemTimeSeries[i:j,5:]        

        for ii in range(nWindow):

            for jj in range(ii+1):
                for kk in range(nSystemParams):

                    systemDeviationView[i,ii,kk] += systemSlicedView[jj,kk]
            
            for kk in range(nSystemParams):    
                systemDeviationView[i,ii,kk] = systemSlicedView[0,kk] - systemDeviationView[i,ii,kk]/(ii+1)

            SDperMolperTimestepSum = 0
            for kk in range(nDimensions):
                com_disp = comSlicedView[ii,0,kk] - comSlicedView[0,0,kk]
                for jj in range(nMolecules):
                    disp = molSlicedView[ii,jj,kk] - molSlicedView[0,jj,kk]
                    val = disp - com_disp
                    SDperMolperTimestepSum += val * val / nMolecules
                
            molSumView[i,ii] = SDperMolperTimestepSum

    for ii in range(nMSD):
        for jj in range(nWindow):
            for kk in range(nSystemParams):
                weightedMolSumView[ii,jj,kk] = molSumView[ii,jj] * systemDeviationView[ii,jj,kk]

    for ii in range(nWindow):
        SDperWindowSum = 0
        for j in range(nSystemParams):
            weightedSDperWindowSumView[j] = 0
        
        for i in range(nMSD):
            SDperWindowSum += molSumView[i,ii]
        
        for j in range(nSystemParams):
            for i in range(nMSD):
                weightedSDperWindowSumView[j] += weightedMolSumView[i,ii,j]


        MSDview[ii] = SDperWindowSum / nMSD
        for j in range(nSystemParams):
            weightedMSDview[ii,j] = weightedSDperWindowSumView[j] / nMSD
    sys.stdout.write('\n')
    combinedMSD = np.concatenate((MSD.reshape(MSD.shape[0],1),weightedMSD),axis=1)
    return combinedMSD