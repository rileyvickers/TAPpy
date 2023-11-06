import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
plt.rcParams.update({'font.size': 16})

def plotDiffusivity(domain,plotFile):
    nSubplots = len(domain.MSDDFDict)
    xSubplots = int(np.floor(np.sqrt(nSubplots)))
    ySubplots = int(np.ceil(nSubplots/xSubplots))
    fig,axes = plt.subplots(xSubplots,ySubplots, figsize=(ySubplots*5,xSubplots*5))
    try:
        axes = axes.flatten()
    except:
        axes =[axes]
    for key,ax in zip(domain.MSDDFDict,axes):
        MSDDF = domain.MSDDFDict[key]
        x = MSDDF.index.to_numpy().reshape(-1,1)
        y = MSDDF.to_numpy()
        regression = LinearRegression(fit_intercept=False).fit(x,y)

        
        ax.plot(x,y)
        ax.plot(x,regression.predict(x))
        ax.grid()
        
        ax.annotate(r'R$^2$: {0:.2f}'.format(regression.score(x,y)),xy=(0,0.85),xycoords='axes fraction')
        # ax.annotate(r'MSD [$\AA ^2$] = {0:.3e} [$\AA ^2$/fs] $\times$ t [fs]'.format(regression.coef_[0]),xy=(0,0.80),xycoords='axes fraction')
        ax.annotate(r'D = {0:.3e} cm$^2$/s '.format(regression.coef_[0][0]*1e-8*1e-8/1e-15/6),xy=(0,0.8),xycoords='axes fraction')

        ax.set_ylabel(r'MSD [$\AA ^2$]')
        ax.set_xlabel(r't [fs]')
        ax.set_title(key)
        fig.tight_layout()
        fig.savefig(plotFile,dpi=400)