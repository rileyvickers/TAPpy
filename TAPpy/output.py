import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
plt.rcParams.update({'font.size': 16})

def func(x,a,b):
    return a * (x ** b)

def plotDiffusivity(domain,plotFile):
    nSubplots = len(domain.MSDDict)
    xSubplots = int(np.floor(np.sqrt(nSubplots)))
    ySubplots = int(np.ceil(nSubplots/xSubplots))
    fig,axes = plt.subplots(xSubplots,ySubplots, figsize=(ySubplots*5,xSubplots*5))
    try:
        axes = axes.flatten()
    except:
        axes =[axes]
    for key,ax in zip(domain.MSDDict,axes):
        MSD = domain.MSDDict[key]
        y = MSD.reshape(-1,1)
        x=[]
        for i in range(len(MSD)):
            x.append(i*(domain.timesteps[1]-domain.timesteps[0])*domain.timestepLength)
        x=np.array(x).reshape(-1,1)
        regression = LinearRegression(fit_intercept=False).fit(x,y)

        
        ax.plot(x,y,color='blue')
        ax.plot(x,regression.predict(x),color='red',ls='-.')
        ax.grid()
        diffusivity = regression.coef_[0][0]*1e-8*1e-8/1e-15/6
        Rsquared = regression.score(x,y)
        
        ax.annotate(r'R$^2$: {0:.2f}'.format(Rsquared),xy=(0,0.85),xycoords='axes fraction')
        # ax.annotate(r'MSD [$\AA ^2$] = {0:.3e} [$\AA ^2$/fs] $\times$ t [fs]'.format(regression.coef_[0]),xy=(0,0.80),xycoords='axes fraction')
        ax.annotate(r'D = {0:.3e} cm$^2$/s '.format(diffusivity),xy=(0,0.8),xycoords='axes fraction')

        ax.set_ylabel(r'MSD [$\AA ^2$]')
        ax.set_xlabel(r't [fs]')
        ax.set_title(key)
        fig.tight_layout()
        fig.savefig(plotFile,dpi=400)

        return diffusivity, Rsquared


def plotLinearDiffusivity(domain,plotFile,saveFile=None):
    nSubplots = len(domain.MSDDict)
    xSubplots = int(np.floor(np.sqrt(nSubplots)))
    ySubplots = int(np.ceil(nSubplots/xSubplots))
    fig,axes = plt.subplots(xSubplots,ySubplots, figsize=(ySubplots*5,xSubplots*5))
    try:
        axes = axes.flatten()
    except:
        axes =[axes]
    for key,ax in zip(domain.MSDDict,axes):
        MSD = domain.MSDDict[key]
        y = MSD.reshape(-1,1)
        x=[]
        for i in range(len(MSD)):
            x.append(i*(domain.timesteps[1]-domain.timesteps[0])*domain.timestepLength)
        x=np.array(x).reshape(-1,1)
        regression = LinearRegression(fit_intercept=True).fit(x[int(len(x)/2):],y[int(len(y)/2):])

        diffusivity = regression.coef_[0][0]*1e-8*1e-8/1e-15/6
        Rsquared = regression.score(x[int(len(x)/2):],y[int(len(y)/2):])
        
        ax.plot(x,y,color='blue')
        ax.plot(x[int(len(x)/2):],regression.predict(x[int(len(x)/2):]),color='red',ls='-.')
        ax.grid()
        
        ax.annotate(r'R$^2$: {0:.2f}'.format(Rsquared),xy=(0,0.85),xycoords='axes fraction')
        # ax.annotate(r'MSD [$\AA ^2$] = {0:.3e} [$\AA ^2$/fs] $\times$ t [fs]'.format(regression.coef_[0]),xy=(0,0.80),xycoords='axes fraction')
        ax.annotate(r'D = {0:.3e} cm$^2$/s '.format(diffusivity),xy=(0,0.8),xycoords='axes fraction')

        ax.set_ylabel(r'MSD [$\AA ^2$]')
        ax.set_xlabel(r't [fs]')
        ax.set_title(key)
        fig.tight_layout()
        fig.savefig(plotFile,dpi=400)
        
        return diffusivity, Rsquared

def plotSubdiffusivity(domain,plotFile,savefile=None):
    nSubplots = len(domain.MSDDict)
    xSubplots = int(np.floor(np.sqrt(nSubplots)))
    ySubplots = int(np.ceil(nSubplots/xSubplots))
    fig,axes = plt.subplots(xSubplots,ySubplots, figsize=(ySubplots*5,xSubplots*5))
    try:
        axes = axes.flatten()
    except:
        axes =[axes]
    for key,ax in zip(domain.MSDDict,axes):
        MSD = domain.MSDDict[key]
        y = MSD
        x=[]
        for i in range(len(MSD)):
            x.append(i*(domain.timesteps[1]-domain.timesteps[0])*domain.timestepLength)
        x=np.array(x)
        regression,regressionCov = curve_fit(func,x,y)
        yPred = func(x,regression[0],regression[1])
        regressionErr = np.sqrt(np.diag(regressionCov))
        ax.plot(x,y,color='blue')
        ax.plot(x,yPred,color='red',ls='-.')
        ax.grid()

        beta = regression[0]
        alpha = regression[1]
        RMSE = np.sqrt(((y-yPred) * (y-yPred)).sum()/len(y))

        
        ax.annotate(r'RMSE: {0:.2f}'.format(RMSE) + r' $\AA^2$',xy=(0,0.85),xycoords='axes fraction')
        # ax.annotate(r'MSD [$\AA ^2$] = {0:.3e} [$\AA ^2$/fs] $\times$ t [fs]'.format(regression.coef_[0]),xy=(0,0.80),xycoords='axes fraction')
        ax.annotate(r'$\beta$ =' + '{0:.3e}'.format(beta) + r' $\AA^2$/fs$^\alpha$',xy=(0,0.8),xycoords='axes fraction')
        ax.annotate(r'$\alpha$ =' + ' {0:.3f}'.format(alpha),xy=(0,0.75),xycoords='axes fraction')
        ax.set_ylabel(r'MSD [$\AA ^2$]')
        ax.set_xlabel(r't [fs]')
        ax.set_title(key)
        fig.tight_layout()
        fig.savefig(plotFile,dpi=400)

        return beta,alpha,RMSE