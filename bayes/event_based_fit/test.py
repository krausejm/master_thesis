import cmdstanpy as sp
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
import matplotlib.patches as mpatches
import seaborn as sns
import arviz as az
import scipy.stats as stats
import ROOT as r
import time
plt.rcParams["xtick.minor.visible"] =  True
plt.rcParams["ytick.minor.visible"] =  True
plt.rcParams["mathtext.fontset"]="cm"
plt.rcParams['errorbar.capsize'] = 3
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.family']='serif'
plt.rcParams['font.size']=22
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.linestyle']=':'
plt.rcParams['grid.color']='black'
plt.rcParams['axes.prop_cycle'] = cycler(color=['black', 'red', 'blue', 'green'])


df=pd.read_csv("test.txt",sep="\t")
df.columns=['pol','phi','weight']
#these are prompt peak events
prmpt=df.loc[df['weight']==1]
prmpt=prmpt.reset_index(drop=True)
#sideband
side=df.loc[df['weight']!=1]
side=side.reset_index(drop=True)
nprmpt=len(prmpt)
total_nside=(len(side))
#weights used in the data, including 1--> exterminate 1
weights=pd.unique(df['weight'].values)
weights=np.array(weights[weights!=1])
#sideband events are identified by the corresponding weights
nside=np.array([len(df.loc[df['weight']==weights[i]])for i in range(7)])
#fraction of signal in prompt peak, note that the weights are already negative (TO BE CHECKED W/ FARAHS DATA!)
f=(nprmpt+np.sum(nside*weights))/nprmpt

stan_data={
    'N':nprmpt, #no. of prompt peak events and corresponding pol and phi values
    'phi_prmpt':list(prmpt['phi'].values),
    'pol_prmpt':list(prmpt['pol'].values),
    'M':total_nside,#no. of sideband events and..
    'phi_side':list(side['phi'].values),
    'pol_side':list(side['pol'].values),
    'f':f #fraction of signal in prmpt peak
}
model=sp.CmdStanModel(stan_file='toyMC_stan.stan')
model.compile()
fitobj=model.sample(data=stan_data,iter_sampling=1000)
#samples=fitobj.draws_pd()
#print(samples)


