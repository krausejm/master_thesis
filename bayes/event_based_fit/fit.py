import cmdstanpy as sp
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import arviz as az
import scipy.stats as stats
import arviz as az

def fit(nsamples):
    cols=[f'ebin{i:02d}costbin{j:02d}' for i in range(11) for j in range(12)]
    diagnostics_df=pd.DataFrame(columns=cols,index=['sigma_median','mcse','rhat'])
    sigma_df=pd.DataFrame(columns=cols)
    for i in range(11):#no. of energy bins
        for j in range(12):#no. of angle bins
            #read data
            df=pd.read_csv(f"ebin{i:02d}/ebin{i:02d}costbin{j:02d}.txt",sep="\t")
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
            nside=np.array([len(df.loc[df['weight']==weights[i]])for i in range(len(weights))])
            #fraction of signal in prompt peak, note that the weights are already negative (TO BE CHECKED W/ FARAHS DATA!)
            f=(nprmpt+np.sum(nside*weights))/nprmpt
            #print(f)
            stan_data={
                'N':nprmpt, #no. of prompt peak events and corresponding pol and phi values
                'phi_prmpt':list(prmpt['phi'].values),
                'pol_prmpt':list(prmpt['pol'].values),
                'M':total_nside,#no. of sideband events and..
                'phi_side':list(side['phi'].values),
                'pol_side':list(side['pol'].values),
                'f':f #fraction of signal in prmpt peak
            }
            #now the stan model and mcmc
            model=sp.CmdStanModel(stan_file='toyMC_stan.stan')
            model.compile()
            fitobj=model.sample(data=stan_data,iter_sampling=nsamples,inits=0,output_dir='stan_trash')
            summary=fitobj.summary()
            samples=fitobj.draws_pd()
            #get mcmc diagnostics
            median=summary['50%']['sigma']
            mcse=(az.mcse(np.transpose(fitobj.draws(concat_chains=False)[:,:,7]),method='median'))
            rhat=az.rhat(np.transpose(fitobj.draws(concat_chains=False)[:,:,7]))
            tmp_list=[median,mcse,rhat]
            currbin=f"ebin{i:02d}costbin{j:02d}"
            print(currbin)
            diagnostics_df[currbin]=tmp_list
            sigma_df[currbin]=samples['sigma']
    return diagnostics_df, sigma_df
        
def fit_bin(nsamples,i,j): #fit only one bin
    #read data
    df=pd.read_csv(f"ebin{i:02d}/ebin{i:02d}costbin{j:02d}.txt",sep="\t")    
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
    nside=np.array([len(df.loc[df['weight']==weights[i]])for i in range(len(weights))])
    #fraction of signal in prompt peak, note that the weights are already negative (TO BE CHECKED W/ FARAHS DATA!)
    f=(nprmpt+np.sum(nside*weights))/nprmpt
    print(f)
    stan_data={
        'N':nprmpt, #no. of prompt peak events and corresponding pol and phi values
        'phi_prmpt':list(prmpt['phi'].values),
        'pol_prmpt':list(prmpt['pol'].values),
        'M':total_nside,#no. of sideband events and..
        'phi_side':list(side['phi'].values),
        'pol_side':list(side['pol'].values),
        'f':f #fraction of signal in prmpt peak
    }
    print(nprmpt, total_nside)
    #now the stan model and mcmc
    model=sp.CmdStanModel(stan_file='toyMC_stan.stan')
    model.compile()
    fitobj=model.sample(data=stan_data,iter_sampling=nsamples,inits=0)
    summary=fitobj.summary()
    samples=fitobj.draws_pd()
    return samples,summary

dfs=fit(nsamples=5000)
diagnostics=dfs[0]
sigma=dfs[1]
diagnostics.to_csv('new_diagnostics.csv')
sigma.to_csv('new_sigma.csv')
