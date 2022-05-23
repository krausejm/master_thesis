import cmdstanpy as sp
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler
import matplotlib.patches as mpatches
import arviz as az
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


def fit(nsamples,nbins,start): #define starting index
    cols=[f'toybin{i:04d}' for i in range(start,start+nbins)]
    diagnostics_df=pd.DataFrame(columns=cols,index=['sigma_median','mcse','rhat'])
    sigma_df=pd.DataFrame(columns=cols)
    sigma_2pi0_df=pd.DataFrame(columns=cols)
    for i in range(start,start+nbins):#no. of toy bins
        print(f"Fitting toy MC bin no.{i}")
        #read data
        df=pd.read_csv(f"toybins/toybin{i:04d}.txt",sep="\t")
        df_2pi0=pd.read_csv(f"toybins/2pi0bin.txt",sep='\t')
        #df=pd.read_csv(f"new_toy_MC.txt",sep="\t")
        df.columns=['pol','phi','weight']
        df_2pi0.columns=['sigma','dsigma','f_s','f_b']
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
            'f':f, #fraction of signal in prmpt peak
            'f_s': df_2pi0['f_s'].values[0], #fraction of etap events
            'f_b': df_2pi0['f_b'].values[0], #fraction of 2pi0 events
            'sigma_2pi0_meas': df_2pi0['sigma'].values[0],#measured value of sigma_2pi0 
            'dsigma_2pi0_meas': df_2pi0['dsigma'].values[0] #w/ stat error
        }
        print(nprmpt, total_nside)
        #now the stan model and mcmc
        model=sp.CmdStanModel(stan_file='toyMC_stan.stan')
        model.compile()
        fitobj=model.sample(data=stan_data,iter_sampling=nsamples,inits=0,output_dir='./stan_trash',show_progress=True)
        summary=fitobj.summary()
        samples=fitobj.draws_pd()
        #get mcmc diagnostics
        median=summary['50%']['sigma']
        mcse=(az.mcse(np.transpose(fitobj.draws(concat_chains=False)[:,:,7]),method='median'))
        rhat=(summary['R_hat']['sigma'])
        tmp_list=[median,mcse,rhat]
        currbin=f"toybin{i:04d}"
        diagnostics_df[currbin]=tmp_list
        sigma_df[currbin]=samples['sigma']
        sigma_2pi0_df[currbin]=samples['sigma_2pi0']
    return diagnostics_df, sigma_df, sigma_2pi0_df, summary

def fit_bin(nsamples,binnr): #fit only one bin
    #read data
    df=pd.read_csv(f"toybins/toybin{binnr:04d}.txt",sep="\t")
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
        'f':f, #fraction of signal in prmpt peak
        'f_s': f_s, #fraction of etap events
        'f_b': f_b, #fraction of 2pi0 events
        'sigma_2pi0_meas': sigma_2pi0,#measured value of sigma_2pi0 
        'dsigma_2pi0_meas': dsigma_2pi0 #w/ stat error
    }
    print(nprmpt, total_nside)
    #now the stan model and mcmc
    model=sp.CmdStanModel(stan_file='toyMC_stan.stan')
    model.compile()
    fitobj=model.sample(data=stan_data,iter_sampling=nsamples,inits=0)
    summary=fitobj.summary()
    samples=fitobj.draws_pd()
    return samples,summary

dfs=fit(nsamples=5000,nbins=1000,start=0)
diagnostics=dfs[0]
sigma=dfs[1]
sigma_2pi0=dfs[2]
diagnostics.to_csv('toy_diagnostics.csv')
sigma.to_csv('toy_sigma.csv')
sigma_2pi0.to_csv('toy_sigma2pi0.csv')


#sigma_df=pd.read_csv('toy_sigma.csv',index_col=0)
#diagnostics_df=pd.read_csv('toy_diagnostics.csv',index_col=0)

#rel_err=[(diagnostics_df[f'toybin{i:04d}']['mcse']/diagnostics_df[f'toybin{i:04d}']['sigma_median']) for i in range(len(diagnostics_df.columns))]
#plt.hist(rel_err,histtype='step',color='midnightblue')
#plt.xlabel('$\\frac{\Delta M(p(\Sigma|y))}{M(p(\Sigma|y)}$')
#plt.grid(False)
#plt.savefig('./plots/toyMC_mcse_hist.pdf',format='pdf',dpi=1000,bbox_inches='tight')

#all_sigmas=[]
#for i in range(len(sigma_df.columns)):
#    all_sigmas.extend((np.array(sigma_df[f'toybin{i:04d}'].values)))
