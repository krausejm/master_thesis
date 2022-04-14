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
import ctypes
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

import arviz as az
import warnings
warnings.filterwarnings('ignore')
#this is my framework for bayesian fits and pretty plots
class bayesian_fit():
    def __init__(self,data,stanfile,nsamples=1000):
        self.data=data
        self.stan_data={'n':len(self.data['n'].values),
             'phi':list(self.data['phi'].values),
             'asym':list(self.data['asym'].values),
           'dasym':list(self.data['dasym'].values)}
        self.stanfile=stanfile
        self.x=self.data['phi']
        self.y=self.data['asym']
        self.dy=self.data['dasym']
        self.model = sp.CmdStanModel(stan_file=stanfile)
        self.fitobj = self.model.sample(data=self.stan_data,iter_sampling=nsamples)
    #funtion w/o phase
    def f(self,x,a):
            return a*np.cos(2*np.pi/180.*(-45-x))
    #function with phase
    def f1(self,x,a,b,c):
            return a*np.cos(x+c)+b
    #chi2 parameters of fitted data
    def get_chisqpars(self,function):
            popt, pcov = curve_fit(function,self.data['phi'].values,self.data['asym'].values,
                                   sigma=self.data['dasym'].values)
            return popt,pcov
    #print stan code
    def stan_code(self):
        return self.model.code()
    #get pd dataframe of samples
    def samples_df(self):
        return self.fitobj.draws_pd()
    #calculate waic for given samples and log likelihoods
    def waic(self):
        samples=self.fitobj.draws_pd()
        pwaic=np.array([np.var(samples[f'log_lik[{i+1}]'], ddof=1)for i in range(len(self.data['phi']))])
        lpd=np.array([np.log(np.mean(np.exp(samples[f'log_lik[{i+1}]'])))for i in range(len(self.data['phi']))])
        elpd=lpd-pwaic
        se=np.sqrt(len(self.x)*np.var(elpd,ddof=1))
        return sum(elpd), se, elpd
    #single method for pvalue
    def pval(self):
        samples=self.fitobj.draws_pd()
        y_arr=[samples[f'y_tilde[{i+1}]']for i in range(len(self.data['phi']))]
        pvalue=[len([y for y in y_arr[i] if y>=self.data['asym'][i]])/len(y_arr[i])for i in range(len(self.data['phi']))]
        return pvalue
    #plot posterior predictive distributions
    def plot_ppd(self,save,i,display,fitfunc,width=.5):
        samples=self.fitobj.draws_pd()
        #cosmetics
        fig, (ax, ax1)=plt.subplots(nrows=2,sharex=True,dpi=1000,gridspec_kw={'height_ratios': [4, 1]},num=1,clear=True)
        ax.set_ylabel('$p(A_{rep}|A)$')
        ax.set_xlim([-180,180])
        ax.xaxis.set_ticks_position("top")
        ax.grid(which='minor',color='grey',lw='.4')
        #get results from chisq fit and plot them
        xx=np.linspace(np.min(self.data['phi']),np.max(self.data['phi']),200)
        params=self.get_chisqpars(fitfunc)[0]
        yfit=fitfunc(xx,*params)
        ax.errorbar(x=self.data['phi'],y=self.data['asym'],yerr=self.data['dasym'],fmt='.',label='Data points',color='darkorange')
        #ax.plot(data['x'],data['y'],'x',label='Data points')
        ax.plot(xx,yfit,'r-',label='$\chi^2$ fit',color='darkorange')
        #violinplot for the replicated samples
        vp=ax.violinplot([samples[f'y_tilde[{i+1}]']for i in range(len(self.data['phi']))],np.array(self.x),
                         showmeans=False, showextrema=False, showmedians=False,widths=width)
        #cosmetics
        for b in vp['bodies']:
            # get the center
            m = np.mean(b.get_paths()[0].vertices[:, 0])
            # modify the paths to not go further right than the center
            b.get_paths()[0].vertices[:, 0] = np.clip(b.get_paths()[0].vertices[:, 0], m, np.inf)
            b.set_alpha(.5)
            b.set_color('midnightblue')

        #compute p value as measure of goodness of fit
        y_arr=[samples[f'y_tilde[{i+1}]']for i in range(len(self.data['phi']))]
        pvalue=[len([y for y in y_arr[i] if y>=self.data['asym'][i]])/len(y_arr[i])
                for i in range(len(self.data['phi']))]
        pvalue_lower=[len([y for y in y_arr[i] if y-self.data['dasym'][i]>=self.data['asym'][i]])/len(y_arr[i])
                for i in range(len(self.data['phi']))]
        pvalue_upper=[len([y for y in y_arr[i] if y+self.data['dasym'][i]>=self.data['asym'][i]])/len(y_arr[i])
                for i in range(len(self.data['phi']))]
        pval_errors=[np.abs(np.array(pvalue_lower)-np.array(pvalue)),np.abs(np.array(pvalue_upper)-np.array(pvalue))]
        #plot pvalue
        #ax1.grid(which='minor',color='grey',lw='.4')
        ax1.errorbar(x=self.data['phi'],y=pvalue,yerr=pval_errors,fmt='.',color='midnightblue')
        ax1.axhline(y=0.5, color='darkorange', linestyle='--',label='optimal value')
        ax1.set_xlabel('$\phi$ / deg')
        ax1.set_ylabel('$T(A_{rep}>A)$')
        ax1.set_ylim([-0.1,1.1])
        #cosmetics and legend
        plt.subplots_adjust(wspace=0, hspace=0)
        lines,labels = ax.get_legend_handles_labels() 
        patch = mpatches.Patch(color='midnightblue', label='$A_{rep}$',alpha=.5)
        lines.append(patch)
        lines1,labels1=ax1.get_legend_handles_labels()
        lines+=lines1
        tmp=lines[0]
        lines[0]=lines[1]
        lines[1]=tmp
        plt.legend(handles=lines,bbox_to_anchor=(1,2))
        fig.suptitle(f'toybin{i:04d}',y=1.1)
        if (save==True):
            plt.savefig(f'./ppd_checks/toybin{i:04d}.pdf',format='pdf',bbox_inches='tight',dpi=1000)
        if (display==True):
            plt.show()
    # plot posterior distributions of desired parameters, indicate chi2 fit value(s) and error(s)
    def plot_posterior(self,params,func,save=False):
        samples=self.fitobj.draws_pd()
        chi2pars=self.get_chisqpars(function=func)
        fig, (ax)=plt.subplots(ncols=1,nrows=len(params),dpi=1000)
        #print(len(params))
        
        #plot data
        ax.grid(False)
        ax.errorbar(x=chi2pars[0][i],y=0,xerr=np.sqrt(chi2pars[1][i,i]),label='$\chi^2$-fit',color='blue')
        #plot posterior kde distribution bc pretty
        sns.distplot(samples[params[i]],hist=True,kde=True,ax=ax,kde_kws={'bw':0.35})
        kde_curve = ax.lines[0]
        x = np.array(kde_curve.get_xdata())
        y = np.array(kde_curve.get_ydata())
        maxpos = y.argmax()
        maxx=x[maxpos]
        counts, bins = np.histogram(samples[params[i]],bins=20)
        mean=np.mean(samples[params[i]])
        sd=np.std(samples[params[i]],ddof=1)
        #mids = 0.5*(bins[1:] + bins[:-1])
        #probs = counts / np.sum(counts)
    
        #mean = np.sum(probs * mids)  
        #sd = np.sqrt(np.sum(probs * (mids - mean)**2))
        #indicate chi2 values
        #ax.axvline(x=maxx,ymin=0,ymax=1,linewidth='.5',label='MPV',color='red')
        ax.axvline(x=mean,ymin=0,ymax=1,linewidth='.5',label='$\mu\pm 1\sigma$',color='red',linestyle='dashed')
        ax.axvline(x=mean-sd,ymin=0,ymax=1,linewidth='.5',color='red',linestyle='dashed')
        ax.axvline(x=mean+sd,ymin=0,ymax=1,linewidth='.5',color='red',linestyle='dashed')
        #cosmetics and legend
        ax.legend(bbox_to_anchor=(1,1.2),fontsize=10)
        #cosmetics
        ax.set_xlabel('$\Sigma$',fontsize=10)
        ax.set_ylabel('Frequency',fontsize=10)
        ax.tick_params(axis='y', which='both',left=False,right=False,labelleft=False)
        ax.tick_params(axis='x', which='both',labelsize=10)
        ax.set_xlim([-1,1])
        #plt.subplots_adjust(hspace=1)
        if(save):
            fig.savefig(f'./posterior_{len(params)}_params.pdf',format='pdf',dpi=1000)

    # function to plot the posterior distributions as is     
    def plot_trace(self,params,save=False):
        samples=self.fitobj.draws_pd()
        fig, (ax)=plt.subplots(ncols=2,nrows=len(params),dpi=1000)
        #print(ax.shape)
        ls=['solid','dotted','dashed','dashdot']
        for i in range(len(params)):
            ax[i][0].grid(visible=False,axis='y')
            ax[i][0].set_xlabel(params[i],fontsize=10)
            ax[i][0].set_ylabel('Frequency',fontsize=10)
            ax[i][1].set_xlabel('iteration',fontsize=10)
            ax[i][1].set_ylabel(params[i],fontsize=10)
            for k in range(4):
                sns.kdeplot(samples[params[i]][k*1000:(k+1)*1000],ax=ax[i][0],linestyle=ls[k],color='black',
                            linewidth=.5,bw=0.35)
                ax[i][1].plot(np.arange(0,len(samples[params[i]])/4),
                              samples[params[i]][k*1000:(k+1)*1000],ls=ls[k],color='black',alpha=.7,linewidth=.5)
                for j in (0,1):
                    ax[i][j].tick_params(axis='y', which='both',left=False,right=False,labelleft=False)
                    #ax[i][j].set_title(params[i])
                    ax[i][j].tick_params(axis='both', which='major', labelsize=10)
        plt.subplots_adjust(hspace=1)
        #fig.set_size_inches(8.29,1*len(params))
def find_y_for_x(x, y, value):
    for i in range(len(x)):
        if x[i]<= value and x[i+1]>= value:
            return y[i]
            break
        
#calculate savage dickey density ratio = Bayes factor
def plot_sddr(samples,fixpar,saveas):#this is the model with one parameter more
    fig, ax =plt.subplots(dpi=1000)
    #define prior data and get posterior using kde
    x=np.linspace(-.5,.5,1000)
    def prior(x):
        return stats.norm.pdf(x,0,np.pi)
    priorval=prior(fixpar)
    sns.distplot(samples,hist=True,kde=True,ax=ax,kde_kws={'bw':0.35},label='$p(c|y,M_1)$')
    kde_curve = ax.lines[0]
    xx = kde_curve.get_xdata()
    y = kde_curve.get_ydata()
    postval=find_y_for_x(xx,y,fixpar)
    #plot
    ax.plot(x,prior(x),label='$\pi(c|M_1)$')
    ax.vlines(fixpar,0,np.max(y),linestyle='--',color='peachpuff',label='BF={:.2f}'.format(postval/priorval))
    ax.legend(bbox_to_anchor=(1,1))
    fig.savefig(saveas,format='pdf',dpi=1000,bbox_inches='tight')
def elpd_diff(fit,fit1):
    elpd=np.array(fit.waic()[2])
    elpd1=np.array(fit1.waic()[2])
    return np.abs(sum(elpd)-sum(elpd1)),np.sqrt(len(fit.x)*np.var(elpd-elpd1))
def fitfunc(x,sigma):
    return sigma*np.cos(2*(-45-x)*np.pi/180.)

def fit(nsamples,nbins,start): #define starting index
    cols=[f'toybin{i:04d}' for i in range(start,start+nbins)]
    index=['sigma_median','mcse','rhat']
    for i in range(12):
        index.append(f'pval{i:02d}')
    
    diagnostics_df=pd.DataFrame(columns=cols,index=index)
    sigma_df=pd.DataFrame(columns=cols)
    sigma_chi2_df=pd.DataFrame(columns=cols,index=['sigma','error'])

    for i in range(start,start+nbins):#no. of toy bins
        #read data
        df=pd.read_csv(f"toybins/toybin{i:04d}.txt",sep="\t")
        #df=pd.read_csv(f"new_toy_MC.txt",sep="\t")
        df.columns=['pol','setting','phi']
        #now write data to histo
        p45=df.loc[df['setting']==+45]
        m45=df.loc[df['setting']==-45]
        
        hp45=r.TH1F(f"hp45_{i}",";#phi;",12,-180,180)
        hm45=r.TH1F(f"hm45_{i}",";#phi;",12,-180,180)

        for k in p45['phi']:
            hp45.Fill(k)
        for k in m45['phi']:
            hm45.Fill(k)

        #now build event yield asymmetry
        #write histos to arrays
        p45_arr=[]
        m45_arr=[]
        for k in range(12):
            p45_arr.append(hp45.GetBinContent(k+1))
            m45_arr.append(hm45.GetBinContent(k+1))
            
        p45_arr=np.array(p45_arr)
        m45_arr=np.array(m45_arr)
        
        
        
        #first, normalize event yields
        norm_p_err=ctypes.c_double(0.)
        norm_m_err=ctypes.c_double(0.)
        
        norm_p=hp45.IntegralAndError(0,-1,norm_p_err,"")
        norm_m=hm45.IntegralAndError(0,-1,norm_m_err,"")
        
        p45_arr*=1./norm_p
        m45_arr*=1./norm_m
        #c=r.TCanvas()

        

        #enumerator of asymmetry
        enum=p45_arr-m45_arr
        #nominator of asymmetry
        nom=0.25*p45_arr+0.3*m45_arr
        final_e=[]
        for k in range(12):
            n_bot=hp45.GetBinContent(k+1)
            n_bot_err=hp45.GetBinError(k+1)
            #print(n_bot,n_bot_err)
            n_par=hm45.GetBinContent(k+1)
            n_par_err=hm45.GetBinError(k+1)
            pol_bot=0.3
            pol_par=0.25
            tmp_enum=((n_bot*n_par*norm_m*(pol_bot+pol_par)*norm_p_err.value)**2
                    +(n_bot*n_par*norm_p*(pol_bot+pol_par)*norm_m_err.value)**2
                    +(n_par*norm_p*norm_m*(pol_bot+pol_par)*n_bot_err)**2
                    +(n_bot*norm_p*norm_m*(pol_bot+pol_par)*n_par_err)**2)
            tmp_nom=np.power(pol_par*n_bot*norm_m+pol_bot*n_par*norm_p,4)
            final_e.append(np.sqrt(tmp_enum/tmp_nom))
        
        #enum.Divide(nom)
        #set error
        #save the data in a dataframe!
        #create arrays with y,dy
        asym=enum/nom
        phibins=[]
        for k in range(12):
            phibins.append(hp45.GetXaxis().GetBinCenter(k+1))
        
        stan_data={
            'n':len(asym), #no. of bins
            'phi':phibins, # bin centres of phi
            'asym':asym, #value of asymmetry 
            'dasym':final_e,#error of asymmetry
        }
        df_stan=pd.DataFrame(stan_data)
        
        #fit with mcmc and get samples
        fit=bayesian_fit(data=df_stan,stanfile='stan_file.stan',nsamples=nsamples)
        sigma=fit.samples_df()['sigma']
        #write results
        currbin=f"toybin{i:04d}"
        sigma_df[currbin]=sigma
        #write chi2 results
        popt, pcov = curve_fit(fitfunc,df_stan['phi'].values,df_stan['asym'].values,
                                   sigma=df_stan['dasym'].values)
        sigma_chi2_df[currbin]=[popt[0],np.sqrt(pcov[0][0])]
        #get fit diagnostics, i.e. pvalue
        pvalue=fit.pval()
        #get mcmc diagnostics, i.e. r_hat and median+mcse
        summary=fit.fitobj.summary()
        median=(summary['50%']['sigma'])
        mcse=(az.mcse(np.transpose(fit.fitobj.draws(concat_chains=False)[:,:,7]),method='median'))
        #print(fit.fitobj.draws_pd())
        rhat=(summary['R_hat']['sigma'])
        dlist=[median,mcse,rhat]
        dlist.extend(pvalue)
        #write diagnostics
        diagnostics_df[currbin]=dlist
        #make ppd plots
        #fit.plot_ppd(save=True,i=i,display=False,fitfunc=fit.f,width=50)
        print(f'fitted toyMC exp. no. {i}')

    return diagnostics_df, sigma_df, sigma_chi2_df
    #return sigma_chi2_df    
res=fit(nsamples=1000,nbins=5000,start=5000)
res[0].to_csv('diagnostics_1.csv')
res[1].to_csv('sigma_1.csv')
res[2].to_csv('sigma_chi2_1.csv')
