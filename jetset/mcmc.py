import numpy as np
import scipy as sp
from scipy import stats

class MCMC (object):

    def __init__(self,fit_par,range_rel,range_sigma,jump_step_rel,jump_step_sigma,proposal='uniform'):
        
        self.burn_in_size
        
        self.trials_size
        
        self.fit_par=fit_par

        self.par_range
        
        self.jump_step
        
        if proposal=='uniform':
            self.proposal_prob=stats.uniform.pdf
            self.proposal_val=stats.uniform.rvs
        elif proposal=='normal':
            self.proposal_prob=stats.norm.pdf
            self.proposal_val=stats.norm.rvs
        else:
            print "proposal distribution not in allowed"
            
            return     
        
    
    
    def run(self):
        
        LEN=0
        
        TRIAL=0

        oldLike=self.func_stat()
        
        self.init_par()
        
        while TRIAL <=self.MCMC_TRIALS and LEN<self.Lenght:
            
            TRIAL+=1
                       
            #Draw parameters 
            self.random_walk()
          
            self.update_prior()
        
            if self.prior21>0:
                
                self.UpdatePrefac()
        
                newLike=self.func_stat()
        
                ratio=np.exp(newLike-oldLike)*(self.Prefac)
                
                aprob=min([1.0,ratio])
        
                u = stats.uniform.rvs(loc=0.0,scale=1.0)
      
                #Decide trial acceptance
                if  u < aprob:
                    LEN+=1
                    print"TRIAL=%d accetto ratio=%e u=%e par[1]=%s -LogLike=%e "%(TRIAL,ratio,u,self.par[1],-newLike)
                    oldLike=newLike
          
    
    def init_MCMC_pars(self):
        for pi in xrange(len(self.fit_par)): 
            if self.fit_par[pi].err is not None:
                delta_par=self.fit_par[pi].err*self.range_sigma
                step_par=self.fit_par[pi].err*self.jump_step_sigma
            else:
                delta_par=self.fit_par[pi].val*self.range_rel
                step_par=self.fit_par[pi].val*self.jump_step_rel
                
            self.fit_par[pi].MCMC_range_min=self.fit_par[pi].val - delta_par
            self.fit_par[pi].MCMC_range_max=self.fit_par[pi].val + delta_par
            self.fit_par[pi].MCMC_jump_step=step_par
            self.fit_par[pi].MCMC_val=self.fit_par[pi].val
            self.fit_par[pi].MCMC_trial=None
            self.fit_par[pi].MCMC_sample=sp.zeros(self.trials_size)
    
    
    def random_walk(self,size):
        for pi in xrange(len(self.fit_par)): 
            self.fit_par[pi].MCMC_trial=self.draw_proposal(self.parfit_par[pi])
           

            
    def draw_proposal(self,par):
        loc=par.MCMC_range_min
        scale=par.MCMC_range_max-par.MCMC_range_min
        return self.proposal_val(loc=loc,scale=scale)
    
    
    def update_prior(self):
       
        self.prior12=1
        self.prior21=1
        
        for pi in xrange(len(self.fit_par)): 
    
            self.prior12 =self.prior12 * self.prior_fact(self.fit_par[pi],self.fit_par[pi].MCMC_val,self.fit_par[pi].MCMC_trial)

            self.prior21 =self.prior21 * self.prior_fact(self.fit_par[pi],self.fit_par[pi].MCMC_trial,self.fit_par[pi].MCMC_val)    
    
    
    def prior_fact(self,_prime,x,scale):
        
        loc=par.MCMC_range_min
        scale=par.MCMC_range_max-par.MCMC_range_min
        return self.proposal_prob(x,loc=loc,scale=scale) 