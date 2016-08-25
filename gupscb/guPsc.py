
# coding: utf-8

# # MCMC analysis of GU Psc b

# In[1]:

#import modules and enable inline plots
import numpy as np
from scipy import stats #For Kernel Density Estimation
from scipy.linalg.lapack import dpotrf
from scipy.linalg.lapack import dpotrs
from scipy.optimize import leastsq
import matplotlib  #ploting
#matplotlib.use("Agg") 
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import ScalarFormatter
#get_ipython().magic('matplotlib inline')


# In[2]:

import mcmcroutines as mcmc


# ## read in data

# In[3]:

import sys
if len(sys.argv) < 2 :
    print('Usage: FileNumber, StarNumber')
    print(' ')
    print(' Filenumber 1-3: ')
    print('  1: Gu_variab_20141011_J_normlightcurve_GU_comp.txt')
    print('  2: Gu_variab_20141010_J_normlightcurve_GU_comp.txt')
    print('  3: Gu_variab_20131222_J_normlightcurve_GU_comp.txt')
    print(' ')
    print(' Star Number: ')
    print('  1: GU Psc b')
    print('  2-8: comparision stars')
    print(' ')
    sys.exit()


nf=int(sys.argv[1])
nd=int(sys.argv[2])

if nf == 1:
    filename='Gu_variab_20141011_J_normlightcurve_GU_comp.txt'
elif nf == 2:
    filename='Gu_variab_20141010_J_normlightcurve_GU_comp.txt'
elif nf == 3:
    filename='Gu_variab_20131222_J_normlightcurve_GU_comp.txt'
else:
    print('Filenumber is out of range [1-3]',nf)

data=[]
#nd=1 #data column - range is 1-8 with 1=guPscb
f = open(filename, 'r')
for line in f:
    line = line.strip() #get rid of the \n at the end of the line
    columns = line.split() #break into columns
    #data.append([float(i) for i in columns])
    data.append([float(columns[0]),float(columns[nd]),float(columns[9]),float(columns[10]),float(columns[11])])
f.close()
data=np.array(data)
data[:,2]=data[:,2]-np.mean(data[:,2])#+1.0 #remove offset from ELCs
data[:,3]=data[:,3]-np.mean(data[:,3])#+1.0
data[:,4]=data[:,4]-np.mean(data[:,4])#+1.0


# In[4]:

plt.figure(figsize=(10,7)) 
plt.plot(data[:,0],data[:,1],color='b',label='Target Star')
plt.plot(data[:,0],data[:,2],color='r',label='ELC-1')
plt.plot(data[:,0],data[:,3],color='g',label='ELC-2')
plt.plot(data[:,0],data[:,4],color='orange',label='ELC-3')
plt.legend()
plt.show()


# ## Our Model (cosine + PCA)

# In[5]:

def model(pars,data):
    "Our Model"
    
    #Fitted Model parameters
    a=pars[0]       #amplitude of cosine (relative flux)
    p=pars[1]       #period of cosine (hours)
    t0=pars[2]     #phase offset (radians)
    bk=np.zeros(3)  #ELC scales
    bk[0]=pars[3]     #ELC-1
    bk[1]=pars[4]     #ELC-2
    bk[2]=pars[5]     #ELC=3
    mean=pars[6]    #Mean (zero point offset)
    
    #Fixed Model parameters
    #p=4.7 #period of cosine (days)
    tpi=2.0*np.pi #two-pi
    
    m=a*np.cos(tpi/p*(data[:,0]-t0))
    for i in range(len(bk)):
        m=m+bk[i]*data[:,i+2]
    m=m+mean
    
    return m;


# ## our likelihood model for uncorrelated noise

# In[6]:

def loglikelihood(func,pars,data):
    "log-likeihood function"
    
    #The next two are part of the noise-model (see the log-likelihood function)
    sig=pars[7]     #point-to-point scatter
    
    m=func(pars,data) #get model
    n=len(data[:,1])
    if n < 1:
        ll=-1.0e30  #set bad value for no data.
    else:
        ll=-0.5*(n*np.log(2*np.pi)+n*(np.log(sig*sig))         +sum((m-data[:,1])*(m-data[:,1])/(sig*sig)))
    
    return ll;
        


# ## Our Prior

# In[8]:

def lprior(pars):
    "Simple prior with log(pr)=0 with valid parameters, otherwise return large small value"
        
    badlpr=-1.0e30    #bad lpr value
    lpr=0.0           #default return value
    
    if pars[0] < 0: #we want a positive amplitude
        lpr=badlpr
        
    if pars[1] > 10:  #Keep period between 0 and 10 hours
        lpr=badlpr
    if pars[1] < 0:
        lpr=badlpr
        
    if pars[2] > 2*np.pi:  #Keep phi betwee 0 and 2*pi
        lpr=badlpr
    if pars[2] < 0:
        lpr=badlpr
        
    if pars[3] > 10:  #broad priors for Bk
        lpr=badlpr
    if pars[3] < -10:
        lpr=badlpr
    if pars[4] > 10:
        lpr=badlpr
    if pars[4] < -10:
        lpr=badlpr
    if pars[5] > 10:
        lpr=badlpr
    if pars[5] < -10:
        lpr=badlpr
        
    if pars[6] > 10:  #broad prior on mean
        lpr=badlpr
    if pars[6] < -10:
        lpr=badlpr
        
    if pars[7] > 10:  #point-to-point scatter - must be positive
        lpr=badlpr
    if pars[7] < 0:
        lpr=badlpr
    
    #this last parameter is only used for the correlated noise model
    if len(pars) > 8:
        if pars[8] > 10:  #correlated noise amplitude - must be positive
            lpr=badlpr
        if pars[8] < 0:
            lpr=badlpr
        
        
    return lpr;


# ## Initial Guess for Parameters and get beta for MCMC

# In[9]:

#       A  period  t0  b1    b2   b3   mean   sig   
#       0    1     2   3     4     5    6      7    
label=['A','Per','t0','B1','B2','B3','mean','sig']
colour=['r','yellow','b','g','purple','brown','black','orange']
pars=[0.021, 5.854, 4.000, 0.301,-0.215,-0.080, 1.000, 0.041]
#pars=[0.021, 5.854, 1.806, 0.301,-0.215,-0.080, 0.000, 0.041]
beta=[0.005, 0.500, 0.100, 0.100, 0.100, 0.100, 0.005, 0.010]
niter=2000    #number of chains to generate for testing acceptance rates
burnin=200
corscale=mcmc.betarescale(pars,data,beta,niter,burnin,model,loglikelihood,lprior,mcmc.mhgmcmc,imax=40)
betanew=beta*corscale #apply our new beta
print(betanew)


# Use Metropolis-Hastings to create a buffer for DEMCMC

# In[10]:

betanew=beta*corscale #apply our new beta
nbuffer=10000 #Size of our deMCMC buffer
burnin=500   #burn-in for M-H-G
niter=nbuffer+burnin 
chain,accept=mcmc.genchain(pars,data,betanew,niter,model,loglikelihood,lprior,mcmc.mhgmcmc)
buffer=np.copy(chain[burnin:,:])


# Generate Chains with DEMCMC

# In[11]:

corbeta=1.0 #scale for correlated jumps
parin=buffer[nbuffer-1,:] #we can start with the last state from our buffer.
niter=50000
chain,accept=mcmc.genchain(parin,data,betanew,niter,model,loglikelihood,lprior,mcmc.demhmcmc, buffer=buffer,corbeta=corbeta)


# Generate Final Chains

# In[12]:

buffer=np.copy(chain)
nbuffer=len(buffer[:,0])
corbeta=0.5 #scale for correlated jumps
parin=buffer[nbuffer-1,:] #we can start with the last state from our buffer.
niter=50000
chain,accept=mcmc.genchain(parin,data,betanew,niter,model,loglikelihood,lprior,mcmc.demhmcmc, buffer=buffer,corbeta=corbeta)
chain2,accept=mcmc.genchain(parin,data,betanew,niter,model,loglikelihood,lprior,mcmc.demhmcmc, buffer=buffer,corbeta=corbeta)
chain3,accept=mcmc.genchain(parin,data,betanew,niter,model,loglikelihood,lprior,mcmc.demhmcmc, buffer=buffer,corbeta=corbeta)


# ## Examine Chains

# In[13]:

mcmc.plotchains(chain,label,colour,burnin)


# ## Calculate Confidence Intervals

# In[14]:

npars=len(chain[1,:])
mm=np.zeros(npars)
for i in range(0,npars):
    perc = np.percentile(chain[burnin:,i],[0.3, 50.0, 99.7])
    mm[i]=perc[1]
    print('%s = %.3f +%.3f %.3f (3 Sigma)' %(label[i],perc[1],perc[2]-perc[1],perc[0]-perc[1]) )


# ## Overlay model chains with raw data

# In[15]:

mcmc.plotmodels(data,chain,model,burnin)


# ## Show cosine component ontop of ELC corrected lightcurve (median values adopted)

# In[16]:

plt.figure(figsize=(10,7)) 
tpars=np.copy(mm)
tpars[0]=0.0
mtest=model(tpars,data)
res=data[:,1]-mtest
plt.scatter(data[:,0],res,s=100.0,color='b',label='ELC Corrected')
plt.xlabel('Time (hours)')
plt.ylabel('Flux')

chainlen=len(chain[:,0])
for i in range(0,200): 
    nchain=int(np.random.rand()*(chainlen-burnin)+burnin) 
    tpars=np.copy(chain[nchain,:])
    tpars[3:7]=0
    mtest=model(tpars,data)
    plt.plot(data[:,0],mtest,color='r',alpha=0.1)

plt.plot(data[:,0],mtest,color='r',alpha=0.1,label='Cosine model')
plt.legend()
plt.show()


# ## Check convergence Tests

# In[17]:

npt=len(data[:,0]) #length of data set is need to get degrees of freedom
grtest=mcmc.gelmanrubin(chain,chain2,chain3,burnin=burnin,npt=npt)
print('Gelman-Rubin Convergence:')
print('parameter  Rc')
for i in range(0,len(chain[1,:])):
    print('%8s  %.3f' %(label[i],grtest[i]))


# Acceptance rates

# In[18]:

mcmc.calcacrate(accept,burnin,label)


# ## Triangle Plot

# In[19]:

nbin=30
mcmc.triplot(chain2,burnin,label,colour,nbin,ntick=4)


# In[20]:

pars1=np.copy(mm)
ll1=loglikelihood(model,pars1,data)
print(ll1)
pars2=np.copy(mm)
pars2[0]=0
ll2=loglikelihood(model,pars2,data)
print(ll2)
print(ll1-ll2)


def fun(x,data,amp):
    
    parsin=np.zeros(len(x)+1)
    parsin[0]=amp
    parsin[1:]=np.copy(x)
    
    f=model(parsin,data)-data[:,1]
    
    return f;
    
x=mm[1:]
amp=mm[0]
pfit, pcov, infodict, errmsg, success = leastsq(fun,x,args=(data,amp),full_output=1)
pars1=np.copy(mm); pars1[0]=amp; pars1[1:]=np.copy(pfit)
ll1=loglikelihood(model,pars1,data)
print(ll1)

x=mm[1:]
amp=0.006
pfit, pcov, infodict, errmsg, success = leastsq(fun,x,args=(data,amp),full_output=1)
pars2=np.copy(mm); pars2[0]=amp; pars2[1:]=np.copy(pfit)
ll2=loglikelihood(model,pars2,data)
print(ll2)

print(ll1-ll2)



