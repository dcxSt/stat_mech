#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# In[2]:


# load data
MCS = 550
N = 1200
BS = "37x37"
sim_type = "glauber"
STEP = "0003"
name = "data_N{}_MCS{}_STEP{}_{}_boardsize{}_all".format(N,MCS,STEP,sim_type,BS)
print(name)
df = pd.read_csv("{}.csv".format(name),header=None)
# df = pd.read_csv("data_glauber2.csv",header=None)
df.head()


# In[3]:


len(df) # number of data points


# In[6]:


# functions for finding tau
scaled_exp = lambda t,tau,a:a*np.exp(-t/tau)
Tc = 2.2691853


# For each temperature value, find the time constant, add it to the list. We find the time constant &tau; by optimizing the parameters &tau; and c which fit the tail of the average magnetization at each time of each simulation of the ising model for that temperature. (tail corresponds to long time behaviour)

# In[21]:


# for each temperature value find the time constant, add it to list

tau = []
cutoff = 200
batch = 2
note = "_cut_high_T"
# k=10
temps = [i for i in df[0]]

plt.figure(figsize=(7,5))
for i in range(len(df)):
    temp = df.iloc[i][0]
    y = np.array(df.iloc[i][cutoff+1:])
    x = np.arange(cutoff+1,cutoff+1+len(y))
    plt.plot(x,y,label="T/Tc = {}".format(round(temp/Tc,3)))
    popt,pcov = curve_fit(scaled_exp,x,y,p0=[850,0.6])
    tau.append(popt[0])
plt.xlabel("time in discrete steps",fontsize=16)
plt.ylabel("average magnetization",fontsize=16)
plt.title("Simulation Averages \nbatch {} {}".format(batch,sim_type),fontsize=25)
# plt.legend()
plt.tight_layout()
plt.savefig("time_domain_simulation_glauber_batch{}{}.png".format(batch,cutoff,note))
plt.show()
print(len(df))


# In[22]:


# get rid of the outlier
print(len(tau),len(temps))
tau_new = []
temps_new = []
for i,j in zip(tau,temps):
    # gets rid of outlier, high T, low T etc.
    if i>50:# and j<2.5:
        tau_new.append(i)
        temps_new.append(j)
tau = tau_new
temps = temps_new
print(len(tau),len(temps))


# Find best fit parameters for A and &mu;

# In[23]:


f = lambda temp, A, mu:A*(temp-Tc)**(-mu)
popt,pcov = curve_fit(f,temps,tau,p0=[7.7,1.82])
A,mu = popt
print("A={}, mu={}".format(round(A,3),round(mu,3)))
print("pcov",pcov)
perr = np.sqrt(np.diag(pcov))
print("error = ",perr)


# In[26]:


plt.figure(figsize=(7,5))
plt.plot(temps,tau,"x",label="time constant tau")
plt.plot([Tc]*3,np.linspace(min(tau),max(tau),3),"-",label="critical temperature")
temp_x = np.linspace(2.355,2.53,1000)
plt.plot(temp_x,f(temp_x,A,mu),"-",label="curve fit\nA = {} ± {} \nmu = {} ± {}".format(round(A,3),round(perr[0],2),round(mu,3),round(perr[1],3)))
plt.legend()
plt.title("Result, (batch {} - {})\nN={} MCS={} BS={} cutoff={}".format(batch,sim_type,N,MCS,BS,cutoff),fontsize=22)
plt.xlabel("temperature",fontsize=16)
plt.ylabel("time constant",fontsize=16)
plt.grid()
plt.tight_layout()
plt.savefig("result_plot_{}_cutoff{}{}.png".format(name,cutoff,note))
plt.show()

