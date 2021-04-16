#!/usr/bin/env python
# coding: utf-8

# In[1]:


import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# In[2]:


MCS = 600
N = 2000
BS = "57x57"
sim_type = "metropolis"
df = pd.read_csv("data_N{}_MCS{}_STEP0002_{}_boardsize{}.csv".format(N,MCS,sim_type,BS),header=None)
df.head()


# In[3]:


len(df)


# In[4]:


exp = lambda t,tau:np.exp(-t/tau)
scaled_exp = lambda t,tau,a:a*np.exp(-t/tau)
Tc = 2.2691853


# In[5]:


tau = []
cutoff = 250
batch = 7
note = "_cut_high_T"
k=1
# k=10
temps = [i for i in df[0]]
temps = temps[:-k]

plt.figure(figsize=(7,5))
for i in range(len(df)-k):
    temp = df.iloc[i][0]
    y = np.array(df.iloc[i][cutoff+1:])
    x = np.arange(cutoff+1,cutoff+1+len(y))
    plt.plot(x,y,label="T/Tc = {}".format(round(temp/Tc,3)))
    popt,pcov = curve_fit(scaled_exp,x,y,p0=[850,0.6])
    tau.append(popt[0])
plt.xlabel("time in discrete steps",fontsize=16)
plt.ylabel("average magnetization",fontsize=16)
plt.title("Simulation Averages \nbatch {} metropolis".format(batch),fontsize=25)
# plt.legend()
plt.tight_layout()
# plt.savefig("time_domain_simulation_metropolis_batch{}{}.png".format(batch,cutoff,note))
plt.show()
print(len(df)-k)


# In[6]:


# get rid of the outlier
print(len(tau),len(temps))
tau_new = []
temps_new = []
for i,j in zip(tau,temps):
    # gets rid of outlier, high T, low T etc.
    if i>50 and j>2.30:
        tau_new.append(i)
        temps_new.append(j)
tau = tau_new
temps = temps_new
print(len(tau),len(temps))


# In[7]:


f = lambda temp, A, mu:A*(temp-Tc)**(-mu)
popt,pcov = curve_fit(f,temps,tau,p0=[7.7,1.82])
A,mu = popt
print("A={}, mu={}".format(round(A,3),round(mu,3)))
print("pcov",pcov)
perr = np.sqrt(np.diag(pcov))
print("error = ",perr)


# In[8]:


plt.figure(figsize=(7,5))
plt.plot(temps,tau,"x",label="time constant tau")
plt.plot([Tc]*3,np.linspace(min(tau),max(tau),3),"-",label="critical temperature")
temp_x = np.linspace(2.355,2.50,1000)
plt.plot(temp_x,f(temp_x,A,mu),"-",label="curve fit\nA = {} ± {} \nmu = {} ± {}".format(round(A,3),round(perr[0],2),round(mu,3),round(perr[1],3)))
plt.legend()
plt.title("Result, (batch {} - metropolis)\nN={} MCS={} BS={} cutoff={}".format(N,batch,MCS,BS,cutoff),fontsize=22)
plt.xlabel("temperature",fontsize=16)
plt.ylabel("time constant",fontsize=16)
plt.grid()
plt.tight_layout()
# plt.savefig("result_plot_{}_overnight_large_board_cutoff{}{}.png".format(batch,cutoff,note))
plt.show()

