#!/usr/bin/env python
# coding: utf-8

# In[209]:


import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# In[210]:


df = pd.read_csv("data_N1200_MCS400_STEP001_metropolis_boardsize37x37.csv",header=None)
df.head()


# In[211]:


len(df)


# In[212]:


exp = lambda t,tau:np.exp(-t/tau)
scaled_exp = lambda t,tau,a:a*np.exp(-t/tau)
Tc = 2.2691853


# In[231]:


tau = []
cutoff = 100
k=11
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
plt.title("Simulation Averages (batch 3 - metro)",fontsize=25)
plt.legend()
plt.tight_layout()
# plt.savefig("time_domain_simulation_metropolis_batch3.png")
plt.show()


# In[227]:


f = lambda temp, A, mu:A*(temp-Tc)**(-mu)
popt,pcov = curve_fit(f,temps,tau)
A,mu = popt
print("A={}, mu={}".format(round(A,3),round(mu,3)))
print("pcov",pcov)
perr = np.sqrt(np.diag(pcov))
print("error = ",perr)


# In[228]:


plt.figure(figsize=(7,5))
plt.plot(temps,tau,"x",label="time constant tau")
plt.plot([Tc]*3,np.linspace(min(tau),max(tau),3),"-",label="critical temperature")
temp_x = np.linspace(2.38,2.50,1000)
plt.plot(temp_x,f(temp_x,A,mu),"-",label="curve fit\nA = {}\nmu = {}".format(round(A,2),round(mu,3)))
plt.legend()
plt.title("Result, (batch 3 - metropolis)\nN=1200 MCS=400",fontsize=22)
plt.xlabel("temperature",fontsize=16)
plt.ylabel("time constant",fontsize=16)
plt.grid()
plt.tight_layout()
plt.savefig("result_plot_3.png")
plt.show()


# In[229]:


plt.figure(figsize=(7,5))
plt.plot(temps,tau,"x",label="time constant tau")
plt.plot([Tc]*3,np.linspace(min(tau),max(tau),3),"-",label="critical temperature")
temp_x = np.linspace(2.35,2.50,1000)
plt.plot(temp_x,f(temp_x,A,mu),"-",label="curve fit\nA = {}\nmu = {}".format(round(A,2),round(mu,3)))
plt.legend()
plt.title("Result, test data\nN=500 MCS=300",fontsize=22)
plt.xlabel("temperature",fontsize=16)
plt.ylabel("time constant",fontsize=16)
plt.grid()
plt.tight_layout()
# plt.savefig("result_plot_2xxxx.png")
plt.show()


# In[ ]:




