# python code
import numpy as np
import random
import math
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

sqrt = math.sqrt
# generate such a random walk with 8n steps
rw = lambda n:np.multiply(np.random.choice([-1,1],8*n),np.tile(np.array([1,1,1,1,1,1,1,3]),n))
def getWalk(n):
    steps = rw(n)
    walk = np.zeros(steps.shape)
    for i in range(len(steps)):
        walk += np.concatenate([np.zeros(i),steps])[:len(steps)]
    rmsTrace = [np.linalg.norm(walk[:i+1])/sqrt(i+1) for i in range(len(walk))]
    rmsTrace = np.array(rmsTrace)
    return steps,walk,rmsTrace
    
print("exit the figure for program to continue, for statistically reliable result, wait like 5 mins for a bunch of simulations to take place after that!")

# plot some figures
plt.figure(figsize=(12,6))
for i in range(3):
    steps,walk,rmsTrace=getWalk(2000)
    plt.subplot(2,3,i+1)
    plt.plot(walk)
    plt.title("Walk {}".format(i+1))
    plt.subplot(2,3,i+4)
    plt.plot(rmsTrace)
    plt.title("ROOT MEAN SQUARE {}".format(i+1))
plt.tight_layout()
plt.savefig("randomWalk1DPlotsRMS.png")
plt.show()

# simulate 100 rw for n=16000, put them in matrix and take averate of rms
import datetime as dt
import os
n = 2000
k = 200
rmsMatrix = np.zeros((k,8*n))
for i in range(k):
    if i%5==0:print(i,end="\t")
    steps,walk,rmsTrace=getWalk(n)
    rmsMatrix[i] = rmsTrace
        

rmsMean = np.mean(rmsBigMatrix,axis=0)
np.save("rmsMean_{}.npy".format(str(dt.datetime.now())),rmsMean)

        
def inClassRMS(t):
    return sqrt(t)

def angelaFit(t,eps):
    return np.sqrt((1+2*eps))*np.sqrt(t)

def donaldFit(t,delta):
    return t**(delta+0.5)

# we're concerned about for large n, so let's get rid of first 200 entries of rmsMean
rmsMean = np.load("rmsBigMean_2021-01-31 07:26:34.424188.npy")
t = np.arange(len(rmsMean))
rmsMean = rmsMean[200:]
t = t[200:]
eps, pcovA = curve_fit(angelaFit,t,rmsMean)
delta, pcovD = curve_fit(donaldFit,t,rmsMean)

plt.figure(figsize=(8,8))
plt.plot(t,rmsMean,label="rms mean",alpha=0.6)
plt.plot(t,angelaFit(t,eps),label="Angela Fit eps={}\nError={}".format(round(eps[0],3),pcovA[0][0]),alpha=0.7)
plt.plot(t,donaldFit(t,delta),label="Donald Fit delta={}\nError={}".format(round(delta[0],3),pcovD[0][0]),alpha=0.7)
plt.legend()
plt.title("scipy curve fits")
plt.xlabel("steps t")
plt.ylabel("rms")
plt.savefig("Pset1Q3AngelaDonaldCurveFits.png")
