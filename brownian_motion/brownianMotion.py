import numpy as np
import math
import matplotlib.pyplot as plt
import random

245d19d2797c38a58aed7cff76f92fec65eb4d13


sin = math.sin
cos = math.cos
arctan = math.atan

# global variables
N = 1500
PI = math.pi
MI = 1
MO = 250
R = 0.03
MU,SIGMA = 0.02,0.01 # for the velocities of particles
STEPS = 2000

# quick helper functions
norm = lambda x:  math.sqrt(sum([i**2 for i in x]))
unitvec = lambda x:  x / norm(x)
polarToCart = lambda r,theta: np.array([math.cos(theta)])
rotateCartVec = lambda vec: np.array([-vec[1],vec[0]])
polarToCart = lambda vec2: np.array([vec2[0]*cos(vec2[1]),vec2[0]*sin(vec2[1])])

def init():# returns the initial values of q,v,pos,vel, uses global variables
    q = np.array([np.random.rand(N),np.random.rand(N)]).T
    polarV = np.array([np.random.normal(MU,SIGMA,N),np.random.random(N)*PI*2]).T
    v = np.array([polarToCart(j) for j in polarV])
    pos,vel = np.array([0.5,0.5]),np.array([0.0,0.0])
    return q,v,pos,vel

def update(q,v,pos,vel):# pos,vel is pollen mollecule, others are little particles
    qNew = q + v
    vNew = v + np.zeros(v.shape) # deep copy
    deltaP = np.array([0.0,0.0])
    nCollisions = 0 
    for idx,(qi,qiNew,vi) in enumerate(zip(q,qNew,v)):
        if norm(qi - pos%1) < R: # the modular thing is because - for plotting purposes, i don't display trajectory
            # of the pollen mollecule as being on a torus, even though it is
            nCollisions += 1
            
            # find incoming angle (phi) [angle in ball], and angle of collision (theta)
            u = unitvec(qi-qiNew)
            if u[0]>0:
                theta = arctan(u[1]/u[0])
            else:
                theta = arctan(u[1]/u[0])+PI
            ballNormal = unitvec(pos-qiNew)
            if ballNormal[0]>0:
                phi = arctan(ballNormal[1]/ballNormal[0])
            else:
                phi = arctan(ballNormal[1]/ballNormal[0])
                
            # find the velocity of incomming particle in the ball frame
            viBallFrame = vi - vel
            
            # approximation M >> m for collision, bounces back with same velocity
            ballNormal = unitvec(pos-qiNew)
            ballTanj = rotateCartVec(ballNormal)
            viNewBallFrame = ballTanj * np.dot(viBallFrame,ballTanj) - ballNormal * np.dot(viBallFrame,ballNormal)
            
            deltaP += (viBallFrame - viNewBallFrame)*MI # this is the imparted momentum
            viNew = viNewBallFrame + vel # bring it back to the lab frame
            
            # update arrays
            vNew[idx] = viNew
            qNew[idx] = q[idx] + 2*R*unitvec(viNew) # + vel + 2*ballTanj * np.dot(viBallFrame,ballTanj) # not super sure about this
            # factor chosen so that the particle to escape
            
        # lets put everyone on a torus
        if qiNew[0] > 1:
            qNew[idx][0] = qiNew[0] - 1
        elif qiNew[0] < 0:
            qNew[idx][0] = qiNew[0] + 1
        if qiNew[1] > 1:
            qNew[idx][1] = qiNew[1] - 1
        elif qiNew[1] < 0:
            qNew[idx][1] = qiNew[1] + 1
            
        # update position
        velNew = vel + deltaP / MO
        posNew = pos + velNew
        
    return qNew,vNew,posNew,velNew,nCollisions
   
   
if __name__=="__main__":
	q,v,pos,vel = init()

	# compute the energy of the system
	energy = round(sum([0.5*norm(vi)**2 for vi in v]),3)

	pollenPath = [pos]
	pollenVel = [vel]
	collisions = []

	for i in range(k):
	    q,v,pos,vel,nCollisions = update(q,v,pos,vel)
	    pollenPath.append(pos)
	    pollenVel.append(vel)
	    collisions.append(nCollisions)
	    
        if i%(STEPS//5)==0:
            plt.figure()
            plt.subplots(figsize=(14,7))
			plt.subplot(1,2,1)
			plt.plot(q.T[0],q.T[1],"b.")
			plt.title("pollen (red) in a sea of particles")
			plt.plot(pos[0]%1,pos[1]%1,"ro",label="N={} E={} tStep={}".format(N,energy,i))
			plt.legend()
			plt.subplot(1,2,2)
			plt.title("the velocities of the particles (random gaussian)")
			plt.plot(v.T[0],v.T[1],".")
			plt.show()
		
		if i==0:
			plt.subplots(figsize=(14,7))
			plt.subplot(1,2,1)
			plt.plot(q.T[0],q.T[1],"b.")
			plt.title("pollen (red) in a sea of particles")
			plt.plot(pos[0]%1,pos[1]%1,"ro",label="N={} E={} tStep={}".format(N,energy,i))
			plt.legend()
			plt.subplot(1,2,2)
			plt.title("the velocities of the particles (random gaussian)")
			plt.plot(v.T[0],v.T[1],".")
			plt.savefig("pollensimulation.png")
			plt.show()

	# plot the pollen path and it's velocity trace
	pollenPath = np.array(pollenPath)
	pollenVel = np.array(pollenVel)
		
	plt.subplots(figsize=(14,7))
	plt.subplot(1,2,1)
	plt.plot(pollenPath.T[0],pollenPath.T[1],"r-")
	plt.title("pollen path")
	plt.subplot(1,2,2)
	plt.plot(pollenVel.T[0],pollenVel.T[1],"k-")
	plt.title("pollen veolcity path (should be tethered)")
	plt.savefig("pollen_path.png")
	plt.show()

	# plot dist of collisions (should be poissonian)
	plt.figure(figsize=(8,8))
	plt.title("distribution of number of collisions per timestep")
	# plt.hist(collisions)#,rwidth=0.9)
	plt.hist(collisions, bins=np.arange(max(collisions)+1)-0.5, ec="k",label="N={} E={} STEPS={} MO/MI={} R={}".format(N,energy,STEPS,int(MO/MI),R))
	plt.legend()
	plt.savefig("collisions_per_timestep_distribution.png")
	plt.show()

