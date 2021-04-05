import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def rw_mean(k,t):
	# k is the number of walks, t is the duration of each walk
	x = np.zeros(k)
	for i in range(t):
		x+=np.random.choice((-1,2),k)
	return x.mean()

def rw_delta(k,t):
	x = np.zeros(k)
	for i in range(t):
		x+=np.random.choice((-1,2),k)
	xmean = np.ones(k)*t/2
	dx = x-xmean
	return np.dot(dx,dx)/len(dx)

def linear(x,a):
	return a*x

### find variance emperically
if __name__ == "__main__":
	var = []
	time = [t for t in range(1000,10000,1000)]
	for t in time:
		var.append(rw_delta(1000,t))
	popt = curve_fit(linear,time,var)
	plt.plot(time,[linear(t,popt[0]) for t in time],label='popt={}'.format(popt[0]))
	plt.plot(time,var,'x',label='data')
	plt.title("variance x(t)")
	plt.legend()
	plt.show()


### find mean empirically
if __name__ == "__main__":
	means = []
	times = []
	for t in range(1000,10000,1000):
		times.append(t)
		means.append(rw_mean(1000,t))
	popt,pcov = curve_fit(linear,times,means)
	plt.plot(times,[linear(t,popt[0]) for t in times],label='curve_fit popt={}'.format(popt[0]))
	plt.plot(times,means,'x',label='data')
	plt.legend()
	plt.title('random walk simulation')
	plt.show()
