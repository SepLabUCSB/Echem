import numpy as np
import matplotlib.pyplot as plt


y = np.array([1,1,1.1,1,0.9,1,1,1.1,1,0.9,1,1.1,1,1,0.9,1,1,1.1,1,1,1,1,1.1,0.9,1,1.1,1,1,0.9,
       1,1.1,1,1,1.1,1,0.8,0.9,1,1.2,0.9,1,1,1.1,1.2,1,1.5,1,3,2,5,3,2,1,1,1,0.9,1,1,3,
       2.6,4,3,3.2,2,1,1,0.8,4,4,2,2.5,1,1,1])

x = np.arange(0, len(y), 1)



def thresholding_algo(y, lag, threshold, influence):
    peaks = []
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:lag])
    stdFilter[lag - 1] = np.std(y[0:lag])
    for i in range(lag, len(y)):
        if abs(y[i] - avgFilter[i-1]) > threshold * stdFilter [i-1]:
            peaks.append(i)

            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])
        else:
            signals[i] = 0
            filteredY[i] = y[i]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])

    avgFilter = np.asarray(avgFilter)
    stdFilter = np.asarray(stdFilter)

    return peaks, avgFilter, stdFilter


peaks, avgFilter, stdFilter = thresholding_algo(y, 30, 5, 0)

for i in peaks:
    m = 10
    n = np.where(y == max(y[i-m:i+m+1]))[0][0]
    if i != n:
        peaks.remove(i)

        

        # for (x,y) in list(self.criticalPoints):  
        #     # find largest local step, search +- m points
        #     m = 10
        #     xi = np.where(self.xdata == x)[0][0] #convert to index
        #     n = np.where(delta == max(delta[xi-m:xi+m+1]))[0][0]
        #     if not n == xi:
        #         # for m in self.criticalPoints[(x,y)]:
        #         #     m.remove()
        #         self.criticalPoints.pop((x,y), None)
        #         self.drawPoints(self.ax, self.xdata[n], self.ydata[n])

#%%
plt.plot(x,y)
plt.ylim(0.5, )
