import crossmods
import numpy as np
import matplotlib.pyplot as plt
import time

def main():
    dt = 1/30.0
    vddm = crossmods.Vddm(
            dt=dt, std=0.75,
            damping=1.6,
            tau_threshold=2.3,
            pass_threshold=0.0,
            )
    
    dur = 20.0
    tau0 = 10.0

    dur = tau0
    ts = np.arange(0, dur, dt)
    tau = tau0 - ts
    
    tau_block = tau - 3.0

    grid = crossmods.Grid1d(-3.0, 3.0, 100)
    
    pdf = vddm.blocker_decisions(grid, tau, tau_block)
    #pdf = vddm.decisions(grid, tau)
    samplesize = 1000
    sample = np.random.choice(ts, size=samplesize, p=np.array(pdf.ps)/np.sum(pdf.ps))
    sample[np.random.rand(samplesize) < pdf.uncrossed] = np.inf
    
    #plt.hist(sample[np.isfinite(sample)], bins=np.arange(ts[0], ts[-1], 0.2), density=True)
    #plt.plot(ts, np.array(pdf.ps)/(1 - pdf.uncrossed))
    #plt.show()

    #print(list(map(pdf, sample)))
    niter = 100
    

    start = time.perf_counter()
    for i in range(niter):
        #loglik = vddm.loglikelihood(grid, tau, sample)
        #loglik = vddm.decisions(grid, tau).loglikelihood(sample)
        loglik = vddm.blocker_decisions(grid, tau, tau_block).loglikelihood(sample)
    
    print((time.perf_counter() - start)/niter*1000)
    plt.plot(ts, pdf.ps)
    plt.twinx()
    plt.plot(ts, tau)
    plt.axvline(ts[np.argmin(np.abs(tau_block))], color='black')
    plt.show()


if __name__ == '__main__': main()
