import crossmods
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import scipy.stats

dt = 1/30
dur = 20.0
ts = np.arange(0, dur, dt)
tau = 3.0 - ts
tau_b = tau - 2.0

params = dict(
    thm=np.log(2.0),
    ths=np.sqrt(1/6),
    lagm=np.log(1.0),
    lags=np.sqrt(1/6),
    pass_th=0.0,
        )


thdist = scipy.stats.lognorm(s=params['ths'], scale=np.exp(params['thm']))

plt.plot(ts, thdist.pdf(ts))
plt.show()
model = crossmods.LognormalTdm(**params)
#bpdf = model.blocker_decisions(tau, tau_b - np.inf, dt)
pdf = model.decisions(tau, dt)

plt.plot(ts, pdf.ps)
plt.twinx()
plt.plot(ts, tau)
print(pdf.uncrossed)
#plt.plot(ts, bpdf.ps)

plt.show()

"""
blocked = crossmods.LognormalTdm(**{**params, 'thm': np.inf}).decisions(tau_b, dt)
early = crossmods.LognormalTdm(**{**params, 'passm': -np.inf}).decisions(tau, dt)

unblocked = np.cumsum(np.array(blocked.ps)*dt)
early = np.cumsum(np.array(early.ps)*dt)

plt.plot(ts, unblocked)
plt.plot(ts, early)
plt.twinx()
plt.plot(ts, tau)
plt.plot(ts, tau_b)
plt.show()
"""

"""
wtf = crossmods.LognormalTdm(2.0, 0.3, -np.inf, 0.5, 0.0, 0.01).decisions(tau_b, dt)
#plt.plot(ts, wtf.ps)
plt.plot(np.cumsum(np.array(wtf.ps)*dt))
print(wtf.uncrossed)
plt.twinx()
plt.plot(ts, tau_b, color='black')
plt.show()
"""

#model = crossmods.LognormalTdm(np.inf, 0.3, 0.5, 0.5, 0.0, 0.1)
#pdf = model.decisions(tau_b, dt)
print(pdf.uncrossed)
sample = np.random.choice(ts, size=1000, p=np.array(pdf.ps)/(1 - pdf.uncrossed)*dt)
plt.hist(sample, bins=np.arange(ts[0], ts[-1], 0.2), density=True)
plt.plot(ts, np.array(pdf.ps))
#plt.plot(ts, np.cumsum(np.array(pdf.ps)*dt))
plt.twinx()
plt.plot(ts, tau, color='black')
plt.plot(ts, tau_b, '--', color='black')
plt.show()
niters = 100

start = time.perf_counter()
for i in range(niters):
    pdf = model.blocker_decisions(tau, tau_b, dt)
    loglik = pdf.loglikelihood(sample)

print(f"{(time.perf_counter() - start)/niters*1000:.2f} ms per call")
print(loglik)

#print(loglik)
#sample = np.random.choice(ts, size=100, p=np.array(pdf.ps)/(1 - pdf.uncrossed)*dt)


plt.show()
