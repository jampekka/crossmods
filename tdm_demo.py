import crossmods
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dt = 1/30
dur = 20.0
ts = np.arange(0, dur, dt)
tau = 5.0 - ts

model = crossmods.LognormalTdm(2.0, 0.3, 0.5, 0.5, 0.0, 0.1)
pdf = model.decisions(tau, dt)
print(pdf.uncrossed)
plt.plot(ts, pdf.ps)
plt.twinx()
plt.plot(ts, tau, color='black')
plt.show()
