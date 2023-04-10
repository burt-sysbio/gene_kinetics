#%%
from utils import gamma_mixture

import numpy as np
x = np.arange(0,10,0.1)
y = gamma_mixture(x, alpha1 = 1, beta1 = 0.1, alpha2 = 10, beta2= 10)

import matplotlib.pyplot as plt

plt.plot(x,y)


# %%
