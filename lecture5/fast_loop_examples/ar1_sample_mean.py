import numpy as np
import time
from numba import jit

@jit
def ar1_sample_mean(N, alpha, beta, s):
    x = beta / (1 - alpha)
    sm = 0.0
    for i in range(N):
        x = beta + alpha * x + s * np.random.randn()
        sm += x
    return sm / N

N = 10000000
alpha = 0.9
beta = 1.0
s = 1.0

t = time.time()
result = ar1_sample_mean(N, alpha, beta, s)
elapsed = time.time() - t

print("mean = {}".format(result))
print("elapsed time = {}".format(elapsed))
