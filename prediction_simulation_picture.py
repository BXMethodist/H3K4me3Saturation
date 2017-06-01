from simulation_generator import get_random_peak
from predict import optimize_allocs
from visualizationUtils import plot_predict
import numpy as np
import pandas as pd

def get_simulated_variants():
    variants = []

    variant = np.zeros(1000)
    variant[100:200] = get_random_peak(100, 100)
    variant[800:900] = get_random_peak(100, 100)
    variants.append(variant)

    variant = np.zeros(1000)
    variant[400:500] = get_random_peak(100, 100)
    # variant[700:900] = get_random_peak(100, 100)
    variants.append(variant)

    return np.asarray(variants)


variants = get_simulated_variants()

data = np.zeros(1000)
data[600:700] = get_random_peak(100, 100)
# data = variants[0]*0.8 + variants[1] * 0.2
allocs = optimize_allocs(data, variants)

print allocs

predicted = 0
for j in range(len(variants)):
    signal = (np.sum(data)/np.sum(variants[j])) * variants[j]
    # print np.sum(signal), np.sum(data)
    predicted += signal * allocs[j]

error = np.absolute(data-predicted).sum()/np.sum(predicted)
corr = np.corrcoef(predicted, data)[0,1]
print corr
plot_predict(data, variants, allocs)