import numpy as np

def calculate_bin_limits(bins):
    bins_low = bins - 0.5*(np.roll(bins, -1) - bins)
    bins_low[-1] = bins[-1] - 0.5*(bins[-1] - bins[-2])
    bins_high = bins + 0.5*(bins - np.roll(bins, 1))
    bins_high[0] = bins[0] + 0.5*(bins[1] - bins[0]) 
    
    return (bins_low, bins_high)

