import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import argparse
import time

import read_input
import binning

parser = argparse.ArgumentParser(prog="hari: histogram arbitrary rebinning intelligent")
input_bins_group = parser.add_mutually_exclusive_group()
output_bins_group = parser.add_mutually_exclusive_group()
parser.add_argument("histogram", metavar=("HISTOGRAM_FILE"), help="Single-column histogram file")
input_bins_group.add_argument("-ib", "--input_bins", metavar=("BIN_FILE"), help="Single-column file that contains the bin centers of the input")
output_bins_group.add_argument("-ob", "--output_bins", metavar=("BIN_FILE"), help="Single-column file that contains the bin centers of the output")
input_bins_group.add_argument("-c", "--calibration", metavar=("CALIBRATION_FILE"), help="Use file with calibration parameters for polynomial")
parser.add_argument("-d", "--deterministic", help="Rebin deterministically", action="store_true")
output_bins_group.add_argument("-f", "--binning_factor", help="Rebinning factor")
parser.add_argument("-i", "--conserve_integral", help="Conserve sum over h[i]*dx[i] instead of sum over h[i]", action="store_true")
parser.add_argument("-l", "--limits", nargs=2, metavar=("LOWER_LIMIT", "UPPER_LIMIT"), help="Set limits of the plot range")
output_bins_group.add_argument("-n", "--n_bins", help="Number of bins for the output histogram")
parser.add_argument("-o", "--output", help="Set output file name")
parser.add_argument("-p", "--plot", help="Create plots of the rebinned spectrum", action="store_true")
parser.add_argument("-r", "--range", nargs=2, metavar=("START", "STOP"), help="Set range of the output bins")
parser.add_argument("-s", "--seed", help="Set random number seed")
parser.add_argument("-v", "--verbose", help="Print messages during program execution", action="store_true")
args=parser.parse_args()

start = time.time()

#
# Read input histogram
#

if args.verbose:
    print("> Reading input histogram from file", args.histogram)

input_hist = np.loadtxt(args.histogram)
input_hist_size = np.size(input_hist)

#
# Read input bins
#

if args.input_bins:

    if args.verbose:
        print("> Reading input bins from file", args.input_bins)

    input_bins = np.loadtxt(args.input_bins)
    n_input_bins = np.size(input_bins)
    if n_input_bins != input_hist_size:
        print("Error: Number of bins does not match.")
        print("\t Histogram size: ", input_hist_size)
        print("\t Number of bins: ", n_input_bins)
        exit()

elif args.calibration:
    if args.verbose:
        print("> Reading calibration parameters for input histogram from file", args.calibration)
    input_bins = read_input.calibrate(input_hist_size, args.calibration, read_input.remove_suffix_and_path(args.histogram), args.verbose)
    n_input_bins = input_hist_size

else:
    if args.verbose:
        print("> No input bins given, assume bin center == number of bin")
    input_bins = np.arange(0, input_hist_size)
    n_input_bins = input_hist_size

if args.verbose:
    print("> Input Histogram:", n_input_bins, "bins from", input_bins[0], "to", input_bins[-1])

#
# Create output bins
#

if args.range:
    output_bin_range = np.array([float(args.range[0]), float(args.range[1])])

    if output_bin_range[0] < input_bins[0] or output_bin_range[1] < input_bins[-1]:
        print("Warning: Range of output bins [", output_bin_range[0], ",", output_bin_range[1], "] is larger than input bins [", input_bins[0], ",", input_bins[-1], "]")

else:
    if args.verbose:
        print("> No range for output bins given, assume same range as input")
    output_bin_range = np.array([float(input_bins[0]), float(input_bins[-1])])

if args.binning_factor:
    if args.verbose:
        print("> Rebinning factor:", args.binning_factor)
    output_bins = np.linspace(output_bin_range[0], output_bin_range[1], int(n_input_bins/float(args.binning_factor)))

elif args.output_bins:
    if args.verbose:
        print("> Reading output bins from file", args.output_bins)
    output_bins = np.loadtxt(args.bins)
    if output_bins[0] < input_bins[0] or output_bins[-1] < input_bins[-1]:
        print("Warning: Range of output bins [", output_bin_range[0], ",", output_bin_range[1], "] is larger than input bins [", input_bins[0], ",", input_bins[-1], "]")

elif args.n_bins:
    if args.verbose:
        print("> Number of output bins set to:", args.n_bins)
    output_bins = np.linspace(output_bin_range[0], output_bin_range[1], int(args.n_bins))

else:
    if args.verbose:
        print("> No output bins given, assume bin center == number of bin")
    output_bins = np.linspace(input_bins[0], input_bins[-1], n_input_bins)

n_output_bins = np.size(output_bins)

if args.verbose:
    print("> Output Histogram:", n_output_bins, "bins from", output_bins[0], "to", output_bins[-1])

#
# Calculate the lower and upper limits of the bins
#

output_hist = np.zeros(n_output_bins)

input_bins_low, input_bins_high = binning.calculate_bin_limits(input_bins)
output_bins_low, output_bins_high = binning.calculate_bin_limits(output_bins)

#
# Interpolate the input histogram
#

# Calculate bins widths

input_bins_width = input_bins_high - input_bins_low
output_bins_width = output_bins_high - output_bins_low

if args.conserve_integral:
    if args.verbose:
        print("> Interpolating the input histogram. Conserving the integral over h[i]*dx[i].")
    inter = interpolate.InterpolatedUnivariateSpline(input_bins, input_hist)
else:
    if args.verbose:
        print("> Interpolating the input histogram. Conserving the sum over the bin contents h[i].")
    inter = interpolate.InterpolatedUnivariateSpline(np.arange(0., n_input_bins), input_hist)

inter_bins = interpolate.InterpolatedUnivariateSpline(input_bins, np.arange(0., n_input_bins))

#
# Calculate the bin contents of the output histogram
#

if args.conserve_integral:
    for i in range(0, n_output_bins):
        output_hist[i] = max(inter.integral(output_bins_low[i], output_bins_high[i]), 0.)/output_bins_width[i]

else:
    for i in range(0, n_output_bins):
        output_hist[i] = max(inter.integral(inter_bins(output_bins_low[i]), inter_bins(output_bins_high[i])), 0.)

if not args.deterministic:
    if args.verbose:
        print("> Rebinning and preserving the statistical fluctuations")
    if args.seed:
        np.random.seed(int(args.seed))
    output_hist = np.random.poisson(output_hist)
else:
    if args.verbose:
        print("> Rebinning without preserving the statistical fluctuations")

#
# Calculate calibration for output histogram
#

if args.binning_factor or args.n_bins:
    calibration_coefficients = np.array([output_bin_range[0], (output_bin_range[1] - output_bin_range[0])/n_output_bins]) 

#
# Write the result
#

if args.output:
    output_hist_filename = read_input.remove_suffix_and_path(args.output) + "_hist.txt"
    output_bins_filename = read_input.remove_suffix_and_path(args.output) + "_bins.txt"
else:
    output_hist_filename = read_input.remove_suffix_and_path(args.histogram) + "_hist.txt"
    output_bins_filename = read_input.remove_suffix_and_path(args.histogram) + "_bins.txt"

if args.verbose:
    print("> Writing output histogram to", output_hist_filename)
np.savetxt(output_hist_filename, output_hist, fmt='%.1f')
if args.verbose:
    print("> Writing output bins to", output_bins_filename)
np.savetxt(output_bins_filename, output_bins, fmt='%.6e')

if args.binning_factor or args.n_bins:
    if args.output:
        output_cal_filename = read_input.remove_suffix_and_path(args.output) + "_cal.txt"
    else:
        output_cal_filename = read_input.remove_suffix_and_path(args.histogram) + "_cal.txt"
    if args.verbose:
        print("> Writing calibration coefficients to", output_cal_filename)
    output_cal_file = open(output_cal_filename, "w")
    output_cal_file.write(output_hist_filename + ":")
    for c in calibration_coefficients:
        output_cal_file.write("\t")
        output_cal_file.write(str(c))

    output_cal_file.close()

stop = time.time()

#
# Print information
#

if args.verbose:
    print("> Execution took", stop-start, "seconds (without plotting)")

    if args.conserve_integral:
        # Calculate integral over histograms
        input_integral = np.sum(input_hist*input_bins_width)
        output_integral = np.sum(output_hist*output_bins_width)

        print("> Integral over input histogram:", input_integral)
        print("> Integral over output histogram:", output_integral)

        if args.deterministic:
            print("> Change of total histogram integral due to interpolation:", (input_integral - output_integral)/input_integral*100., "%")
        else:
            print("> Change of total histogram integral due to interpolation and random sampling:", (input_integral - output_integral)/input_integral*100. , "%")

    else:
        # Calculate sum of histogram bins
        n_input = np.sum(input_hist)
        n_output = np.sum(output_hist)

        print("> Sum of input histogram bins:", n_input)
        print("> Sum of output histogram bins:", n_output)

        if args.deterministic:
            print("> Change of total histogram content due to interpolation:", (n_output - n_input)/n_input*100., "%")
        else:
            print("> Change of total histogram content due to interpolation and random sampling:", (n_output - n_input)/n_input*100., "%")


#
# Plot the result
#

if args.plot:

    f, ax = plt.subplots(2, sharex=True)

    if args.limits:
        energy_range = np.array([float(args.limits[0]), float(args.limits[1])])
    else:
        energy_range = np.array([np.min(input_bins), np.max(input_bins)])

    ax[0].step(np.extract((input_bins > energy_range[0])*(input_bins < energy_range[1]), input_bins), np.extract((input_bins > energy_range[0])*(input_bins < energy_range[1]), input_hist), where="mid", color="black")
    ax[1].step(np.extract((output_bins > energy_range[0])*(output_bins < energy_range[1]), output_bins), np.extract((output_bins > energy_range[0])*(output_bins < energy_range[1]), output_hist), where="mid", color="green")

    plt.show()
