import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gauss(x, b, c, mu, sigma):
    return b + c/np.sqrt(2.*np.pi*sigma*sigma)*np.exp(-(x - mu)*(x - mu)/(2.*sigma*sigma))

DIR = ""

#
# Create output
#
#bins = np.linspace(0., 6., 1000)
#hist = gauss(bins, 20., 100., 3., 0.3)
#hist = np.random.poisson(hist)
#
#plt.step(bins, hist)
#
#np.savetxt(DIR + "gauss_bins.txt", bins)
#np.savetxt(DIR + "gauss_hist.txt", hist)

#
# Fit output
#
#bins = np.loadtxt(DIR + "gauss_bins.txt")
#hist = np.loadtxt(DIR + "gauss_hist.txt")
#
#popt, pcov = curve_fit(gauss, bins, hist, p0=[1., 100., 1., 1.])
#perr = np.sqrt(np.diag(pcov))
#
#name = "Original spectrum"
#print(name)
#print("Fit results:")
#print("B     =", popt[0], "+-", perr[0], "(", perr[0]/popt[0]*100., "% )")
#print("C     =", popt[1], "+-", perr[1], "(", perr[1]/popt[1]*100., "% )")
#print("MU    =", popt[2], "+-", perr[2])
#print("SIGMA =", popt[3], "+-", perr[3], "(", perr[3]/popt[3]*100., "% )")
#print()
#
#fig, ax = plt.subplots(1)
#ax.set_ylim(0., 200.)
#ax.step(bins, hist, color="black", label="Spectrum", where="mid")
#ax.plot(bins, gauss(bins, popt[0], popt[1], popt[2], popt[3]), color="red", label="Fit")
#ax.legend(fontsize=12)
#text = "b = {:3.3f} ($\pm$ {:3.2f} %)\nc = {:3.3f} ($\pm$ {:3.2f} %)\n$\mu$ = {:3.3f} ($\pm$ {:3.2f} %)\n$\sigma$ = {:3.3f} ($\pm$ {:3.2f} %)".format(popt[0], perr[0]/popt[0]*100., popt[1], perr[1]/popt[1]*100., popt[2], perr[2]/popt[2]*100., popt[3], perr[3]/popt[3]*100.)
#ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=12,
#        verticalalignment='top')

####################################################################

#bins = np.loadtxt(DIR + "gauss_f20_bins.txt")
#hist = np.loadtxt(DIR + "gauss_f20_hist.txt")
#
#popt, pcov = curve_fit(gauss, bins, hist, p0=[1., 2000., 1., 1.])
#perr = np.sqrt(np.diag(pcov))
#
#name = "Factor 20"
#print(name)
#print("Fit results:")
#print("B     =", popt[0], "+-", perr[0], "(", perr[0]/popt[0]*100., "% )")
#print("C     =", popt[1], "+-", perr[1], "(", perr[1]/popt[1]*100., "% )")
#print("MU    =", popt[2], "+-", perr[2])
#print("SIGMA =", popt[3], "+-", perr[3], "(", perr[3]/popt[3]*100., "% )")
#print()
#
#fig, ax = plt.subplots(1)
#ax.set_ylim(0., 3500.)
#ax.step(bins, hist, color="black", label="Spectrum", where="mid")
#ax.plot(bins, gauss(bins, popt[0], popt[1], popt[2], popt[3]), color="red", label="Fit")
#ax.legend(fontsize=12)
#text = "b = {:3.3f} ($\pm$ {:3.2f} %)\nc = {:3.3f} ($\pm$ {:3.2f} %)\n$\mu$ = {:3.3f} ($\pm$ {:3.2f} %)\n$\sigma$ = {:3.3f} ($\pm$ {:3.2f} %)".format(popt[0], perr[0]/popt[0]*100., popt[1], perr[1]/popt[1]*100., popt[2], perr[2]/popt[2]*100., popt[3], perr[3]/popt[3]*100.)
#ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=12,
#        verticalalignment='top')
#
######################################################################
#
#bins = np.loadtxt(DIR + "gauss_f20r_bins.txt")
#hist = np.loadtxt(DIR + "gauss_f20r_hist.txt")
#
#popt, pcov = curve_fit(gauss, bins, hist, p0=[1., 100., 1., 1.])
#perr = np.sqrt(np.diag(pcov))
#
#name = "Factor 20 * 1/20"
#print(name)
#print("Fit results:")
#print("B     =", popt[0], "+-", perr[0], "(", perr[0]/popt[0]*100., "% )")
#print("C     =", popt[1], "+-", perr[1], "(", perr[1]/popt[1]*100., "% )")
#print("MU    =", popt[2], "+-", perr[2])
#print("SIGMA =", popt[3], "+-", perr[3], "(", perr[3]/popt[3]*100., "% )")
#print()
#
#fig, ax = plt.subplots(1)
#ax.set_ylim(0., 200.)
#ax.step(bins, hist, color="black", label="Spectrum", where="mid")
#ax.plot(bins, gauss(bins, popt[0], popt[1], popt[2], popt[3]), color="red", label="Fit")
#ax.legend(fontsize=12)
#text = "b = {:3.3f} ($\pm$ {:3.2f} %)\nc = {:3.3f} ($\pm$ {:3.2f} %)\n$\mu$ = {:3.3f} ($\pm$ {:3.2f} %)\n$\sigma$ = {:3.3f} ($\pm$ {:3.2f} %)".format(popt[0], perr[0]/popt[0]*100., popt[1], perr[1]/popt[1]*100., popt[2], perr[2]/popt[2]*100., popt[3], perr[3]/popt[3]*100.)
#ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=12,
#        verticalalignment='top')
#
######################################################################
#
bins = np.loadtxt(DIR + "gauss_f20r_d_bins.txt")
hist = np.loadtxt(DIR + "gauss_f20r_d_hist.txt")

popt, pcov = curve_fit(gauss, bins, hist, p0=[1., 100., 1., 1.])
perr = np.sqrt(np.diag(pcov))

name = "Factor 20 * 1/20 deterministic"
print(name)
print("Fit results:")
print("B     =", popt[0], "+-", perr[0], "(", perr[0]/popt[0]*100., "% )")
print("C     =", popt[1], "+-", perr[1], "(", perr[1]/popt[1]*100., "% )")
print("MU    =", popt[2], "+-", perr[2])
print("SIGMA =", popt[3], "+-", perr[3], "(", perr[3]/popt[3]*100., "% )")
print()

fig, ax = plt.subplots(1)
ax.set_ylim(0., 200.)
ax.step(bins, hist, color="black", label="Spectrum", where="mid")
ax.plot(bins, gauss(bins, popt[0], popt[1], popt[2], popt[3]), color="red", label="Fit")
ax.legend(fontsize=12)
text = "b = {:3.3f} ($\pm$ {:3.2f} %)\nc = {:3.3f} ($\pm$ {:3.2f} %)\n$\mu$ = {:3.3f} ($\pm$ {:3.2f} %)\n$\sigma$ = {:3.3f} ($\pm$ {:3.2f} %)".format(popt[0], perr[0]/popt[0]*100., popt[1], perr[1]/popt[1]*100., popt[2], perr[2]/popt[2]*100., popt[3], perr[3]/popt[3]*100.)
ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=12,
        verticalalignment='top')
#
######################################################################
#
#bins = np.loadtxt(DIR + "gauss_f20_d_bins.txt")
#hist = np.loadtxt(DIR + "gauss_f20_d_hist.txt")
#
#popt, pcov = curve_fit(gauss, bins, hist, p0=[1., 2000., 1., 1.])
#perr = np.sqrt(np.diag(pcov))
#
#name = "Factor 20 deterministic"
#print(name)
#print("Fit results:")
#print("B     =", popt[0], "+-", perr[0], "(", perr[0]/popt[0]*100., "% )")
#print("C     =", popt[1], "+-", perr[1], "(", perr[1]/popt[1]*100., "% )")
#print("MU    =", popt[2], "+-", perr[2])
#print("SIGMA =", popt[3], "+-", perr[3], "(", perr[3]/popt[3]*100., "% )")
#print()
#
#fig, ax = plt.subplots(1)
#ax.set_ylim(0., 3500.)
#ax.step(bins, hist, color="black", label="Spectrum", where="mid")
#ax.plot(bins, gauss(bins, popt[0], popt[1], popt[2], popt[3]), color="red", label="Fit")
#ax.legend(fontsize=12)
#text = "b = {:3.3f} ($\pm$ {:3.2f} %)\nc = {:3.3f} ($\pm$ {:3.2f} %)\n$\mu$ = {:3.3f} ($\pm$ {:3.2f} %)\n$\sigma$ = {:3.3f} ($\pm$ {:3.2f} %)".format(popt[0], perr[0]/popt[0]*100., popt[1], perr[1]/popt[1]*100., popt[2], perr[2]/popt[2]*100., popt[3], perr[3]/popt[3]*100.)
#ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=12,
#        verticalalignment='top')