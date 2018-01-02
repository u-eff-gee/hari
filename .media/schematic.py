import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

hist = np.array([0., 100., 130., 200., 200., 170., 0.])
bins = np.arange(0., 7.)

hist_simple = np.array([0., 0., 50., 50., 65, 65, 100., 100., 100., 100., 85, 85, 0., 0.])
bins_new = np.arange(-0.25, 6.75, 0.5)

#inter1 = interpolate.InterpolatedUnivariateSpline(bins, hist, k=1)
#inter1_x = np.arange(0., 7., 0.1)
#inter1_y = inter1(inter1_x)
inter3 = interpolate.InterpolatedUnivariateSpline(bins, hist, k=3)
inter3_x = np.arange(0., 7., 0.1)
inter3_y = inter3(inter3_x)

hist_inter = np.zeros(np.size(bins_new))
for i in range(np.size(bins_new) - 2):
    hist_inter[i] = inter3.integral(bins_new[i], bins_new[i+1])

fig, ax = plt.subplots(1)
ax.set_ylim(0., 250.)
ax.set_xlim(-1., 8.)

hist_fluc = np.random.poisson(hist_inter)

#
# Plot simple rebinning
#

#ax.fill_between([2.5, 3.5], 0., [200., 200.], facecolor="grey")
#
#ax.step(bins, hist, where="mid", color="black", label="Original histogram")
#ax.scatter(bins, hist, color="black")
#
#ax.step(bins_new, hist_simple, where="mid", color="red", label="New histogram")
#ax.scatter(bins_new, hist_simple, color="red")
#ax.legend()

#
# Plot interpolated rebinning
#

#fill_x = np.arange(2.5, 3.6, 0.1)
#fill_y = inter3(fill_x)
#ax.fill_between(fill_x, 0., fill_y, facecolor="grey")
#
#ax.step(bins, hist, where="mid", color="black", label="Original histogram")
#ax.scatter(bins, hist, color="black")
#
##ax.plot(inter1_x, inter1_y, "--", color="green", label="Spline interpolation (k = 1)")
#ax.plot(inter3_x, inter3_y, color="green", label="Spline (k = 3)")
#
#ax.step(bins_new, hist_inter, where="mid", color="red", label="New histogram")
#ax.scatter(bins_new, hist_inter, color="red")
#
#ax.legend()

#
# Plot interpolated rebinning with fluctuations
#

fill_x = np.arange(2.5, 3.6, 0.1)
fill_y = inter3(fill_x)
ax.fill_between(fill_x, 0., fill_y, facecolor="grey")

ax.step(bins, hist, where="mid", color="black", label="Original histogram")
ax.scatter(bins, hist, color="black")

#ax.plot(inter1_x, inter1_y, "--", color="green", label="Spline interpolation (k = 1)")
ax.plot(inter3_x, inter3_y, color="green", label="Spline (k = 3)")

ax.step(bins_new, hist_inter, where="mid", color="lightsalmon")

ax.step(bins_new, hist_fluc, where="mid", color="red", label="New histogram")
ax.scatter(bins_new, hist_fluc, color="red")

ax.legend()