"""
Sean Keenan, Heriot-Watt University, Edinburgh"
Mazzera Group Summer Project
Narrow Linewidth Laser Cavity Design
"""

# Import modules
from numpy import array, arange, empty, linspace, round
from math import exp, pi, sqrt, cos, log
from numpy.polynomial.polynomial import polyfit, polyval
import matplotlib.pyplot as mp

# set flag for which plots to generate, 0 for all, 1 for FSR and Linewidth, 2 for FP Spectrum
flag1 = 2
# set flag to save plots and path (0 no, 1 yes)
flag2 = 0
path = '/Users/Message/Documents/University/Projects/Laser Locking/'
# speed of light in vacuum (m/s)
c = 2.9972e8
# reflectivity of mirror 1
r_1 = 0.9985
# reflectivity of mirror 2
r_2 = 0.9985
# square of r
R = sqrt(r_1 * r_2)
# absorption co-efficient
alpha = 0
# radius of curvature (m)
C = 200e-3
# cavity length (m) - must be even number
L = linspace(start=15.0, stop=25.0, num=8) * 1e-2
# desired lambda (m)
wave = array([606e-9, 1550e-9])
# number of data points for model
num_steps = 10000
# array of wavelengths for model
wave_array = array([linspace(start=wave[0] - wave[0] / 1e5 , stop=wave[0] + wave[0] / 1e5, num=num_steps),
                    linspace(start=wave[1] - wave[1] / 1e5, stop=wave[1] + wave[1] / 1e5, num=num_steps)])
# subsequent frequency (Hz)
nu_0 = c / wave
# decay time constant
tau_c = -(2 * L) / (c * log(R ** 2 * (1 - alpha) ** 2))
# output linewidth (Hz)
delta_nu = (2 * pi * tau_c) ** -1
# free spectral range (Hz)
FSR = c / (2 * L)
# finesse from FSR
finesse_0 = FSR / delta_nu
# finesse from reflectivity
finesse = pi * sqrt(R) / (1 - R)
# intensity (arb. units)
I_0 = 100
# effective index of medium
n_eff = 1

# set global plot tick size
mp.rcParams['xtick.labelsize'] = 8
mp.rcParams['ytick.labelsize'] = 8

if flag1 == 0 or flag1 == 1:

    # create greater range of L
    L_new = linspace(L[0], L[-1], num=len(L) * 10)
    # find co-efficients of data
    coeff_1 = polyfit(x=L, y=delta_nu, deg=2)
    coeff_2 = polyfit(x=L, y=FSR, deg=2)
    # generate a fit for data
    fit_1 = polyval(L_new, coeff_1)
    fit_2 = polyval(L_new, coeff_2)
    # create subplots and handles to axis
    fig_0, axs_0 = mp.subplots(nrows=1, ncols=2, sharex='all', sharey='none', constrained_layout='true')
    # plot data to axis 1 and format
    axs_0[0].plot(L, delta_nu / 1e3, 'xr', markersize='8')
    axs_0[0].plot(L_new, fit_1 / 1e3, color='red', linestyle='dashed', label='fit')
    axs_0[0].set_title('Linewidth against Length', fontsize='10')
    axs_0[0].set(ylabel='Frequency (MHz)')
    # plot data to axis 2 and format
    axs_0[1].plot(L, FSR / 1e6, 'xb', markersize='8')
    axs_0[1].plot(L_new, fit_2 / 1e6, color='blue', linestyle='dashed', label='fit')
    axs_0[1].set_title('FSR against Length', fontsize='10')
    # set formatting for both axis
    for ax in axs_0.flat:
        # plots gridlines
        ax.grid()
        # display legend
        ax.legend()
        # set axis labels
        ax.set(xlabel='Length (m)')
        # change tick label format
        ax.ticklabel_format(style='scientific')

    # generate new figure and axis handle
    fig_1, axs_1 = mp.subplots()
    axs_1.grid()
    for index, value in enumerate(delta_nu):
        axs_1.plot(value / 1e3, FSR[index] / 1e6, 'x', markersize='8')

    axs_1.set_title('FSR against Linewidth', fontsize='12')
    axs_1.set(xlabel='Linewidth (MHz)', ylabel='Free Spectral Range (MHz)')
    axs_1.ticklabel_format(style='plain')
    axs_1.legend(L*100, fontsize=8, title='Cavity Length (cm)')

    if flag2 == 1:
        fig_0.savefig(fname=path + 'LwidthFSR_v_Length.pdf', quality=95, format='pdf')
        fig_1.savefig(fname=path + 'FSR_v_Lwidth.pdf', quality=95, format='pdf')

if flag1 == 0 or flag1 == 2:
    # create empty array for delta phase delay term
    delta = empty([len(wave_array), len(L), num_steps], dtype=float)

    # generate data for delta array
    for m in range(len(wave_array)):
        for n in range(len(L)):
            for o in range(len(wave_array[m])):
                delta[m][n][o] = (4 * pi * L[n]) / wave_array[m][o]

    # create empty arrays for model
    denominator = empty([len(L), num_steps], dtype=float)
    I_t = empty([len(wave_array), len(L), num_steps], dtype=float)

    # generate data for model
    for p in range(len(wave_array)):
        for q in range(len(L)):
            numerator = I_0 * (1 - R * exp(-2 * alpha * L[q])) ** 2
            for r in range(len(wave_array[p])):
                denominator = 1 + R ** 2 * exp(-4 * alpha * L[q]) - 2 * R * exp(-2 * alpha * L[q]) * cos(delta[p][q][r])
                I_t[p][q][r] = numerator / denominator

    # fig_2, axs_2 = mp.subplots(nrows=2, ncols=len(L), sharex='none', sharey='row', constrained_layout='true')
    # fig_2.suptitle('Output Intensity of Cavity of Length L', fontsize='12')

    plot_col = ['red', 'blue']
    for s in range(len(wave_array)):
        fig_2, axs_2 = mp.subplots(nrows=len(wave), ncols=int(len(L)/2), sharex='all', sharey='all', constrained_layout=False)
        fig_2.suptitle('Output Intensity of Cavity of Length L', fontsize='12')
        for index, ax in enumerate(axs_2.flatten()): #range(len(L)):
            ax.grid()
            ax.plot(wave_array[s] * 1e9, I_t[s][index][:], color=plot_col[s])
            ax.set(ylabel='Intensity (arb. units)')
            ax.set_title('L =' + str(round(L[index] * 100, 2)) + 'cm', fontsize='10')
            ax.autoscale(enable=True, axis='both', tight=True)
            ax.ticklabel_format(useOffset=True, style='plain')
            ax.label_outer()
        fig_2.text(0.5, 0.01, 'Wavelength (nm)', ha='center')
        if flag2 == 1:
            fig_2.savefig(fname=path + str(wave[s]) + 'FP_output.pdf', quality=95, format='pdf')
