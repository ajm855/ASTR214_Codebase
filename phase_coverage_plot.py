# Phase Coverage Graphing Tool
# Alexander Magnus
import os
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt


def phase(jd, p):
    """
    :param jd: Float. Julian date to be phased.
    :param p: Float. Period to phase by.
    :return: Float. Phase of observation.
    """
    return (jd % p) / p


def add_prediction(plt_ax, p, jd_start, jd_end, colour='r'):
    """
    :param plt_ax: Pyplot axis instance to add the phase range prediction to.
    :param p: Float. Period to phase by.
    :param jd_start: Float. Starting time of Julian date range.
    :param jd_end: Float. Ending time of Julian date range.
    :param colour: String. Color of the prediction line from pyplot colour naming convention.
    :return: Null.
    """
    # Automatically obtain UTC MMDD from JD input (at start of range).
    t = Time(jd_start, format='jd', scale='utc')
    utc_split = t.iso.split('-')
    date_label = utc_split[1] + utc_split[2].split(' ')[0]

    if phase(jd_start, p) < phase(jd_end, p):
        plt_ax.plot([phase(jd_start, p), phase(jd_end, p)], [date_label, date_label], colour)
    else:
        # Handling for phase wrap-around.
        plt_ax.plot([phase(jd_start, p), 1], [date_label, date_label], colour)
        plt_ax.plot([0, phase(jd_end, p)], [date_label, date_label], colour)


def graph_phase_coverage(plt_ax, p, infolder):
    """
    :param plt_ax: Pyplot axis instance to add the covered phase range to.
    :param p: Float. Period to phase by.
    :param infolder: String. Observation folder to get JDs from FITS headers.
    :return: Null.
    """
    phase_list = []
    date_list = []
    for filename in os.listdir(infolder):
        # Calculate phase location of each observation.
        fits_file = fits.open(infolder + '/' + filename)
        phase_list.append(phase(fits_file[0].header['JD'], p))

        # Get UTC MMDD from header for y-axis labels.
        utc_split = fits_file[0].header['DATE-OBS'].split('-')
        date_label = utc_split[1] + utc_split[2][0:2]
        date_list.append(date_label)
        fits_file.close()

    plt_ax.scatter(phase_list, date_list, s=1)


# Initialize plot.
fig1, ax1 = plt.subplots()
plt.title("BetaAur Period Coverage, P=3.96004 days")
plt.xlabel("Phase")
plt.ylabel("Day")
plt.xlim([0, 1])

# Build coverage plot from observation files.
period = 3.96004
obs_folder = 'BetaAurSp'
graph_phase_coverage(ax1, period, obs_folder)

# Add predictions to coverage plot.
add_prediction(ax1, period, 2459879.58333, 2459880.02083)
add_prediction(ax1, period, 2459880.58333, 2459881.02083, 'c')
add_prediction(ax1, period, 2459881.58333, 2459882.02083, 'g')
add_prediction(ax1, period, 2459882.58333, 2459883.02083, 'y')
plt.show()
