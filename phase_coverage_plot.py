# Phase Coverage Graphing Tool
# Alexander Magnus
import os
import matplotlib.pyplot as plt
from suntimes import SunTimes
from datetime import date, datetime, timedelta
from astropy.io import fits
from astropy.time import Time


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

def utc_to_jd(utc_date):
    t = Time(utc_date, scale='utc')
    return t.jd

def find_phase_night(target_phase_at_sunset, p, threshold=0.01):
    """
    :param target_phase_at_sunset: Float. Desired phase at sunset for given object with period p.
    :param p: Float. Period to phase by.
    :param threshold: Float. Max acceptable absolute difference between desired phase and iterated phase during search.
    :return: datetime date. First date found at which the day's phase at sunset is within threshold of target_phase_at_sunset. 
    Returns -1 if nothing is found within one year.
    """
    # Setup Saskatoon sun
    latitude = 52.146973
    longitude = -106.647034
    altitude = 482 # meters
    sun = SunTimes(latitude=latitude, longitude=longitude, altitude=altitude)

    # Start iterating from today, end on the same date next year.
    date_iter = datetime.utcnow().date()
    next_year = date(date_iter.year+1, date_iter.month, date_iter.day)
    
    while date_iter < next_year :
        # Calculate time of sunset in UTC at Saskatoon
        sunset = sun.setutc(date_iter)
        #Calculate phase of object at sunset
        iter_phase_at_sunset = phase(utc_to_jd(sunset), p)

        # If phase is within threshold of desired phase, return date.
        if abs(iter_phase_at_sunset - target_phase_at_sunset) < threshold:
            print("The object's phase will be", iter_phase_at_sunset, "at sunset on", date_iter)
            return date_iter
        date_iter += timedelta(days=1)
    print("Couldn't find required phase within a year.")
    return -1


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
#plt.show()

# Saskatoon
find_phase_night(0.5, 3.96004)

print()

