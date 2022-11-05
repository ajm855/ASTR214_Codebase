# Phase Coverage Graphing Tool
# Alexander Magnus
import os
import matplotlib.pyplot as plt
import numpy as np
from suntimes import SunTimes
from datetime import date, datetime, timedelta
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz


def phase(jd, p):
    """
    :param jd: Float. Julian date to be phased.
    :param p: Float. Period to phase by.
    :return: Float. Phase of observation.
    """
    return (jd % p) / p

def utc_to_jd(utc_date):
    """
    :param utc_date: Julian date to be converted.
    :return: Float. Converted Julian date.
    """
    t = Time(utc_date, scale='utc')
    return t.jd

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

    plt_ax.scatter(phase_list, date_list, s=1, zorder=4)

def find_phase_at_sunset(target_phase_at_sunset, p, sky_obj, threshold=0.01, verbose=True):
    """
    :param target_phase_at_sunset: Float. Desired phase at sunset for given object with period p.
    :param p: Float. Period to phase by.
    :param threshold: Float. Max acceptable absolute difference between desired phase and iterated phase during search.
    :return: datetime date. First date found at which the day's phase at sunset is within threshold of target_phase_at_sunset. 
    Returns -1 if nothing is found within one year.
    """
    # Set up Saskatoon sun for sunset time calculation.
    latitude = 52.146973
    longitude = -106.647034
    altitude = 482 # meters
    sun = SunTimes(latitude=latitude, longitude=longitude, altitude=altitude)

    # Set up observatory location and object information for alitude/elevation calculation.
    obj = sky_obj
    saskatoon = EarthLocation(lat= latitude*u.deg, lon=longitude*u.deg, height=altitude*u.m )

    # Start iterating from today, end on the same date next year.
    date_iter = datetime.utcnow().date()
    next_year = date(date_iter.year+1, date_iter.month, date_iter.day)
    
    while date_iter < next_year :
        # Calculate time of sunset in UTC at Saskatoon.
        sunset = sun.setutc(date_iter)

        #Calculate phase of object at sunset.
        iter_phase_at_sunset = phase(utc_to_jd(sunset), p)

        # If phase is within threshold of desired phase, continue. 
        # Handling for cases where phase is near 0 and 1, threshold needs to "wrap" over.
        wrap_up = iter_phase_at_sunset + threshold - 1
        wrap_down = iter_phase_at_sunset - threshold + 1
        if (
            (abs(iter_phase_at_sunset - target_phase_at_sunset) < threshold) or 
            (wrap_up > 0 and wrap_up > target_phase_at_sunset) or 
            (wrap_down < 1 and wrap_down < target_phase_at_sunset)
        ):

            # Calculate altitude of given object at sunset at observatory location.
            time = sunset
            obj_altaz = obj.transform_to(AltAz(obstime=time, location=saskatoon))
            obj_alt = float(f"{obj_altaz.alt}"[:-4])
            
            # If altitude is greater than 20.0 degrees, return JD of sunrise and sunset.
            if obj_alt > 20.0:
                tomorrow = date_iter + timedelta(days=1)
                sunrise = sun.riseutc(tomorrow)
                iter_phase_at_sunrise = phase(utc_to_jd(sunrise), p)

                if verbose:
                    print(date_iter, tomorrow)
                    print("At sunset, the object's phase will be", iter_phase_at_sunset, "occuring on", sunset,"UTC. It's altitude at this time is", obj_alt, "deg."  )
                    print("At sunrise the following morning, the object's phase will be", iter_phase_at_sunrise, "occuring on", sunrise, "UTC.")
                return [utc_to_jd(sunset), utc_to_jd(sunrise)]
        # Iterate date.
        date_iter += timedelta(days=1)
    
    # No sunset within phase threshold and elevation constraints for this observatory over the next year. 
    print("Couldn't find required phase within threshold and elevation constraints within a year.")
    return [-1,-1]

def phase_coverage_tonight(p, sky_obj, threshold=0.01, verbose=True):
    """
    :param target_phase_at_sunset: Float. Desired phase at sunset for given object with period p.
    :param p: Float. Period to phase by.
    :param threshold: Float. Max acceptable absolute difference between desired phase and iterated phase during search.
    :return: datetime date. First date found at which the day's phase at sunset is within threshold of target_phase_at_sunset. 
    Returns -1 if nothing is found within one year.
    """
    # Set up Saskatoon sun for sunset time calculation.
    latitude = 52.146973
    longitude = -106.647034
    altitude = 482 # meters
    sun = SunTimes(latitude=latitude, longitude=longitude, altitude=altitude)

    # Set up observatory location and object information for alitude/elevation calculation.
    obj = sky_obj
    saskatoon = EarthLocation(lat= latitude*u.deg, lon=longitude*u.deg, height=altitude*u.m )
    
    today = datetime.utcnow().date()
    tomorrow = today + timedelta(days=1)

    sunset = sun.setutc(today)
    sunrise = sun.riseutc(tomorrow)
    phase_at_sunset = phase(utc_to_jd(sunset), p)
    phase_at_sunrise = phase(utc_to_jd(sunrise), p)

    obj_altaz = obj.transform_to(AltAz(obstime=sunset, location=saskatoon))
    obj_alt = float(f"{obj_altaz.alt}"[:-4])

    print("At sunset, the object's phase will be", phase_at_sunset, "occuring on", sunset,"UTC. It's altitude at this time is", obj_alt, "deg."  )
    print("At sunrise the following morning, the object's phase will be", phase_at_sunrise, "occuring on", sunrise, "UTC.")
    return [utc_to_jd(sunset), utc_to_jd(sunrise)]


# Initialize plot.
fig1, ax1 = plt.subplots()
ax1.minorticks_on()
plt.title("BetaAur Period Coverage, P=3.96004 days")
plt.xlabel("Phase")
plt.xlim([0, 1])
plt.xticks(np.linspace(0, 1, 11))
plt.ylabel("Date, MMDD (UTC)")

plt.grid(visible=True)

# Build coverage plot from observation files.
period = 3.96004
obs_folder = 'BetaAurSp'
graph_phase_coverage(ax1, period, obs_folder)

sunrise, sunset = phase_coverage_tonight(period, SkyCoord.from_name("beta aurigae"))
add_prediction(ax1, period, sunrise, sunset)

# Add predictions to coverage plot.
#add_prediction(ax1, period, 2459879.58333, 2459880.02083)
#add_prediction(ax1, period, 2459880.58333, 2459881.02083, 'c')
#add_prediction(ax1, period, 2459881.58333, 2459882.02083, 'g')
#add_prediction(ax1, period, 2459882.58333, 2459883.02083, 'y')

# Search for desired phase at sunset, add predictions to plot.
#results = find_phase_at_sunset(0.1, 3.96004, SkyCoord.from_name("beta aurigae"))
#add_prediction(ax1, period, results[0], results[1])
plt.show()

# Testing for altitude calculation.
"""
ba = SkyCoord.from_name('beta aurigae')
saskatoon = EarthLocation(lat= 52.146973*u.deg, lon=-106.647034*u.deg, height=482*u.m )
time = Time.now()

ba_altaz = ba.transform_to(AltAz(obstime=time, location=saskatoon))
print(f"BetaAur's Altitude = {ba_altaz.alt:.4}")
print(ba_altaz.alt)
print(float(f"{ba_altaz.alt}"[:-4]))
"""