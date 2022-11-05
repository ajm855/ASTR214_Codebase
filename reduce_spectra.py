from astropy.io import fits
import os
import numpy as np

def read_calib_files():

    # Create a list of filenames of each calibration type
    biases = []
    darks = []
    flats = []
    lights = []

    # Search only in the current directory
    for filename in os.listdir('.'):
        try:
            fits_file = fits.open(filename)
            calib_type = ''

            # Handling for different fits header conventions
            try:
                try_read = fits_file[0].header['OBSTYPE'].strip().lower().split(' ')[0]
                calib_type = try_read
            except:
                try_read = fits_file[0].header['IMAGETYP'].strip().lower().split(' ')[0]
                calib_type = try_read

            if calib_type == 'flat':
                flats.append(filename)
            elif calib_type == 'dark':
                darks.append(filename)
            elif calib_type == 'bias':
                biases.append(filename)
            elif calib_type == 'light':
                lights.append(filename)
            fits_file.close()

        except Exception as e:
            pass
            #print(filename, ':', e)

    if len(biases) == 0: print("No bias frames found.")
    if len(darks) == 0: print("No dark frames found.")
    if len(flats) == 0: print("No flat frames found.")

    return biases, darks, flats, lights

def build_masters():
    biases, darks, flats, lights = read_calib_files()
    x = fits.open(biases[0])
    y = fits.getdata(biases[0])

    os.makedirs('./Calib Masters', exist_ok=True)

    bias_exists = len(biases) > 0
    dark_exists = len(darks) > 0
    flats_exist = len(flats) > 0

    bias_master = []
    dark_master = []
    flat_master = []

    if bias_exists:
        print('Creating master bias (median combine)...')

        bias_concat = [fits.getdata(image) for image in biases]
        bias_master = np.median(bias_concat, axis=0)

        hdu = fits.PrimaryHDU(bias_master)
        hdu.writeto('./Calib Masters/bias_master.fits', overwrite=True)
    else:
        print('Biases not found. Skipping bias master generation.')
    
    if dark_exists:
        print('Creating master dark (median combine)...')
        dark_concat = [fits.getdata(image) for image in darks]

        if bias_exists:
            print('     Subtracting biases...')
            dark_concat -= bias_master
        
        dark_master = np.median(dark_concat, axis=0)

        hdu = fits.PrimaryHDU(dark_master)
        hdu.writeto('./Calib Masters/dark_master.fits', overwrite=True)
    else:
        print('Darks not found. Skipping dark master generation.')
    
    if flats_exist:
        print('Creating master flat (median combine)...')
        flat_concat = [fits.getdata(image) for image in flats]
        
        if bias_exists:
            print('     Subtracting biases...')
            flat_concat -= bias_master
        if dark_exists:
            print('     Subtracting darks...')
            flat_concat -= dark_master

        flat_master = np.median(flat_concat, axis=0)

        hdu = fits.PrimaryHDU(flat_master)
        hdu.writeto('./Calib Masters/flat_master.fits', overwrite=True)
    else:
        print('Flats not found. Skipping flat master generation.')
    
    print('Calibration master generation process complete.')

    return bias_master, dark_master, flat_master, lights

def calibrate_light_frames():
    bias_master, dark_master, flat_master, lights = build_masters()

    os.makedirs('./Reduced Observations', exist_ok=True)

    for lightfile in lights:

        print('Processing:', lightfile)
        light_data = fits.getdata(lightfile)

        
        
        if len(bias_master) > 0:
            light_data = light_data - bias_master
        if len(dark_master) > 0:
            light_data = light_data - dark_master
        if len(flat_master) > 0:
            light_data = light_data / flat_master

        hdu = fits.PrimaryHDU(light_data)
        hdu.writeto('./Reduced Observations/'+"REDUCED-"+lightfile, overwrite=True)
        
calibrate_light_frames()







