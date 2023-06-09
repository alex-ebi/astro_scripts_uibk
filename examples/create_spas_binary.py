"""
Example for creating a spas binary file from a spectrum.
"""
from pkg_resources import resource_filename
from astro_scripts_uibk import io_asu, spas
import numpy as np


def main():
    # Create microgrid:
    # If you only want to fit the abundance and you allready know the stellar parameters,
    # create a microgrid for each abundance in your grid.
    # Do this for every abundance.
    filename = resource_filename("astro_scripts_uibk", "test_data/T9000G1.50X7.0M02HE.100C8.41_flux")
    surface_spec = io_asu.read_surface_spectrum(filename)

    xmin = 3.653555e+03
    xmax = 3.657808e+03

    t_eff = 9000
    log_g = 1.5
    c = 8.41

    # if needed, crop/slice the spectrum to your desired wavelength range
    wave, flux = spas.crop_spectrum(surface_spec[0], surface_spec[1], xmin, xmax)
    spec = np.array([wave, flux])

    spas.microgrid(t_eff, log_g, c, spec, filename + '1.bin.tmp')  # in the real application use the file ending '.bin'

    # Create binary
    # If you want to make a spas binary grid with different stellar parameters, use 'spas.create_binary'.
    # We start with an example grid with two effective temperatures and two logg's.
    # This makes a grid with four grid points.

    wave = spec[0]  # Extract wavelength array. This array is used for all spectra in the grid.

    # Now we make a grid of fluxes for each grid point.
    # In this case we just replicate our dummy spectrum so we have is four times
    flux_array = np.repeat([spec[1]], 4, axis=0)  # Expand fluxes in microgrid
    print('Wave array:\n', wave)
    print('Flux array:\n', flux_array)

    # Now we make our T-eff and log_g grid.
    # It is very important that the grid is set up in the way stated below.
    # T_eff stays constant until all different log_g values are cycled.
    # Then we take the next T_eff value and log_g gets cycled again.
    t_eff_array = [9000, 9000, 10000, 10000]
    log_g_array = [1.5, 1.8, 1.5, 1.8]

    spas.create_binary(wave, flux_array, t_eff_array, log_g_array, c, filename+'2.bin.tmp')


if __name__ == '__main__':
    main()
