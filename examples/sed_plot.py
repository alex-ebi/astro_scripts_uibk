"""
This is an example on how to make and SED fit.
"""

from pkg_resources import resource_filename
from astro_scripts_uibk import query, photometry, io_asu

# Define path to atlas SED file
sed_path = resource_filename("astro_scripts_uibk", "test_data/hd7902_sed.flux")

# Define star name
star_name = 'HD7902'

# Get magnitudes of the star via the query module
mag_df = query.sed_mags(star_name)

print(mag_df)

# If you have just the magnitudes, you can add the fluxes to our dataframe:
# Beware that the filters must be included in the directory "filter_profiles"
mag_df = query.apply_flux(mag_df)

print(mag_df)

sed = io_asu.read_atlas9_flux(sed_path)

print(sed)

photometry.sed_plotter(sed, mag_df, plot_filtered_sed=True, rv=3.16, ebv=0.58)
