"""
Simple usage example for PyReduce.
Loads a sample dataset, and runs the full extraction for ECHELLE_UIBK.
Before running this example, pull astro-pyreduce and insert the ECHELLE_UIBK config files in the pyreduce repo and
install pyreduce.
The frames do NOT have to be flipped!

Location of config files: examples/pyreduce

Clone PyReduce:
git clone https://github.com/AWehrhahn/PyReduce

"""
from pyreduce.configuration import get_configuration_for_instrument
from pyreduce.reduce import Reducer
from pyreduce.util import start_logging
from pathlib import Path
from pprint import pprint
from astropy.io import fits
from astropy.coordinates import SkyCoord


# from pyreduce.combine_frames import combine_frames


def add_header_info(filelist: list, star_name: str):
    coord = SkyCoord.from_name(star_name)

    for fits_file in filelist:
        with fits.open(fits_file, 'update') as f:
            hdr = f[0].header
            dict_keys = list(hdr.keys())
            if 'OBJECT' not in dict_keys:
                hdr['OBJECT'] = star_name
            if 'RA--TAN' not in dict_keys:
                hdr['RA--TAN'] = coord.ra.deg
            if 'DEC-TAN' not in dict_keys:
                hdr['DEC-TAN'] = coord.dec.deg
            # Instrument Location
            if 'GEOLON' not in dict_keys:
                hdr['GEOLON'] = 11.400375
            if 'GEOLAT' not in dict_keys:
                hdr['GEOLAT'] = 47.259659
            if 'GEOELEV' not in dict_keys:
                hdr['GEOELEV'] = 617


instrument = 'ECHELLE_UIBK'  # Innsbruck echelle
target = 'polaris'  # This Keyword is used to search for science frames

# Define your working directory
w_dir = Path(r"C:\Users\Public\Documents\Echelle")
# We define the path to the output directory
output_dir = w_dir / "output"

# Collect paths of calibration and science frames
bias_paths = []
flat_paths = []
lamp_paths = []
science_paths = []

for path in (w_dir / 'raw').glob('*.fits'):
    if path.match('*bias*'):
        bias_paths.append(str(path))
    elif path.match('*flat*'):
        flat_paths.append(str(path))
    elif path.match('*lamp*'):
        lamp_paths.append(str(path))
    elif path.match(f'*{target}*'):
        science_paths.append(str(path))
    else:
        continue

# Optional: Combine science frames to reduce signal to noise
# You can also reduce every Science frame individually
# Not done yet
# master_science_path = output_dir / f'master_{target}.fits'
# master_science = combine_frames(lamp_paths, instrument, '')
# fits.writeto(master_science_path, master_science[0].data, overwrite=True)
# science_paths = [master_science_path]

# Truncate the science file list to one file
science_paths = science_paths[:1]

# If not specified, add necessary header info
add_header_info(science_paths, target)

# Print paths
print('Bias paths')
pprint(bias_paths)

print('Flat paths')
pprint(flat_paths)

print('Lamp paths')
pprint(lamp_paths)

print('Science paths')
pprint(science_paths)

# Now we load the instrument configuration of ECHELLE_UIBK.
# This way you can change config parameters.
config = get_configuration_for_instrument(instrument, plot=1)

# Define which steps you want to plot

# config['bias']['plot'] = False
# config['flat']['plot'] = False
# config['orders']['plot'] = False
# config['norm_flat']['plot'] = False
# config['curvature']['plot'] = False
# config['wavecal']['plot'] = False
# config['wavecal_master']['plot'] = False
# config['science']['plot'] = False
# config['continuum']['plot'] = False

# Redefine config parameters if you want to try things out
# The possible values can be found in examples/pyreduce/settings_ECHELLE_UIBK.json

# Examples:
# Degree of polynomials to detect orders
# config["orders"]["degree"] = 5

# Config for curvature detection
# config["curvature"]["peak_width"] = .5

# Since we can't find the files ourselves (at least not without defining the criteria we are looking for)
# We need to manually define which files go where
files = {"bias": bias_paths, "flat": flat_paths, "orders": flat_paths, "curvature": lamp_paths,
         "science": science_paths, "wavecal_master": lamp_paths}

# (optional) We need to define the log file
log_file = output_dir / "log_file.txt"
start_logging(log_file)

# Define other parameter for PyReduce
night = ""
mode = ""
steps = (
    "bias",
    "flat",
    "orders",
    "norm_flat",
    "curvature",
    # "scatter",
    "wavecal",
    "science",
    # "continuum",
    "finalize",
    # "wavecal_master",
    # # "freq_comb",
)

# Call the PyReduce algorithm
reducer = Reducer(
    files,
    str(output_dir),
    target,
    instrument,
    mode,
    night,
    config,
    # order_range=order_range,
    # skip_existing=False,
)
data = reducer.run_steps(steps=steps)
