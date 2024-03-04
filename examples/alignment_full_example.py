import numpy as np
import astro_scripts_uibk as asu
from importlib.resources import files
import os
import pandas as pd
from astropy.io import fits
import random
import multiprocessing as mp

# Booleans for major working steps
download_spectra = True
perform_alignment = True
plot_results = True

# define path to spectra directory
my_spec_dir = files('astro_scripts_uibk') / 'example_data/spectra'

# Downloading the standard reduced EDIBLES sample from the ESO Science Archive
if download_spectra:
    download_script_path = files('astro_scripts_uibk') / 'example_data/spectra/download_script_uves_edibles_targets.sh'
    os.system(f'chmod +x {download_script_path}')  # chmod it
    os.chdir(download_script_path.parent)  # change working directory for os
    os.system(f'{download_script_path}')  # execute; already downloaded files are omitted by the bash script

my_index_path = files('astro_scripts_uibk') / 'example_data/spec_index.xlsx'


def read_spec(spec_file):
    hdul = fits.open(spec_file)
    hdr = hdul[0].header
    wave = hdul[1].data['WAVE'][0]
    flux = hdul[1].data['FLUX_REDUCED'][0]
    spec = np.array([wave, flux])
    return spec, hdr


def index_orders(spec_dir, index_path):
    setting_paths = spec_dir.glob('*.fits')

    df = pd.DataFrame()

    for path in setting_paths:
        spec, hdr = read_spec(spec_dir / path)

        x_limits = [min(spec[0]), max(spec[0])]

        obs_time_str = hdr['DATE-OBS']

        star_name = hdr['OBJECT']

        # simbad_star_name = star_name_dict.get(star_name)
        # if simbad_star_name is not None:
        #     star_name = simbad_star_name

        set_list = [star_name, obs_time_str, str(path), x_limits[0], x_limits[1]]
        row = pd.Series(set_list)
        df = pd.concat((df, row), axis=1, ignore_index=True)

    df = df.T

    df.rename(columns={0: 'star_name', 1: 'obs_date', 2: 'spec_path', 3: 'x_min', 4: 'x_max'}, inplace=True)

    # df = drop_spectra(df, exclude_ref_sight_lines)

    df.to_excel(index_path)

    return df


# Make an index file for the spectra
index_orders(my_spec_dir, my_index_path)

# Apply DIB alignment ======================================================================

cluster_output_path = files('astro_scripts_uibk') / 'example_data/clustering_results'
my_plot_dir = files('astro_scripts_uibk') / 'example_data/clustering_plots'

my_literature_dib_path = files("astro_scripts_uibk") / 'example_data/literature_dibs_2024.xlsx'

plot_query = False
test_run = False

data_index = pd.read_excel(my_index_path, index_col=0)

max_dist = 0.5  # Define maximal distance between Query and Subject for a match

lit_dibs = pd.read_excel(my_literature_dib_path, index_col=[0])

single_cloud_sightlines = ['HD23180', 'HD24398', 'HD144470', 'HD147165', 'HD147683', 'HD149757', 'HD166937', 'HD170740',
                           'HD184915', 'HD185418', 'HD185859', 'HD203532']

query_keys = [6379]  # run DIB alignment only for this DIB

sm_ratio = 0.1  # Define the width of the smoothing kernel relative to the DIB width


def io_function(spec_path):
    """
    Function to read the input spectra for DIB alignment.
    The format has to be in wave numbers and not wavelengths.

    Parameters
    ----------
    spec_path : str
        Path of the spectrum.

    Returns
    -------
    np.array
        Numpy array of spectrum.
    """
    spec, _ = read_spec(spec_path)
    spec[0] = asu.transformations.angstrom_to_wavenumber(spec[0])
    return spec


if __name__ == '__main__' and perform_alignment:
    for query_key in query_keys:
        spec_aligner = asu.alignment.SpectralAligner(lit_dibs, query_key, data_index,
                                                     io_function=io_function, spec_dir=my_spec_dir,
                                                     test_run=test_run, max_dist=max_dist,
                                                     plot_query=plot_query, smoothing_range_ratio=sm_ratio)

        iter_vars = spec_aligner.iter_vars
        node_count = 10

        print('Iter_var length:', len(iter_vars))
        if not test_run:
            random.shuffle(iter_vars)

        if plot_query or test_run:
            # Use loop for troubleshooting - better error traceback
            for iter_var in iter_vars:
                spec_aligner.one_query_analysis(*iter_var)

        else:
            # Step 1: Init multiprocessing.Pool()
            pool = mp.Pool(node_count)

            # Step 2: `pool.apply` the `howmany_within_range()`
            results = pool.starmap(spec_aligner.one_query_analysis, iter_vars)

            # Step 3: Don't forget to close
            pool.close()

            print('Results:', results)
            if len(results) > 0:
                result_df = pd.concat(results, ignore_index=True)
                print('Result DataFrame:', result_df)

                # List of "peaks" of the lowest distances between Query and Subject
                result_df.to_pickle(cluster_output_path / f'peaks_{query_key}.p')

# Plot some of the results
ang_range = [6088, 6092]  # Only plot the clusters in this Angstrom range; can be set to None

if __name__ == '__main__' and plot_results:
    for query in query_keys:
        asu.alignment.auto_plot_clusters(io_function=io_function, spec_dir=my_spec_dir,
                                         result_file=cluster_output_path / f'peaks_{query}.p',
                                         query_key=query,
                                         match_dist_cut=0.5,
                                         pearson_r_cut=0.90,
                                         sort_by='pearson_r',
                                         ascending=False,
                                         padding_factor=0.4,
                                         plot_dir=my_plot_dir,
                                         mem_range=0.5,
                                         eps=0.2,
                                         show=False,
                                         min_cluster_size=50,
                                         sm_ratio=sm_ratio,
                                         smooth=True,
                                         single_cloud_sightlines=single_cloud_sightlines,
                                         ang_range=ang_range
                                         )
