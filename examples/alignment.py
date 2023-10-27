from importlib.resources import files
import multiprocessing as mp
import astro_scripts_uibk as asu
import pandas as pd
import random
from pathlib import Path
from astropy.io import fits
import numpy as np


# Define some functions made for the WINERED example_data set
def read_spec(abs_path, return_header=False):
    abs_path = Path(abs_path)
    hdul = fits.open(abs_path)
    flux = hdul[0].data
    hdr = hdul[0].header
    wave = asu.io_asu.wave_from_dispersion(flux, start_wave=hdr['CRVAL1'], dispersion=hdr['CDELT1'],
                                           crpix=hdr['CRPIX1'])

    if return_header:
        return hdr, np.array([wave, flux], dtype=float)
    else:
        return np.array([wave, flux], dtype=float)


def order_from_string(str_in: str):
    print(str_in)
    order_start = str_in.find('_sum_')
    order_str = str_in[order_start + 6:order_start + 8]
    if order_start == -1:
        order_start = str_in.find('T_')
        order_str = str_in[order_start + 3:order_start + 5]
    # order
    order = int(order_str)
    print(order)

    return order


def index_orders(spec_dir=None,
                 setting='wide', index_path=None):
    setting_paths = spec_dir.rglob('*.fits')
    setting_paths = [item.relative_to(spec_dir) for item in setting_paths if not item.name.startswith('._')]
    setting_paths = [item for item in setting_paths if not item.match('*norm.fits')]

    star_name_dict = {'CygOB2No.12': 'Cyg OB2 12',
                      'CygOB2No3': 'Cyg OB2 3',
                      'CygOB2No5': 'Cyg OB2 5',
                      'CygOB2No8A': 'Cyg OB2 8A',
                      'CygOB2No9': 'Cyg OB2 9',
                      'CygOB2No10': 'Cyg OB2 10',
                      'CygOB2No11': 'Cyg OB2 11',
                      'Westerlund 1 W 33': 'Cl* Westerlund 1 W 33'}

    df = pd.DataFrame()

    for path in setting_paths:
        order = order_from_string(path.name)

        hdr, spec = read_spec(spec_dir / path, return_header=True)

        x_limits = [min(spec[0]), max(spec[0])]

        obs_time_str = hdr['ACQTIME1']

        star_name = hdr['OBJECT']

        simbad_star_name = star_name_dict.get(star_name)
        if simbad_star_name is not None:
            star_name = simbad_star_name

        set_list = [setting, order, star_name, obs_time_str, str(path), x_limits[0], x_limits[1]]
        row = pd.Series(set_list)
        df = pd.concat((df, row), axis=1, ignore_index=True)

    df = df.T

    df.rename(columns={0: 'setting', 1: 'order', 2: 'star_name', 3: 'obs_date', 4: 'spec_path', 5: 'x_min', 6: 'x_max'},
              inplace=True)

    # df = drop_spectra(df, exclude_ref_sight_lines)

    df.to_excel(index_path)

    return df


# Apply DIB alignment ======================================================================

my_index_path = Path('path/to/order_spec_index.xlsx')
cluster_output_path = Path('path/to/output_dir')
my_spec_dir = Path('path/to/spectra')
my_plot_dir = Path('path/to/my/plot_dir')

my_literature_dib_path = files("astro_scripts_uibk") / 'example_data/literature_dibs_2023_nir.xlsx'

plot_query = False
test_run = False

index_orders(spec_dir=my_spec_dir, index_path=my_index_path)  # Run this function to index your spectra

data_index = pd.read_excel(my_index_path, index_col=0)

max_dist = 0.5  # Define maximal distance between Query and Subject for a match

lit_dibs = pd.read_excel(my_literature_dib_path, index_col=[0])

query_keys = lit_dibs.index
query_keys = [10780]  # uncomment if you want to run only one query.

sm_ratio = 0.2


def io_function(spec_path):
    spec = read_spec(spec_path)
    spec[0] = asu.transformations.angstrom_to_wavenumber(spec[0])
    return spec


if __name__ == '__main__':
    for query_key in query_keys:
        spec_aligner = asu.alignment.SpectralAligner(lit_dibs, query_key, data_index,
                                                     io_function=io_function, spec_dir=my_spec_dir,
                                                     test_run=test_run, max_dist=max_dist,
                                                     plot_query=plot_query)

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
if __name__ == '__main__':
    # query_list = [10780, 13175]
    lit_dibs = pd.read_excel(my_literature_dib_path, index_col=[0])
    query_list = lit_dibs.index
    query_list = [10780]

    print(query_list)
    for query in query_list:
        asu.alignment.auto_plot_clusters(io_function=read_spec, spec_dir=my_spec_dir,
                                         result_file=cluster_output_path / f'peaks_{query}.p',
                                         query_key=query,
                                         match_dist_cut=0.5,
                                         pearson_r_cut=0.7,
                                         sort_by='pearson_r',
                                         ascending=False,
                                         padding_factor=0.4,
                                         plot_dir=my_plot_dir,
                                         mem_range=0.5,
                                         eps=0.2,
                                         show=False,
                                         min_cluster_size=10,
                                         sm_ratio=sm_ratio,
                                         smooth=True
                                         )
