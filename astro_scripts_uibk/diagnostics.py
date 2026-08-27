import astro_scripts_uibk as asu
import numpy as np

def signal_to_noise(spectrum: np.array, wave_center=None, wave_range=None, delta_wave=1):
    if wave_center is not None:
        wave_range = np.array([wave_center - delta_wave, wave_center + delta_wave])
    if wave_range is None:
        raise ValueError("Either wave_center or wave_range must be provided.")
    spec_slice = asu.spectrum_reduction.crop_spectrum(spectrum, wave_range[0], wave_range[1])
    signal = np.nanmean(spec_slice[1])
    noise = np.nanstd(spec_slice[1])

        
    return signal / noise