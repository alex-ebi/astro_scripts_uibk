import dataclasses
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import SpanSelector
from astro_scripts_uibk import spectrum_reduction, pub_plot
from matplotlib.backend_bases import MouseButton
from scipy.interpolate import CubicSpline


@dataclasses.dataclass
class EWContainer:
    ew: float
    ew_error: float


def ew_measurement(spectrum: np.array, resolution=100000, cont_range=1):
    """
    Widget for measuring the equivalent width of spectra.
    A linear continuum is calculated by selecting two continuum ranges of width cont_range.
    Then, an equivalent width is calculated using this linear continuum by direct integration over the spectrum.

    Parameters
    ----------
    spectrum : np.array
        Input spectrum.
    resolution : float
        Resolution of the spectrograph.
    cont_range : float
        Width of the continuum ranges.

    Returns
    -------
    (float, float)
        Equivalent width and error of equivalent width.
    """
    ew_set = EWContainer
    f, ax = plt.subplots(figsize=[10, 6])

    span1 = plt.axvspan(min(spectrum[0]) - cont_range * .5, min(spectrum[0]) + cont_range * .5,
                        alpha=.5, color='orange')
    span2 = plt.axvspan(max(spectrum[0]) - cont_range * .5, max(spectrum[0]) + cont_range * .5,
                        alpha=.5, color='orange')

    def onselect(x, y):
        s1min = x - cont_range * .5
        s1max = x + cont_range * .5
        span1.set_xy(np.array([[s1min] * 2 + [s1max] * 2, span1.xy[:-1, 1]]).T)

        s2min = y - cont_range * .5
        s2max = y + cont_range * .5
        span2.set_xy(np.array([[s2min] * 2 + [s2max] * 2, span2.xy[:-1, 1]]).T)

        # calculate normalized spectrum and weights of continuum regions (which is (S/N))
        norm_spec, w = spectrum_reduction.normalize_mean_flux(spectrum, x, y, cont_range=cont_range, return_weight=True)
        # crop spectrum to integrated area
        int_spec = spectrum_reduction.crop_spectrum(norm_spec, x, y)
        # EW
        ew_set.ew = np.trapz(1 - int_spec[1], x=int_spec[0])
        # EW error
        spec_disp = np.mean([x, y]) / resolution  # calculate spectral dispersion
        ew_range = y - x
        ew_set.ew_error = np.sqrt(2 * ew_range * spec_disp) / w

    plt.plot(spectrum[0], spectrum[1])

    rs = SpanSelector(ax, onselect, 'horizontal', props=dict(alpha=0.5, facecolor="tab:blue"))
    rs.set_active(True)
    plt.show()

    return ew_set.ew, ew_set.ew_error


@dataclasses.dataclass
class SpectrumContainer:
    spectrum: np.array
    weight: float = None


def normalize(spectrum: np.array, cont_range=1, return_weight=False):
    """
    Widget for normalizing spectra.
    A linear continuum is calculated by selecting two continuum ranges of width cont_range.
    On quitting the normalization window, the normalized spectrum, and optionally a weight, i.e. S/N ratio is returned.

    Parameters
    ----------
    spectrum : np.array
        Input spectrum.
    cont_range : float
        Width of the continuum ranges.
    return_weight : bool
        If True, the function additionally returns the S/N ratio calculated in the continuum ranges.

    Returns
    -------
    np.array or (np.array, float)
        Normalized spectrum or (Normalized spectrum, S/N)
    """
    out_spec = SpectrumContainer
    f, ax = plt.subplots(figsize=[10, 6])

    span1 = plt.axvspan(min(spectrum[0]) - cont_range * .5, min(spectrum[0]) + cont_range * .5,
                        alpha=.5, color='orange')
    span2 = plt.axvspan(max(spectrum[0]) - cont_range * .5, max(spectrum[0]) + cont_range * .5,
                        alpha=.5, color='orange')

    spec_plot, = plt.plot(spectrum[0], spectrum[1])

    def onselect(x, y):
        s1min = x - cont_range * .5
        s1max = x + cont_range * .5
        span1.set_xy(np.array([[s1min] * 2 + [s1max] * 2, span1.xy[:-1, 1]]).T)

        s2min = y - cont_range * .5
        s2max = y + cont_range * .5
        span2.set_xy(np.array([[s2min] * 2 + [s2max] * 2, span2.xy[:-1, 1]]).T)

        # calculate normalized spectrum and weights of continuum regions (which is (S/N))
        if return_weight:
            norm_spec, w = spectrum_reduction.normalize_mean_flux(spectrum, x, y, cont_range=cont_range,
                                                                  return_weight=return_weight)
            out_spec.spectrum = norm_spec
            out_spec.weight = w
        else:
            norm_spec = spectrum_reduction.normalize_mean_flux(spectrum, x, y, cont_range=cont_range,
                                                               return_weight=return_weight)
            out_spec.spectrum = norm_spec
        spec_plot.set_ydata(norm_spec[1])

    rs = SpanSelector(ax, onselect, 'horizontal', props=dict(alpha=0.5, facecolor="tab:blue"))
    rs.set_active(True)
    plt.show()
    if return_weight:
        return out_spec.spectrum, out_spec.weight
    else:
        return out_spec.spectrum


def mark_molecfit_ranges(ax, include_list: list = None):
    if include_list is None:
        include_list = []

    def onselect(x, y):
        print(x, y)
        plt.axvspan(x, y, alpha=.5, color='orange')
        include_list.append([x, y])

    def onselect_del(x, y):
        for i, row in enumerate(include_list):
            if row[0] < x < row[1]:
                plt.axvspan(row[0], row[1], alpha=.5, color='blue')
                include_list.pop(i)

    for row in include_list:
        plt.axvspan(row[0], row[1], alpha=.5, color='orange')

    _ = SpanSelector(ax, onselect, 'horizontal', props=dict(alpha=0.5, facecolor="tab:blue"), button=MouseButton(1))
    _ = SpanSelector(ax, onselect_del, 'horizontal', props=dict(alpha=0.5, facecolor="tab:green"),
                     button=MouseButton(3))

    plt.show()

    return include_list

def continuum_fit_func(spectrum: np.array):
    """
    The continuum is fitted using a cubic spline, with manually-selected
    anchor points.
    
    A seperate popup window opens for interactive continuum fitting.

    Clicking adds an anchor point at the position of the cursor. Right click
    to remove the last anchor point added.

    After closing the fitting window, the continuum is applied unless anchor points > 2:
    then no changes made.

    The GUI plots are updated automatically after each window is processed.
    """


    wave = spectrum[0]
    flux = spectrum[1]

    cont_anchors = []
    cont_artists = []
    preview = [None]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(wave, flux, 'k-')
    ax.set_title(f'Continuum normalizer - Left-click will add an anchor; right-click will remove last anchor.')
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Flux')

    def onclick(event, wave=wave, ax=ax, fig=fig, cont_anchors=cont_anchors, preview=preview, cont_artists=cont_artists):
        if event.xdata is None or event.ydata is None:
            return
        if event.button == 1:
            lam, flx = event.xdata, event.ydata
            cont_anchors.append((lam, flx))
            dot, = ax.plot(lam, flx, 'o', color='orange', ms=6)
            cont_artists.append(dot)
            if len(cont_anchors) >= 2:
                anc = sorted(cont_anchors)
                cs = CubicSpline([p[0] for p in anc], [p[1] for p in anc], extrapolate=True)
                if preview[0]:
                    try: preview[0].remove()
                    except: pass
                preview[0], = ax.plot(wave, cs(wave), 'orange', alpha=0.7)
            ax.set_title(f'{len(cont_anchors)} anchors. Close when done')
            fig.canvas.draw()

        elif event.button == 3 and cont_anchors:
            cont_anchors.pop()
            cont_artists.pop().remove()
            if len(cont_anchors) >= 2:
                anc = sorted(cont_anchors)
                cs = CubicSpline([p[0] for p in anc], [p[1] for p in anc], extrapolate=True)
                if preview[0]:
                    preview[0].remove()
                preview[0], = ax.plot(wave, cs(wave), 'orange', alpha=0.7)
            else:
                if preview[0]:
                    preview[0].remove()
                    preview[0] = None
            ax.set_title(f'{len(cont_anchors)} anchors. Close when done')
            fig.canvas.draw()

    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.tight_layout()
    plt.show(block=True)
        # plt.close(fig)

    if len(cont_anchors) >= 2:
        anc = sorted(cont_anchors)
        spline = CubicSpline([p[0] for p in anc], [p[1] for p in anc], extrapolate=True)
        continuum = spline(wave)
        spectrum[1] /= continuum

        if len(spectrum) > 2:
            spectrum[2] /= continuum

    return spectrum, continuum
