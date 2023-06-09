# -*- coding: utf-8 -*-
"""
Handles instrument specific info for the UVES spectrograph

Mostly reading data from the header
"""
import logging
import os.path

import numpy as np
from astropy.io import fits
from dateutil import parser

from .common import Instrument, getter, observation_date_to_night

logger = logging.getLogger(__name__)


class ECHELLE_UIBK(Instrument):
    def add_header_info(self, header, mode, **kwargs):
        """read data from header and add it as REDUCE keyword back to the header"""
        # "Normal" stuff is handled by the general version, specific changes to values happen here
        # alternatively you can implement all of it here, whatever works
        header = super().add_header_info(header, mode)


        # header["e_ra"] /= 15
        # if header["e_jd"] is not None:
        #     header["e_jd"] += header["e_exptime"] / (7200 * 24) + 0.5

        return header

    def get_wavecal_filename(self, header, mode, **kwargs):
        """Get the filename of the wavelength calibration config file"""
        info = self.load_info()
        # specifier = int(header[info["wavecal_specifier"]])

        cwd = os.path.dirname(__file__)
        fname = "{instrument}_2D.npz".format(
            instrument="echelle_uibk"
        )
        
        # Option if modes are added
        # fname = "{instrument}_{mode}_{specifier}nm_2D.npz".format(
        #     instrument="echelle_uibk", mode=mode.lower(), specifier=specifier
        # )
        fname = os.path.join(cwd, "..", "wavecal", fname)
        return fname

    def get_mask_filename(self, mode, **kwargs):
        i = self.name.lower()
        m = mode.lower()
        fname = f"mask_{i}.fits"
        cwd = os.path.dirname(__file__)
        fname = os.path.join(cwd, "..", "masks", fname)
        return fname
