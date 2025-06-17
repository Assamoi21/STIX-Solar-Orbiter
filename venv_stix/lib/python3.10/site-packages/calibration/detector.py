import astropy.units as u
import numpy as np
from types import SimpleNamespace

from sunpy.io.special import read_genx
from stixcore.calibration.transmission import Transmission


def get_srm():
    r"""
    Return the spectromber response matrix (SRM) by combing the attenuation with the detetor respoonse matrix (DRM)

    Returns
    -------

    """
    drm_save = read_genx('/Users/shane/Projects/STIX/git/stix_drm_20220713.genx')
    drm = drm_save['SAVEGEN0']['SMATRIX'] * u.count/u.keV/u.photon
    energies_in = drm_save['SAVEGEN0']['EDGES_IN'] * u.keV
    energies_in_width = np.diff(energies_in)
    energies_in_mean = energies_in[:-1] + energies_in_width/2
    trans = Transmission()
    tot_trans = trans.get_transmission(energies=energies_in_mean)
    energies_out = drm_save['SAVEGEN0']['EDGES_OUT'] * u.keV
    energies_out_width = drm_save['SAVEGEN0']['EWIDTH'] * u.keV
    energies_out_mean = drm_save['SAVEGEN0']['EMEAN'] * u.keV
    attenuation = tot_trans['det-0']
    srm = (attenuation.reshape(-1, 1) * drm * energies_out_width) / 4  # avg grid transmission

    res = SimpleNamespace(drm=drm, srm=srm, attenuation=attenuation, energies_in=energies_in,
                          energies_in_width=energies_in_width, energies_in_mean=energies_in_mean,
                          energies_out=energies_out, energies_out_width=energies_out_width,
                          energies_out_mean=energies_out_mean, area=1*u.cm)
    return res


def get_pixel_srm():
    pass
