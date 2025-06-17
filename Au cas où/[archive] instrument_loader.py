"""
The following code is for instrument specific classes each using their own methods to create and edit their
`_loaded_spec_data` attrbutes.

Tips that have been followed:
    * None of the instrument loaders should have public attributes; ie., all attributes should be preceded with `_`
    * Only obvious and useful methods and setters should be public, all else preceded with `_`

Credits to Natalia Bajnokova, Kris Cooper & Dan Ryan.
For original codes & spectral files info, see https://github.com/KriSun95/sunxspex/tree/master/sunxspex/sunxspex_fitting
"""

import numpy as np
import get_file_info

__all__ = ["StixLoader"]


# Get a default class for the instrument specfic loaders
# Once the instrument specific loaders inherit from this then all they really have to do is get the spectral
#    data they want to fit in the correct dictionary form and assigned to `self._loaded_spec_data`.

class StixLoader:
    """
    Loader specifically for STIX spectral data.

    StixLoader Specifics
    ----------------------
    Has methods to plot time series and perform time selection on the data. A background time can be added or removed
    and can fit the event data with the model+background (recommended) or fit a model to data-background using the
    `data2data_minus_background` setter with False or True, respectively.

    We assume that the background (if calculated) and the event emission is calculated from the same sized area. If this
    is False then the background effective exposure time should be multiplied by the ratio of the background area and
    the event area.
    I.e., background_effective_exposure * (background_area / event_area) as described in [1].
    This may be automated (16/03/2022).

    [1] https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html

    Properties
    ----------
    data2data_minus_background : bool
            Returns True if the data to be fitted has been converted to event time minus background.

    end_background_time : `astropy.Time`
            Returns the set backgroun end time.

    end_event_time : `astropy.Time`
            Returns the set event end time.

    start_background_time : `astropy.Time`
            Returns the set background start time.

    start_event_time : `astropy.Time`
            Returns the set event start time.

    Setters
    -------
    data2data_minus_background : bool
            Change the way the data is fitted; either event time minus background is fitted with the model (True) or
            event time is fitted with background+model (False, recommended and default behaviour).

    end_background_time : str, `astropy.Time`, None
            Sets the end time for the background time spectrum. None (default behaviour) sets this to None and removes
            or doesn't add a background component. The `start_background_time` setter can be assigned separately but
            must also be set to produce a background time.
            See `select_time` method to set both end and/or start time(s) in just one line.

    end_event_time : str, `astropy.Time`, None
            Sets the end time for the event time spectrum. None (default behaviour) sets this to the last time the
            loaded spectrum has data for.
            See `select_time` method to set both end and/or start time(s) in just one line.

    start_background_time : str, `astropy.Time`, None
            Sets the start time for the background time spectrum. None (default behaviour) sets this to None and removes
            or doesn't add a background component. The `end_background_time` setter can be assigned separately but must
            also be set to produce a background time.
            See `select_time` method to set both end and/or start time(s) in just one line.

    start_event_time : str, `astropy.Time`, None
            Sets the start time for the event time spectrum. None (default behaviour) sets this to the first time the
            loaded spectrum has data for.
            See `select_time` method to set both end and/or start time(s) in just one line.

    Methods
    -------
    lightcurve : energy_ranges (list of length=2, lists of length=2 lists, or None), axes (axes object or None)
            Plots the STIX time profile in the energy range given and on the axes provided. Default behaviour
            (energy_ranges=None) is to include all energies.

    select_time : start (str, `astropy.Time`, None), end (str, `astropy.Time`, None), background (bool)
            Set the start and end times for the event (background=False, default) or background (background=True) data.
            Both and end and a start time needs to be defined for the background, whereas the event time is assumed to
            commence/finish at the first/final data time if the start/end time is not given.

    spectrogram : axes (axes object or None) and any kwargs are passed to matplotlib.pyplot.imshow
            Plots the STIX spectrogram of all the data.

    Attributes
    ----------
    _channel_bins_inds_perspec : 2d array
            Array of channel bins (columns) per spectrum (rows).

    _construction_string : string
            String to show how class was constructed.

    _count_rate_perspec : 2d array
            Array of count rates per channel bin (columns) and spectrum (rows).

    _count_rate_error_perspec 2d array
            Array of count rate errors per channel bin (columns) and spectrum (rows).

    _counts_perspec : 2d array
            Array of counts per channel bin (columns) and spectrum (rows).

    _counts_err_perspec : 2d array
            Array of count error per channel bin (columns) and spectrum (rows).

    _end_background_time : `astropy.Time`
            End time for the defined background.
            Default: None

    _end_event_time : `astropy.Time`
            End time for the defined event.
            Default: Last time in loaded data.

    _full_obs_time : [`astropy.Time`, `astropy.Time`]
            Start and end time of the data loaded in.

    _lightcurve_data : {"mdtimes":_ts, "lightcurves":_lcs, "lightcurve_error":_lcs_err, "energy_ranges":energy_ranges}
            Arrays used to plots the lightcurves if the lightcurve method has been run.

    _loaded_spec_data : dict
            Instrument loaded spectral data.

    _lvt_perspec : 2d array
            Array of livetimes per channel bin (columns) and spectrum (rows).

    _spectrogram_data : {"spectrogram":_spect, "extent":_ext}
            Arrays used to plots the spectrogram if the spectrogram method has been run.

    _start_background_time
            Starting time for the defined background.
            Default: None

    _start_event_time
            Starting time for the defined event.
            Default: First time in loaded data.

    _time_bins_perspec : 2d array
            Array of time bins per spectrum.

    _time_fmt
            Format for astropy to convert teh times with.
            Default: "isot"

    _time_scale
            Scale for astropy to convert teh times with.
            Default: "utc"

    """
    __doc__ = """Parameters
                             ----------
                             pha_file : string
                                     The PHA file for the spectrum to be loaded.

                             arf_file, rmf_file : string
                                     The ARF and RMF files associated with the PHA file(s). If none are given (e.g, with
                                     NuSTAR data) it is assumed that these are in the same directory with same filename
                                     as the PHA file(s) but with extensions '.arf' and '.rmf', respectively.

                             srm_file : string
                                     The file that contains the spectral response matrix for the given spectrum.

                             srm_custom : 2d array
                                     User defined spectral response matrix. This is accepted over the SRM created from
                                     any ARF and RMF files given.

                             custom_channel_bins, custom_photon_bins : 2d array
                                     User defined channel bins for the columns and rows of the SRM matrix.
                                     E.g., custom_channel_bins=[[1,1.5],[1.5,2],...]

                             Attributes
                             ----------
                             _construction_string : string
                                     String to show how class was constructed.

                             _loaded_spec_data : dict
                                     Loaded spectral data.
                             """

    def __init__(self, pha_file, srm_file=None, srm_custom=None, custom_channel_bins=None, custom_photon_bins=None,
                 **kwargs):
        """Construct a string to show how the class was constructed (`_construction_string`) and set the
        `_loaded_spec_data` dictionary attribute."""

        self._lightcurve_data = []
        self._spectrogram_data = []
        self._full_obs_time = []
        self._channel_bins_inds_perspec = []
        self._time_bins_perspec = []
        self._lvt_perspec = []
        self._counts_perspec = []
        self._counts_err_perspec = []
        self._count_rate_perspec = []
        self._count_rate_error_perspec = []

        self._construction_string = f"StixLoader(pha_file={pha_file},srm_file={srm_file},srm_custom={srm_custom}," \
                                    f"custom_channel_bins={custom_channel_bins}," \
                                    f"custom_photon_bins={custom_photon_bins},**{kwargs})"
        self._loaded_spec_data = self.load_info(pha_file, srm_file, srm=srm_custom, channel_bins=custom_channel_bins,
                                                photon_bins=custom_photon_bins)

        self._time_fmt, self._time_scale = "isot", "utc"
        self._start_background_time, self._end_background_time = None, None
        self._start_event_time, self._end_event_time = self._full_obs_time[0], self._full_obs_time[1]

        # used to give the user a warning if incompatible times are set
        self.__warn = True

    @staticmethod
    def getspec(f_pha):
        """ Return all STIX data needed for fitting.

        Only here so that a `CustomLoader` can overwrite for the same `_load1spec` method.

        Parameters
        ----------
        f_pha : str
                String for the STIX spectral file under investigation.

        Returns
        -------
        A 2d array of the channel bin edges (channel_bins), 2d array of the channel bins (channel_bins_inds),
        2d array of the time bins for each spectrum (time_bins), 2d array of livetimes/counts/count rates/count
        rate errors per channel bin and spectrum (lvt/counts/cts_rates/cts_rate_err, respectively).
        """
        return get_file_info.get_stix_file_info(f_pha)

    @staticmethod
    def getsrm(f_srm):
        """ Return all STIX SRM data needed for fitting.

        SRM units returned as counts ph^(-1) cm^(2).

        Only here so that `StixLoader` can overwrite for the same `_load1spec` method.

        Parameters
        ----------
        f_srm : str
                String for the STIX SRM spectral file under investigation.

        Returns
        -------
        A 2d array of the photon and channel bin edges (photon_bins, channel_bins), number of sub-set channels
        in the energy bin (ngrp), starting index of each sub-set of channels (fchan), number of channels in each
        sub-set (nchan), 2d array that is the spectral response (srm).
        """
        return get_file_info.get_stix_srm_file_info(f_srm)

    def load_info(self, f_pha, f_srm, srm=None, channel_bins=None, photon_bins=None):
        """ Loads all the information in for a given spectrum.

        Parameters
        ----------
        f_pha, f_srm : string
                Filenames for the relevant spectral files.

        srm : 2d array
                User defined spectral response matrix. This is accepted over the SRM from f_srm.
                Default: None

        photon_bins, channel_bins: 2d array
                User defined channel bins for the rows and columns of the SRM matrix.
                E.g., custom_channel_bins=[[1,1.5],[1.5,2],...]
                Default: None

        Returns
        -------
        Dictionary of loaded in spectral information in the form {"photon_channel_bins":channel_bins,
                                                                  "photon_channel_mids":np.mean(channel_bins, axis=1),
                                                                  "photon_channel_binning":channel_binning,
                                                                  "count_channel_bins":channel_bins,
                                                                  "count_channel_mids":np.mean(channel_bins, axis=1),
                                                                  "count_channel_binning":channel_binning,
                                                                  "counts":counts,
                                                                  "count_error":count_error,
                                                                  "count_rate":count_rate,
                                                                  "count_rate_error":count_rate_error,
                                                                  "effective_exposure":eff_exp,
                                                                  "srm":srm,
                                                                  "extras":{"pha.file":f_pha,
                                                                            "srm.file":f_srm,
                                                                            "counts=data-bg":False}
                                                                  }.
        """
        # need effective exposure and energy binning since likelihood works on counts, not count rates etc.
        obs_channel_bins, self._channel_bins_inds_perspec, self._time_bins_perspec, self._lvt_perspec, \
            self._counts_perspec, self._counts_err_perspec, self._count_rate_perspec, self._count_rate_error_perspec \
            = self.getspec(f_pha)

        # now calculate the SRM or use a custom one if given
        if srm is None:
            # needs an srm file load it in
            srm_photon_bins, srm_channel_bins, srm, srm_energy_bands = self.getsrm(f_srm)
            # make sure the SRM will only produce counts to match the data
            data_inds2match = np.where((obs_channel_bins[0][:, 0] <= srm_channel_bins[:, 0]) &
                                       (srm_channel_bins[:, 1] <= obs_channel_bins[0][-1, -1]))
            srm = srm[:, data_inds2match[0]]
        else:
            srm_photon_bins = None

        photon_bins = srm_photon_bins if type(photon_bins) == type(None) else photon_bins
        photon_binning = np.diff(photon_bins).flatten()

        # from the srm file #channel_binning = np.diff(channel_bins).flatten()
        channel_bins = obs_channel_bins if type(channel_bins) == type(None) else channel_bins

        # default is no background and all data is the spectrum to be fitted
        self._full_obs_time = [self._time_bins_perspec[0, 0], self._time_bins_perspec[-1, -1]]
        counts = np.sum(self._data_time_select(stime=self._full_obs_time[0], full_data=self._counts_perspec,
                                               etime=self._full_obs_time[1]), axis=0)
        counts_err = np.sqrt(np.sum(self._data_time_select(stime=self._full_obs_time[0],
                                                           full_data=self._counts_err_perspec,
                                                           etime=self._full_obs_time[1]) ** 2, axis=0))
        # sum errors in quadrature, Poisson still sqrt(N)

        _livetimes = np.mean(self._data_time_select(stime=self._full_obs_time[0], full_data=self._lvt_perspec,
                                                    etime=self._full_obs_time[1]), axis=0)
        # to convert a model count rate to counts, so need mean
        eff_exp = np.diff(self._full_obs_time)[0].to_value("s") * _livetimes

        channel_binning = np.diff(obs_channel_bins[0], axis=1).flatten()
        count_rate = counts / eff_exp / channel_binning  # count rates from here are counts/s/keV
        count_rate_error = counts_err / eff_exp / channel_binning  # was np.sqrt(counts)/eff_exp/channel_binning

        # what spectral info you want to know from this observation
        return {"photon_channel_bins": photon_bins,
                "photon_channel_mids": np.mean(photon_bins, axis=1),
                "photon_channel_binning": photon_binning,
                "count_channel_bins": channel_bins,
                "count_channel_mids": np.mean(channel_bins, axis=1),
                "count_channel_binning": channel_binning,
                "counts": counts,
                "count_error": counts_err,
                "count_rate": count_rate,
                "count_rate_error": count_rate_error,
                "effective_exposure": eff_exp,
                "srm": srm,
                "extras": {"pha.file": f_pha,
                           "srm.file": f_srm,
                           "counts=data-bg": False}
                }  # this might make it easier to add different observations together

    def _data_time_select(self, stime, full_data, etime):
        """ Index and return data in time range stime<=data<=etime.

        If stime/etime is None then they are assumed to be the first/last data time.

        Parameters
        ----------
        stime, etime : `astropy.Time`
                The start and end range (inclusive) of the data required.

        full_data : array
                Array to be indexed along axis 0 such that the data returned in within the time range.

        Returns
        -------
        Indexed array.
        """
        if (type(stime) == type(None)) and (type(etime) == type(None)):
            return full_data
        # FIXME: All other cases listed after have errors because of missing reference in full_data[position]
        elif type(stime) == type(None):
            return full_data[np.where(self._time_bins_perspec[:, 1] <= etime)]
        elif type(etime) == type(None):
            return full_data[np.where(stime <= self._time_bins_perspec[:, 0])]
        else:
            return full_data[
                np.where((stime <= self._time_bins_perspec[:, 0]) & (self._time_bins_perspec[:, 1] <= etime))]

    def _channel_bin_info(self, axis):
        """ Returns the old and new channel bins for the indicated axis (count axis, photon axis or both).
        Parameters
        ----------
        axis : string
                Set to "count", "photon", or "photon_and_count" to return the old and new count
                channel bins, photon channel bins, or both.
        Returns
        -------
        Arrays of old_count_bins, new_count_bins, old_photon_bins, new_photon_bins or Nones.
        """
        old_count_bins, new_count_bins, old_photon_bins, new_photon_bins = None, None, None, None
        if (axis == "count") or (axis == "photon_and_count"):
            old_count_bins = self._loaded_spec_data["extras"]["original_count_channel_bins"]
            new_count_bins = self._loaded_spec_data["count_channel_bins"]
        if (axis == "photon") or (axis == "photon_and_count"):
            old_photon_bins = self._loaded_spec_data["extras"]["orignal_photon_channel_bins"]
            new_photon_bins = self._loaded_spec_data["photon_channel_bins"]
        return old_count_bins, new_count_bins, old_photon_bins, new_photon_bins


if __name__ == '__main__':

    StixLoader(pha_file='./data/stx_spectrum_2021feb14_0140_0155.fits',
               srm_file='./data/stx_srm_2021feb14_0140_0155.fits').\
        load_info(f_pha='./data/stx_spectrum_2021feb14_0140_0155.fits',
                  f_srm='./data/stx_srm_2021feb14_0140_0155.fits')
