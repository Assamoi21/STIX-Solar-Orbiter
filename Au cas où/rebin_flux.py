import numpy as np
import numpy.typing
import get_file_info as get
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model

'''
Perform a flux- (i.e. area-) conserving rebinning of a histogram,
given its edges and a set of new edges.
The __name__ ... section has an example using a power law.
W Setterberg Jan 2023
'''


def flux_conserving_rebin(
        old_edges: np.typing.ArrayLike,
        old_values: np.typing.ArrayLike,
        new_edges: np.typing.ArrayLike,
) -> np.ndarray:
    """
    Rebin a histogram by performing a flux-conserving rebinning.
    The total area of the histogram is conserved.
    Adjacent bins are proportionally interpolated for new edges that do not line up.
    Don't make the new edges too finely spaced;
    don't make a new bin fall inside of an old one completely.
    """
    old_edges = np.array(np.sort(old_edges))
    new_edges = np.array(np.sort(new_edges))
    nd = np.diff(new_edges)
    od = np.diff(old_edges)

    if (new_edges[0] < old_edges[0]) or (new_edges[-1] > old_edges[-1]):
        raise ValueError('New edges cannot fall outside range of old edges.')

    if np.all(nd == od):
        return np.ndarray(old_values)

    orig_flux = od * old_values
    ret = np.zeros(new_edges.size - 1)
    for i in range(ret.size):
        ret[i] = interpolate_new_bin(
            original_area=orig_flux,
            old_edges=old_edges,
            new_left=new_edges[i],
            new_right=new_edges[i + 1]
        )

    return ret


def proportional_interp_single_bin(
        left_edge: float,
        right_edge: float,
        interp: float
) -> tuple[float, float]:
    """
    say what portion of a histogram bin belongs on the left and right
    of an edge to interpolate.
    """
    denom = right_edge - left_edge
    right_portion = (right_edge - interp) / denom
    left_portion = (interp - left_edge) / denom
    return left_portion, right_portion


def bounding_interpolate_indices(
        old_edges: np.ndarray,
        left: float,
        right: float
) -> tuple[int, int]:
    """
    find the indices of the old edges that bound the new left
    and right edges.
    """
    indices = np.arange(old_edges.size)
    new_left = indices[old_edges <= left][-1]
    new_right = indices[old_edges >= right][0]
    return new_left, new_right


def interpolate_new_bin(
        original_area: np.array,
        old_edges: np.array,
        new_left: float,
        new_right: float
) -> float:
    """
    interpolate the new bin value given old edges, new edges,
    and the old flux (aka area).
    """
    oa = original_area
    oe = old_edges

    old_start_idx, old_end_idx = bounding_interpolate_indices(
        oe, new_left, new_right
    )

    # portion of edge bins that get grouped with the new bin
    _, left_partial_prop = proportional_interp_single_bin(
        left_edge=oe[old_start_idx],
        right_edge=oe[old_start_idx + 1],
        interp=new_left
    )
    left_partial_area = left_partial_prop * oa[old_start_idx]

    right_partial_prop, _ = proportional_interp_single_bin(
        left_edge=oe[old_end_idx - 1],
        right_edge=oe[old_end_idx],
        interp=new_right
    )
    right_partial_area = right_partial_prop * oa[old_end_idx - 1]

    partial_slice = slice(old_start_idx + 1, old_end_idx - 1)
    have_bad_slice = (partial_slice.start > partial_slice.stop)
    if have_bad_slice:
        raise ValueError('Your new bins are too fine. Use coarser bins.')

    between_area = oa[partial_slice].sum()

    delta = new_right - new_left
    new_bin_value = (left_partial_area + between_area + right_partial_area) / delta
    return new_bin_value


@custom_model
def fitting_photons(x, x0=1., idx=1.):
    """
    :param x:
    :param idx:
    :param x0:
    :return rebinned_counts:
    """
    srm_file = './stx_srm_2021feb14_0140_0155.fits'
    fit_method = "Power Law"
    # fit = fitting.LevMarLSQFitter()
    photon_bins, channel_bins, srm, energy_bands = get.get_srm(srm_file)

    # maxiter = 10 ** 5

    old_edges = np.linspace(0, 150, num=np.size(x) + 1)     # num=32 + 1
    new_edges = np.linspace(0, 150, num=energy_bands + 1)   # num=20 + 1

    mids = old_edges[:-1] + np.diff(old_edges)/2

    """
    x0 = 10
    idx = 3  # FIXME: If power is negative, fitting cannot be done.
    """

    # model = power_law(x0=x0, idx=idx)
    old_values = (x / x0) ** idx
    new_values = flux_conserving_rebin(
        old_edges=old_edges,
        old_values=old_values,
        new_edges=new_edges
    )
    # print(type(new_values))

    # =========================================== Multiplying the SRM ===========================================
    new_counts = np.matrix(new_values) * srm.T
    # print(type(new_counts))
    # print("New counts : ", np.size(new_counts), new_counts)
    # print(new_counts)

    # print(old_edges)
    # print(old_values)
    # print(new_edges)
    # print(new_values)
    # We need to retransform the matrix "new_counts" in a list(?) to avoid errors
    """
    old_edges = np.linspace(0, 150, num=np.size(new_counts) + 1)    # num=1021
    new_edges = np.linspace(0, 150, num=np.size(x) + 1)             # num=21
    final_rebin = flux_conserving_rebin(
        old_edges=old_edges,
        old_values=new_counts,
        new_edges=new_edges
    )
    """"""
    rebinned_counts = np.zeros(np.size(final_rebin))
    for eb in range(np.size(final_rebin)):  # i in range(20)
        rebinned_counts[eb] = float(new_counts[eb])
    """"""
    rebinned_counts = []
    for i in range(20):  # Need to change the value "20" for a variable
        rebinned_counts.append(float(new_counts[0, i]))
    print("rebinned_counts : ", rebinned_counts)
    """

    rebinned_counts = np.zeros(32)
    for eb in range(32):  # i in range(20)
        rebinned_counts[eb] = float(new_counts[0, eb])

    weight_new_counts = rebinned_counts
    if rebinned_counts is not None:
        for i in range(len(rebinned_counts)):
            if rebinned_counts[i] == 0:
                weight_new_counts[i] = 1
            else:
                weight_new_counts[i] = 1 / rebinned_counts[i]

    """
    print(type(rebinned_counts))
    rebin_final = np.zeros(np.size(rebinned_counts))
    for value in range(len(rebinned_counts)):
        rebin_final[value] = rebinned_counts[value]
    print(type(rebin_final))
    """

    if __name__ == '__main__':

        original_flux = compute_binned_flux(old_edges, old_values)
        new_flux = compute_binned_flux(new_edges, new_values)
        print("New flux : ", new_flux)

        print('new_values: ', np.size(new_values), new_values)
        print(f'original flux:     {original_flux:.4f}')
        print(f'rebinned flux:     {new_flux:.4f}')
        print(f'amount off:        {np.abs(original_flux - new_flux):.2e}')
        print(f'machine precision: {np.finfo(float).eps:.2e}')
        print(';)')

        fig, ax = plt.subplots(figsize=(8, 6), layout='constrained')
        ax.stairs(mids, old_edges, label='original')
        ax.stairs(old_values(old_edges[:-1]), old_edges, label='old power law')
        ax.stairs(new_values, new_edges, label='rebinned')
        ax.stairs(rebinned_counts, new_edges, label='rebinned with srm')
        ax.set(
            xscale='log',
            yscale='log',
            xlabel='$x$',
            ylabel=str(fit_method),
            title='Rebin linear-spaced bins to log-spaced bins'
        )
        ax.set_xlim(4, 150)
        ax.set_ylim(10, 10**6)
        ax.legend()
        plt.show()

    return rebinned_counts  # , weight_new_counts, model


"""
@custom_model
def power_law(x, x0=1., idx=1.):
    return (x / x0) ** idx
"""

if __name__ == '__main__':

    def compute_binned_flux(edges, values):
        return np.sum(values * np.diff(edges))

    fit = fitting.LevMarLSQFitter()

    stix_srm_file = './stx_srm_2021feb14_0140_0155.fits'
    fitting_method = "Power Law"

    rebinned_counts_test = fitting_photons(x0=1., idx=1.)

    pl_fit = fit(rebinned_counts_test, x[1:21], y[1:21], maxiter=10**5)
