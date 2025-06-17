import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

# print("Counts")
# hdulist = fits.open("C:/Users/Assamoua Koman/Downloads/solo_L1_stix-sci-xray-spec_20230319T175504-20230320T000014_V02_2303197888-65462.fits")
# for i, hdu in enumerate(hdulist):
#     print(f"HDU {i}:")
#     if hasattr(hdu, 'data') and hdu.data is not None:
#         print("  Columns:", hdu.data.dtype.names)
#     else:
#         print("  No data")

test1 = fits.open("C:/Users/Assamoua Koman/Downloads/STIX_2025/solo_L1A_stix-sci-spectrogram-2102140001_20210214T014006-20210214T015515_008648_V01.fits")
test2 = fits.open("C:/Users/Assamoua Koman/Downloads/STIX_2025/stx_srm_2021feb14_0140_0155.fits")

test3 = fits.open("C:/Users/Assamoua Koman/Downloads/solo_L1_stix-sci-xray-spec_20230319T175504-20230320T000014_V02_2303197888-65462.fits")
test4 = fits.open("C:/Users/Assamoua Koman/Downloads/stx_srm_2303197888.fits")

# hdul = test1
# header5 = hdul[0].header[5] 
# print("Header:", test1[0].header)
# print("Header____________________________________________________:")

def extract_stix_header(hdulist):
    result = {}
    for key, value, comment in hdulist[0].header.cards:
        result[key] = value
    for key, value, comment in hdulist[3].header.cards:
        result[key] = value
    return result

def extract_energie_header(hdulist):
    result = {}
    for key, value, comment in hdulist[3].header.cards:
        result[key] = value
    return result

v = extract_energie_header(test3)
print("Header data:", v)
print("veleur chercher :", extract_energie_header(test1))
print('energie test1 :', test1[2].data['counts'].shape)
print('energie test3 :', test3[2].data['counts'].shape[1])

# header_data = []
# for key, value, comment in test3[0].header.cards:
#         header_data.append((key, value, comment))
# # Créer un tableau pour afficher le header de manière organisée
# header_data = []
# for key, value, comment in test3[0].header.cards:
#     header_data.append((key, value, comment))

# # Convertir en tableau Astropy pour une belle visualisation
# header_table = Table(rows=header_data, names=['Keyword', 'Value', 'Comment'])

# # Afficher le tableau avec toutes les lignes
# header_table.pprint(max_lines=-1, max_width=-1)

# # Alternative avec pandas pour un affichage encore plus propre (si vous avez pandas installé)
# try:
#     import pandas as pd
#     pd.set_option('display.max_rows', None)
#     pd.set_option('display.max_columns', None)
#     pd.set_option('display.width', None)
#     pd.set_option('display.max_colwidth', None)
    
#     df = pd.DataFrame(header_data, columns=['Keyword', 'Value', 'Comment'])
#     print("\nAffichage avec pandas:")
#     print(df)
# except ImportError:
#     print("\nPour un affichage encore meilleur, installez pandas: pip install pandas")

# # Fermer le fichier FITS
# test1.close()



def extract_stix_data(hdulist):
    result = {}

    for hdu in hdulist:
        if not hasattr(hdu, 'columns'):
            continue
        colnames = hdu.columns.names

        # Time & Timedel
        if 'time' in colnames and 'timedel' in colnames:
            result['time'] = hdu.data['time']
            result['timedel'] = hdu.data['timedel']

        # Counts
        if 'counts' in colnames:
            result['counts'] = hdu.data['counts']
        if 'counts_comp_err' in colnames or 'counts_err' in colnames:
            err_col = 'counts_comp_err' if 'counts_comp_err' in colnames else 'counts_err'
            result['counts_err'] = hdu.data[err_col]

        # Triggers (optionnel)
        if 'triggers' in colnames:
            result['triggers'] = hdu.data['triggers']

        # Energy bins
        if 'e_low' in colnames and 'e_high' in colnames:
            result['e_low'] = hdu.data['e_low']
            result['e_high'] = hdu.data['e_high']

        # Version info
        if 'obt_start' in colnames and 'obt_end' in colnames:
            result['obt_start'] = hdu.data['obt_start']
            result['obt_end'] = hdu.data['obt_end']

    # Vérifications de base
    required_keys = ['counts', 'counts_err', 'e_low', 'e_high', 'time', 'timedel']
    for key in required_keys:
        if key not in result:
            print(f"⚠️  Attention : {key} non trouvé dans le FITS.")
    return result

def extract_srm_data(hdulist):
    result = {}

    for hdu in hdulist:
        if not hasattr(hdu, 'columns'):
            continue
        colnames = hdu.columns.names

        # Matrix
        if 'MATRIX' in colnames:
            result['MATRIX'] = hdu.data['MATRIX']

        # Energy bins
        if 'ENERG_LO' in colnames and 'ENERG_HI' in colnames:
            result['ENERG_LO'] = hdu.data['ENERG_LO']
            result['ENERG_HI'] = hdu.data['ENERG_HI']


    # Vérifications de base
    required_keys = ['MATRIX', 'ENERG_LO', 'ENERG_HI']
    for key in required_keys:
        if key not in result:
            print(f"⚠️  Attention : {key} non trouvé dans le FITS.")
    return result

data_dict = extract_stix_data(test3)
srm_data = extract_srm_data(test4)

counts = data_dict['counts']
counts_err = data_dict['counts_err']
e_low = data_dict['e_low']
e_high = data_dict['e_high']
time = data_dict['time']
timedel = data_dict['timedel']

print("Counts data shape:", counts.shape)
print("Counts error data shape:", counts_err.shape)
print("Energy low data shape:", e_low.shape)
print("Energy high data shape:", e_high.shape)  
print("Time data shape:", time.shape)
print("Time delta data shape:", timedel.shape)

e_low_true = srm_data['ENERG_LO']
e_high_true = srm_data['ENERG_HI']
matrix = srm_data['MATRIX']

print("SRM data shape:", matrix.shape)
print("Energy low true data shape:", e_low_true.shape)
print("Energy high true data shape:", e_high_true.shape)


# # Charger les données FITS
# hdulist1 = fits.open("C:/Users/Assamoua Koman/Downloads/STIX_2025/solo_L1A_stix-sci-spectrogram-2102140001_20210214T014006-20210214T015515_008648_V01.fits")
# hdulist = fits.open("C:/Users/Assamoua Koman/Downloads/STIX_2025/stx_srm_2021feb14_0140_0155.fits")

# # Extraction des données
# header1 = hdulist1[2].header
# header3 = hdulist[3].header
# data1 = hdulist[2].data
# data2 = hdulist[3].data
# data11 = hdulist1[2].data
# data12 = hdulist1[3].data

# import numpy as np
# import matplotlib.pyplot as plt
# from astropy.io import fits
# print("-----------------------------------------------")
# print("SRM")
# hdulist = fits.open("C:/Users/Assamoua Koman/Downloads/STIX_2025/stx_srm_2021feb14_0140_0155.fits")
# for i, hdu in enumerate(hdulist):
#     print(f"HDU {i}:")
#     if hasattr(hdu, 'data') and hdu.data is not None:
#         print("  Columns:", hdu.data.dtype.names)
#     else:
#         print("  No data")



# import numpy as np
# import matplotlib.pyplot as plt
# from astropy.io import fits
# from astropy.modeling import FittableModel, Parameter
# from astropy.modeling.fitting import LevMarLSQFitter
# from astropy.constants import k_B
# from scipy.signal import savgol_filter  # <- ajout pour le filtrage
# import astropy.units as u

# # === PARAMETRE POUR LE LISSAGE ===
# apply_smoothing = True  # Mettre False pour désactiver
# savgol_window = 5       # Doit être impair et <= taille du vecteur
# savgol_order = 2        # Ordre du polynôme

# # === FICHIERS ===
# path_counts = "C:/Users/Assamoua Koman/Downloads/STIX_2025/solo_L1A_stix-sci-spectrogram-2102140001_20210214T014006-20210214T015515_008648_V01.fits"
# path_srm = "C:/Users/Assamoua Koman/Downloads/STIX_2025/stx_srm_2021feb14_0140_0155.fits"

# counts_file = fits.open(path_counts)
# srm_file = fits.open(path_srm)

# # === DONNEES DE COMPTAGE ===
# counts_data = counts_file[2].data
# energy_bounds_counts = counts_file[3].data
# counts_all = np.mean(counts_data['counts'], axis=0)
# counts_err_all = np.mean(counts_data['counts_err'], axis=0)
# exposure = np.mean(counts_data['timedel'])
# e_low_det_all = energy_bounds_counts['e_low']
# e_high_det_all = energy_bounds_counts['e_high']

# # === REPONSE INSTRUMENTALE ===
# srm_data = srm_file[1].data
# e_low_true = srm_data['ENERG_LO']
# e_high_true = srm_data['ENERG_HI']
# matrix = srm_data['MATRIX']

# usable_channels = np.arange(min(matrix.shape[1], len(e_low_det_all)))

# counts = counts_all[usable_channels]
# counts_err = counts_err_all[usable_channels]
# e_low_det = e_low_det_all[usable_channels]
# e_high_det = e_high_det_all[usable_channels]

# # === LISSAGE ===
# if apply_smoothing and len(counts) >= savgol_window:
#     counts = savgol_filter(counts, window_length=savgol_window, polyorder=savgol_order)

# # === FONCTION D'INTEGRATION ===
# def integrate_flux(e1, e2, model_func, n_points=10):
#     energies = np.linspace(e1, e2, n_points)
#     fluxes = model_func(energies)
#     return np.trapezoid(fluxes, energies) / (e2 - e1)

# # === MODELE FORWARDFOLDED POWERLAW ===
# class ForwardFoldedPowerLaw(FittableModel):
#     n_inputs = 1
#     n_outputs = 1
#     amplitude = Parameter(default=1e-2)
#     alpha = Parameter(default=2.0)
#     x_0 = 100.0

#     def __init__(self, e_low_true, e_high_true, matrix, exposure, **kwargs):
#         super().__init__(**kwargs)
#         self.e_low_true = e_low_true
#         self.e_high_true = e_high_true
#         self.matrix = matrix
#         self.exposure = exposure

#     def evaluate(self, x, amplitude, alpha):
#         model_func = lambda E: amplitude * (E / self.x_0) ** (-alpha)
#         true_fluxes = np.array([
#             integrate_flux(e1, e2, model_func)
#             for e1, e2 in zip(self.e_low_true, self.e_high_true)
#         ])
#         folded = np.dot(true_fluxes, self.matrix) / self.exposure
#         return folded

# # === NETTOYAGE DES DONNEES ===
# valid = (counts_err > 0) & np.isfinite(counts_err) & np.isfinite(counts)
# counts = counts[valid]
# counts_err = counts_err[valid]
# x_fake = np.zeros_like(counts)
# matrix = matrix[:, valid]
# e_low_det = e_low_det[valid]
# e_high_det = e_high_det[valid]

# edges_det = np.append(e_low_det, e_high_det[-1])
# dE_det = np.diff(edges_det)
# Edges_photon = np.append(e_low_true, e_high_true[-1])

# # === FITTING ===
# model = ForwardFoldedPowerLaw(e_low_true, e_high_true, matrix, exposure)
# fitter = LevMarLSQFitter()
# fitted_model = fitter(model, x_fake, counts / exposure, weights=1.0 / (counts_err / exposure))

# amplitude = fitted_model.amplitude.value
# alpha = fitted_model.alpha.value
# x_0 = fitted_model.x_0
# rate_modeled = fitted_model(x_fake)

# # === FLUX PHOTONIQUE INCIDENT ===
# model_func = lambda E: amplitude * (E / x_0)**(-alpha)
# flux_photons = np.array([
#     integrate_flux(e1, e2, model_func)
#     for e1, e2 in zip(e_low_true, e_high_true)
# ])

# # === PLOTS ===
# E_plot_min = 4.0
# E_plot_max = 40.0
# mask_flux = (Edges_photon[:-1] >= E_plot_min) & (Edges_photon[1:] <= E_plot_max)
# mask_det = (edges_det[:-1] >= E_plot_min) & (edges_det[1:] <= E_plot_max)

# plt.figure()
# plt.step(edges_det[:-1][mask_det], (rate_modeled / dE_det)[mask_det], where='mid',
#          label='Fitted Model', color='blue')
# plt.step(edges_det[:-1][mask_det], ((counts / exposure) / dE_det)[mask_det], where='mid',
#          label='STIX Data (smoothed)' if apply_smoothing else 'STIX Data', color='red')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel("Channel Energy (keV)")
# plt.ylabel(r"Count rate [counts / (s keV)]")
# plt.title("Forward-folded model fit to STIX data")
# plt.legend()
# plt.grid(True, which="both", ls="--", alpha=0.5)
# plt.text(0.05, 0.4,
#          f"Power Law:\n amplitude = {amplitude:.2e}\n alpha = {alpha:.2f} \n",
#          transform=plt.gca().transAxes,
#          fontsize=10,
#          verticalalignment='top',
#          bbox=dict(facecolor='white', alpha=0.7))
# plt.tight_layout()

# plt.figure()
# plt.step(Edges_photon[:-1][mask_flux], flux_photons[mask_flux], where='mid',
#          label='Photon model (before SRM)', color='green')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel("True Energy (keV)")
# plt.ylabel("Photon flux [photons / (s cm² keV)]")
# plt.title("True Photon Flux Model")
# plt.grid(True, which="both", ls="--", alpha=0.5)
# plt.legend()
# plt.text(0.05, 0.4,
#          f"Power Law:\n amplitude = {amplitude:.2e}\n alpha = {alpha:.2f} \n",
#          transform=plt.gca().transAxes,
#          fontsize=10,
#          verticalalignment='top',
#          bbox=dict(facecolor='white', alpha=0.7))
# plt.tight_layout()
# plt.show()
