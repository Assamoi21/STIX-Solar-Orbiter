

# --- Charger les fichiers FITS ---
path_counts = "C:/Users/Assamoua Koman/Downloads/STIX_2025/solo_L1A_stix-sci-spectrogram-2102140001_20210214T014006-20210214T015515_008648_V01.fits"
path_srm = "C:/Users/Assamoua Koman/Downloads/STIX_2025/stx_srm_2021feb14_0140_0155.fits"

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling import FittableModel, Parameter
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.constants import k_B
import astropy.units as u


path_counts = "C:/Users/Assamoua Koman/Downloads/STIX_2025/solo_L1A_stix-sci-spectrogram-2102140001_20210214T014006-20210214T015515_008648_V01.fits"
path_srm = "C:/Users/Assamoua Koman/Downloads/STIX_2025/stx_srm_2021feb14_0140_0155.fits"

# --- Charger les fichiers FITS ---
counts_file = fits.open(path_counts)
srm_file = fits.open(path_srm)

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

# --- Lire les données de comptage ---
counts_data = extract_stix_data(counts_file)

# counts_data = counts_file[2].data
# energy_bounds_counts = counts_file[3].data

counts_err_all = np.mean(counts_data['counts_err'], axis=0)
exposure = np.mean(counts_data['timedel'])
times = counts_data['time']

# times = np.mean(counts_data.time)

print("counts_data", counts_data['counts'].shape)
print("time :", times.shape)
background = np.mean(counts_data['counts'][55:63, :], axis=0)
print("background", background)

new_data = counts_data['counts'] - background  
new_data[new_data < 0] = 0  

# counts_all = np.mean(counts_data['counts'], axis=0)
counts_all = np.mean(new_data, axis=0)

e_low_det_all = counts_data['e_low']
e_high_det_all = counts_data['e_high']

# --- Lire la réponse instrumentale ---

srm_data = extract_srm_data(srm_file)
e_low_true = srm_data['ENERG_LO']
e_high_true = srm_data['ENERG_HI']
matrix = srm_data['MATRIX']

# --- Appliquer le masque canal utilisable ---

usable_channels = np.arange(min(matrix.shape[1], len(e_low_det_all)))

counts = counts_all[usable_channels]
counts_err = counts_err_all[usable_channels]
e_low_det = e_low_det_all[usable_channels]
e_high_det = e_high_det_all[usable_channels]
print("counts : ", counts.shape)
# --- Fonction d'intégration du flux photonique ---

def integrate_flux(e1, e2, model_func, n_points=10):
    energies = np.linspace(e1, e2, n_points)
    fluxes = model_func(energies)
    return np.trapezoid(fluxes, energies) / (e2 - e1)

# --- Modèle personnalisé compatible avec LevMarLSQFitter ---

class ForwardFoldedPowerLaw(FittableModel):
    n_inputs = 1
    n_outputs = 1

    amplitude = Parameter(default=1e-2)
    alpha = Parameter(default=2.0)
    x_0 = 100.0  # énergie pivot en keV, fixe ici

    def __init__(self, e_low_true, e_high_true, matrix, exposure, **kwargs):
        super().__init__(**kwargs)
        self.e_low_true = e_low_true
        self.e_high_true = e_high_true
        self.matrix = matrix
        self.exposure = exposure

    def evaluate(self, x, amplitude, alpha):
        model_func = lambda E: amplitude * (E / self.x_0) ** (-alpha)
        true_fluxes = np.array([
            integrate_flux(e1, e2, model_func)
            for e1, e2 in zip(self.e_low_true, self.e_high_true)
        ])
        folded = np.dot(true_fluxes, self.matrix) / self.exposure
        return folded

class ForwardFoldedGaussienne(FittableModel):
    n_inputs = 1
    n_outputs = 1

    amp = Parameter(default=1e-3)
    mean = Parameter(default=15.0)
    stddev = Parameter(default=1.0)

    def __init__(self, e_low_true, e_high_true, matrix, exposure, **kwargs):
        super().__init__(**kwargs)
        self.e_low_true = e_low_true
        self.e_high_true = e_high_true
        self.matrix = matrix
        self.exposure = exposure

    def evaluate(self, x, amp, mean, stddev):
        model_func = lambda E: amp * np.exp(-0.5 * ((E - mean) / stddev) ** 2)

        # Calculer le flux photonique forward-foldé
        true_fluxes = np.array([
            integrate_flux(e1, e2, model_func)
            for e1, e2 in zip(self.e_low_true, self.e_high_true)
        ])
        folded = np.dot(true_fluxes, self.matrix) / self.exposure
        return folded

class ForwardFoldedVTH(FittableModel):
    n_inputs = 1
    n_outputs = 1

    EM = Parameter(default=1e48)        # Emission Measure (cm^-3)
    T = Parameter(default=1e7)          # Temperature (K)

    def __init__(self, e_low_true, e_high_true, matrix, exposure, **kwargs):
        super().__init__(**kwargs)
        self.e_low_true = e_low_true
        self.e_high_true = e_high_true
        self.matrix = matrix
        self.exposure = exposure

    def evaluate(self, x, EM, T):
        # Constantes
        kT_keV = (k_B.value * T) / 1.60218e-16 / 1e3  # Convert to keV

        def vth_model(E_keV):
            return EM * E_keV**-1 * np.exp(-E_keV / kT_keV)

        true_fluxes = np.array([
            integrate_flux(e1, e2, vth_model)
            for e1, e2 in zip(self.e_low_true, self.e_high_true)
        ])

        folded = np.dot(true_fluxes, self.matrix) / self.exposure
        return folded

# --- Nettoyage des données avant fit ---

# Supprimer les canaux où les erreurs sont nulles ou non-finies
valid = (counts_err > 0) & np.isfinite(counts_err) & np.isfinite(counts)

counts = counts[valid]
counts_err = counts_err[valid]
x_fake = np.zeros_like(counts)  # x doit avoir la même taille

matrix = matrix[:, valid]       # adapter la SRM
e_low_det = e_low_det[valid]
e_high_det = e_high_det[valid]

edges_det = np.append(e_low_det, e_high_det[-1])
dE_det = np.diff(edges_det)

Edges_photon = np.append(e_low_true, e_high_true[-1])

# --- Fitting avec LevMarLSQFitter ---



# === Intervalle pour le FIT
fit_Emin = 10.0  # keV
fit_Emax = 20.0  # keV

# Centre des canaux détecteurs
E_center_det = 0.5 * (e_low_det + e_high_det)

# Masque des canaux utilisés pour le fitting
fit_mask = (edges_det[:-1] >= fit_Emin) & (edges_det[1:] <= fit_Emax)

# Sous-ensembles pour le fitting
x_fit = x_fake[fit_mask]
counts_fit = counts[fit_mask]
counts_err_fit = counts_err[fit_mask]
matrix_fit = matrix[:, fit_mask]

# Créer un modèle avec SRM réduite pour fitting
model_fit = ForwardFoldedPowerLaw(e_low_true, e_high_true, matrix_fit, exposure)

# Fitting
fitter = LevMarLSQFitter()
fitted_model = fitter(model_fit, x_fit, counts_fit / exposure,
                      weights=1.0 / (counts_err_fit / exposure))

# Paramètres récupérés
amplitude = fitted_model.amplitude.value
alpha = fitted_model.alpha.value
x_0 = fitted_model.x_0

# Modèle complet pour affichage sur tout le domaine
model_display = ForwardFoldedPowerLaw(e_low_true, e_high_true, matrix, exposure)
model_display.amplitude = fitted_model.amplitude
model_display.alpha = fitted_model.alpha

# Calcul du modèle simulé complet
rate_modeled_full = model_display(x_fake)

plt.figure()
plt.step(edges_det[:-1], (rate_modeled_full / dE_det), where='post',
         label='Fitted Model (full)', color='blue')
plt.step(edges_det[:-1], (counts / exposure) / dE_det, where='post',
         label='STIX Data', color='red')
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Channel Energy (keV)")
plt.ylabel("Count rate [counts / (s keV)]")
plt.title(f"Fitting effectué sur [{fit_Emin}, {fit_Emax}] keV")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.text(0.05, 0.4,
         f"Power Law:\n amplitude = {amplitude:.2e}\n alpha = {alpha:.2f}",
         transform=plt.gca().transAxes,
         fontsize=10,
         verticalalignment='top',
         bbox=dict(facecolor='white', alpha=0.7))
plt.tight_layout()


# --- Calculer le flux photonique incident ---
model_func = lambda E: amplitude * (E / x_0)**(-alpha)
flux_photons = np.array([
    integrate_flux(e1, e2, model_func)
    for e1, e2 in zip(e_low_true, e_high_true)
])

plt.figure()
plt.step(Edges_photon[:-1], flux_photons, where='mid',
         label='Photon model (before SRM)', color='green')
# plt.step(E_center, flux_photons, where='mid', label='Photon model (before SRM)', color='green')
plt.xscale('log')
plt.yscale('log')
plt.xlabel("True Energy (keV)")
plt.ylabel("Photon flux [photons / (s cm² keV)]")
plt.title("True Photon Flux Model")
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
plt.text(0.05, 0.4,
         f"Power Law:\n amplitude = {amplitude:.2e}\n alpha = {alpha:.2f} \n",
         transform=plt.gca().transAxes,
         fontsize=10,
         verticalalignment='top',
         bbox=dict(facecolor='white', alpha=0.7))
plt.tight_layout()


# plt.show()





# # Créer le modèle avec la SRM et exposition
# model = ForwardFoldedPowerLaw(e_low_true, e_high_true, matrix, exposure)

# # Créer le fitter
# fitter = LevMarLSQFitter()

# # Fit : x est factice, juste pour respecter l'API Astropy
# fitted_model = fitter(model, x_fake, counts / exposure, weights=1.0 / (counts_err / exposure))

# #Récupérer les résultats

# amplitude = fitted_model.amplitude.value
# alpha = fitted_model.alpha.value
# x_0 = fitted_model.x_0

# # --- Calcul du modèle ajusté ---

# rate_modeled = fitted_model(x_fake)


# # --- Calculer le flux photonique incident ---
# model_func = lambda E: amplitude * (E / x_0)**(-alpha)
# # E_center = 0.5 * (e_low_true + e_high_true)
# flux_photons = np.array([
#     integrate_flux(e1, e2, model_func)
#     for e1, e2 in zip(e_low_true, e_high_true)
# ])

# # Choose my interval for plotting
# E_plot_min = 10.0
# E_plot_max = 15.0
# mask_flux = (Edges_photon[:-1] >= E_plot_min) & (Edges_photon[1:] <= E_plot_max)
# mask_det = (edges_det[:-1] >= E_plot_min) & (edges_det[1:] <= E_plot_max)

# # --- Plot ---
# plt.figure()
# # plt.step(edges_det[:-1], (counts / exposure) / dE_det, where='post', label='Data', color='red')
# # plt.step(edges_det[:-1], rate_modeled / dE_det, where='post', label='Fitted Model', color='blue')
# plt.step(edges_det[:-1][mask_det], (rate_modeled / dE_det)[mask_det], where='mid',
#          label='Fitted Model', color='blue')
# plt.step(edges_det[:-1][mask_det], ((counts / exposure) / dE_det)[mask_det], where='mid',
#          label='STIX Data', color='red')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel("Channel Energy (keV)")
# plt.ylabel(r"Count rate [counts / (s keV)]")
# plt.title("Forward-folded model fit to STIX data")
# plt.legend()
# plt.grid(True, which="both", ls="--", alpha=0.5)
# # Ajouter les paramètres du modèle sur le graphique
# plt.text(0.05, 0.4,
#          f"Power Law:\n amplitude = {amplitude:.2e}\n alpha = {alpha:.2f} \n",
#          transform=plt.gca().transAxes,
#          fontsize=10,
#          verticalalignment='top',
#          bbox=dict(facecolor='white', alpha=0.7))
# plt.tight_layout()  


# # --- Calculer le flux photonique incident ---
# # Calculer le flux photonique incident (avant SRM)

# plt.figure()
# plt.step(Edges_photon[:-1][mask_flux], flux_photons[mask_flux], where='mid',
#          label='Photon model (before SRM)', color='green')
# # plt.step(E_center, flux_photons, where='mid', label='Photon model (before SRM)', color='green')
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

# # plt.figure()
# # plt.imshow(matrix, origin='lower', aspect='auto')
# # plt.colorbar(label='Probabilité')
# # plt.xlabel("Canaux détecteurs")
# # plt.ylabel("Canaux photons")

# plt.show()

# ''' TEST with count rate and fluxes'''
# # Surface effective (cm²)
# Area = 6.0  # tu peux la rendre paramétrable

# # Dimensions
# n_channels = len(e_low_det)
# # deltaE = e_high_det - e_low_det

# # Comptes par canal
# mean_counts = counts
# mean_counts_err = counts_err

# # Taux de comptage
# rate = mean_counts / exposure
# rate_err = mean_counts_err / exposure

# # Flux (photons / s / cm² / keV)
# flux = rate / (Area * dE_det)
# flux_err = rate_err / (Area * dE_det)


# unit = 'counts'  # ou 'rate' ou 'counts'

# if unit == 'rate':
#     y_data = rate[mask_det]
#     y_err = rate_err[mask_det]
#     model_y = (rate_modeled / dE_det)[mask_det]
#     y_label = "Count Rate [counts / (s keV)]"
# elif unit == 'counts':
#     y_data = mean_counts[mask_det]
#     y_err = mean_counts_err[mask_det]
#     model_y = (rate_modeled * exposure)[mask_det]
#     y_label = "Counts"
# elif unit == 'flux':
#     y_data = flux[mask_det]
#     y_err = flux_err[mask_det]
#     model_y = (rate_modeled / (Area * dE_det))[mask_det]
#     y_label = "Photon Flux [photons / (s cm² keV)]"
# else:
#     raise ValueError("Choisir unit = 'rate', 'counts' ou 'flux'")

# # --- Plot ---
# plt.figure()
# plt.step(edges_det[:-1][mask_det], y_data, where='mid', label=f'Data ({unit})', color='red')
# plt.step(edges_det[:-1][mask_det], model_y, where='mid',
#          label='Fitted Model (after SRM)', color='blue')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel("Channel Energy (keV)")
# plt.ylabel(y_label)
# plt.title(f"STIX Spectrum - {unit.capitalize()}")
# plt.grid(True, which="both", ls="--", alpha=0.5)
# plt.legend()
# plt.tight_layout()
# plt.show()
