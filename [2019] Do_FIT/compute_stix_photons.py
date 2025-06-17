import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import nnls

# === 1. Charger les données STIX ===
data_file = "./solo_L1A_stix-sci-spectrogram-2102140001_20210214T014006-20210214T015515_008648_V01.fits"
srm_file = "./stx_srm_2021feb14_0140_0155.fits"

# === 2. Ouvrir les fichiers FITS ===
hdul_data = fits.open(data_file)
hdul_srm = fits.open(srm_file)

# === 3. Extraire les comptes ===
counts_all = hdul_data[2].data['counts']  # shape: (temps, canaux)
timedel = hdul_data[2].data['timedel']
counts = np.mean(counts_all, axis=0)  # Moyenne temporelle

# === 4. Extraire les bandes d'énergie (E_LOW, E_HIGH) détectées ===
e_low = hdul_data[3].data['e_low']
e_high = hdul_data[3].data['e_high']
e_meas = 0.5 * (e_low + e_high)

# === 5. Extraire la SRM et les énergies vraies ===
srm_matrix = np.array(hdul_srm[1].data['MATRIX'])  # shape: (nb_energies_reelles, nb_detectées)
e_true_lo = hdul_srm[1].data['ENERG_LO'].flatten()
e_true_hi = hdul_srm[1].data['ENERG_HI'].flatten()
e_true = 0.5 * (e_true_lo + e_true_hi)

# Transposer pour avoir (nb_mesurées, nb_vraies)
if srm_matrix.shape[0] < srm_matrix.shape[1]:
    srm_matrix = srm_matrix.T

print("SRM shape:", srm_matrix.shape)

# === 6. Adapter dimensions si besoin ===
n_channels = min(srm_matrix.shape[1], len(counts))
counts = counts[:n_channels]
srm_matrix = srm_matrix[:, :n_channels]

# === 7. Résolution : inversion pour obtenir le flux de photons ===
# On résout : C = R · F  <=>  F ≈ nnls(R.T, C)
photons, _ = nnls(srm_matrix.T, counts)

# === 8. Interpolation pour adapter counts aux énergies mesurées ===
# Interpoler counts pour qu'ils aient la même dimension que e_meas
counts_interpolated = np.interp(e_meas, e_true, photons)

# === Adapter e_meas à la taille de counts ===
e_meas = e_meas[:n_channels]

# === 9. Tracer les résultats ===
fig, ax1 = plt.subplots(figsize=(8, 5))

# Tracer les flux de photons estimés (interpolés)
ax1.step(e_true, photons, where='mid', color='darkblue', label='Flux de photons estimé')
ax1.set_xlabel("Énergie réelle (keV)")
ax1.set_ylabel("Flux de photons [ph / s / cm² / keV]", color='darkblue')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.tick_params(axis='y', labelcolor='darkblue')

# Ajouter un second axe y pour les comptes mesurés
ax2 = ax1.twinx()
ax2.plot(e_meas, counts, color='red', label='Comptes mesurés', linestyle='--')
ax2.set_ylabel("Comptes mesurés [counts]", color='red')
ax2.tick_params(axis='y', labelcolor='red')

# Ajouter une légende et un titre
fig.tight_layout()
plt.title("Comparaison entre les flux de photons et les comptes mesurés (STIX)")
ax1.grid(True, which='both', linestyle='--', alpha=0.5)
fig.legend(loc="upper right")
plt.show()

# === 10. Clean-up ===
hdul_data.close()
hdul_srm.close()

