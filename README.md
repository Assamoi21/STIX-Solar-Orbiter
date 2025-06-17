# STIX-Solar-Orbiter

# STIX Spectrum Fitting Tool

Ce projet permet d'ajuster les spectres mesurés par l'instrument STIX (Solar Orbiter) en utilisant des modèles spectraux personnalisés et la réponse instrumentale (SRM). Il fournit une interface interactive basée sur Tkinter pour visualiser, modéliser et comparer les données.

## 📦 Fonctionnalités

- Chargement de fichiers STIX FITS.
- Choix de modèles spectraux : PowerLaw1D, BrokenPowerLaw1D, V_TH, etc.
- Forward folding via la matrice SRM.
- Ajustement automatique des modèles.
- Visualisation interactive des résultats (flux, taux, counts).
- Support d’erreurs statistiques sur les données.

## 🖥️ Interface

L’interface Tkinter permet de :
- Charger les fichiers FITS (spectres et SRM).
- Sélectionner des intervalles d’énergie.
- Ajouter et configurer des modèles.
- Lancer le fitting et afficher les résultats.


## ⚙️ Installation

### Prérequis

- Python 3.9+
- Astropy
- Numpy
- Scipy
- Matplotlib
- Tkinter

### Installation

git clone https://github.com/Assamoi21/STIX-Solar-Orbiter.git
cd STIX-Solar-Orbiter
pip install -r requirements.txt


## 📁 Données

Des exemples de fichiers FITS (spectres et SRM) sont inclus dans le dépôt. Vous pouvez les télécharger eventuellement depuis le site officiel.

## 📜 Licence

Ce projet est sous licence MIT.

## 👨‍🔬 Auteurs

    Abdallah Hamini, Assamoua Koman

    Contact : abdallah.hamini@obspm.fr


