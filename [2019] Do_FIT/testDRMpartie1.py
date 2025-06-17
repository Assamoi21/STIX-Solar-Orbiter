import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling.models import custom_model
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit
import scipy.io
import re

## Open data spectrum information

hdulist1 = fits.open("C:/Users/stage.PCGIL52/Desktop/OSPEXinPythonBis/OSPEXinPython/Data/hsi_spectrum_20020220_080000.fits")
header11 = hdulist1[1].header
header13 = hdulist1[3].header
data11 = hdulist1[1].data
data12 = hdulist1[2].data
Rate = data11.RATE
Livetime = data11.LIVETIME
print('livetime', Livetime.shape)
Time = data11.TIME - 2
Time_del = data11.TIMEDEL
E_min = data12.E_MIN
E_max = data12.E_MAX
print('E_min:',E_min)
print('E_max:',E_max)
Area = header13[24]
E_mean = np.mean(E_min)
Accum_Time = np.sum(Time_del)
n = len(E_min)

#print('Rate:',Rate)
#print('Rate.shape: ', Rate.shape)

deltaE = np.zeros(shape=(n))
for i in range(n):
    deltaE[i] = E_max[i] - E_min[i]

#Rate total
CountRate = np.zeros(shape=(n))
for i in range(n):
    CountRate[i] = np.mean(Rate[:, i])

#Counts total
Counts = np.zeros(shape=(n))
for i in range(n):
    Counts[i] = np.mean(Rate[:,i]*Accum_Time)

#Flux total
Flux = np.zeros(shape=(n))
for i in range(n):
    Flux[i] = np.mean(Rate[:,i]) / (Area*deltaE[i]-2) # Why -2 ?

print('Flux: ', Flux.shape)

### Selection of energy band : We will separate all energy band in three parts : E_min[0]-E_min, E_min-E_max
# and Emax-E_min[-1]


energy_band = input("Entrez la bande d'energie de la forme 10-20 : ")
print("Energy band :", energy_band)
print(type(energy_band))

def returnbands(energy_band,E_min,E_max):
    """Take in account different energy band for plotting """
    regex = re.compile('-')
    regex.findall(energy_band)
    print(regex.findall(energy_band))
    if regex.findall(energy_band) ==[]:
        print('Energy bands from 4 to 100 keV')
        E_min_index_min=0
        E_min_index_max=len(E_min)
        E_max_index_min=0
        E_max_index_max=len(E_max)
    elif regex.findall(energy_band) !=[]:
        print('Interval energy')
        energy_band_min=energy_band.split('-')[0]
        energy_band_max=energy_band.split('-')[1]
        print('energy_band_min:',energy_band_min)
        print('energy_band_max:',energy_band_max)
        E_min_index_min = (np.where(E_min == int(energy_band_min)))[0][0]
        E_min_index_max = (np.where(E_min == int(energy_band_max)))[0][0]
        E_max_index_min = (np.where(E_max == int(energy_band_min)))[0][0]
        E_max_index_max = (np.where(E_max == int(energy_band_max)))[0][0]
        #print('E_min_index_min:',E_min_index_min)
        #print('E_min_index_max:',E_min_index_max)
        #print('E_max_index_min:',E_max_index_min)
        #print('E_min_index_max:',E_min_index_max)
    return E_min_index_min,E_min_index_max,E_max_index_min


E_min_index_min,E_min_index_max,E_max_index_min=returnbands(energy_band,E_min,E_max)
print('E_min_index_min,E_min_index_max,E_max_index_min:',E_min_index_min,E_min_index_max,E_max_index_min)
## Need to improve this funtion to take care if the value gaved by the user is not exact (do it after)

## Open instrument response data information
hdulist = fits.open("C:/Users/stage.PCGIL52/Desktop/OSPEXinPythonBis/OSPEXinPython/Data/hsi_srm_20020220_080000.fits" )
hdulist.info()

header1 = hdulist[1].header
header3 = hdulist[3].header
data1 = hdulist[1].data
data2 = hdulist[2].data
Energy_low = data1.ENERG_LO
print('E_min: ', E_min)
print('E_low: ', Energy_low)
# print("energy law shape: ", Energy_low.shape)
Energy_hi = data1.ENERG_HI
print('E_max: ', E_max)
print('E_high: ',Energy_hi)
n_grp = data1.N_GRP
#print('N_grp: ', n_grp)
f_chan = data1.F_CHAN
#print('F_chan: ', f_chan)
n_chan = data1.N_CHAN
#print('N_chan: ', n_chan)

matrix = data1.MATRIX


# The three parts :
if E_min_index_max==len(E_min):
    edges=np.append(E_min,E_max[-1])
else:
    edges = E_min[E_min_index_min:E_min_index_max+2]
print('edges:',edges)
edges0=E_min[:E_min_index_min+1]
#print('edges0:',edges0)
#edges1=np.append(E_min[E_min_index_max+1:],E_max[-1])
#print('edges1:',edges1)

# The energy bands total
edges_bis = np.append(E_min,E_max[-1])
print('edges_bis:',edges_bis)

# the width of the channels
dE = np.diff(edges)
dE0 = np.diff(edges0)
#dE1 = np.diff(edges1)
dE_bis=np.diff(edges_bis)
#print('dE:',dE)
#print('dE0:',dE0)
#print('dE1:',dE1)
#print('dE_bis:',dE_bis)
nbis=len(E_min[E_min_index_min:E_min_index_max+1])
nbis0=len(E_min[:E_min_index_min])
#nbis1=len(E_min[E_min_index_max+1:])
print('nbis:',nbis)
#print('nbis0:',nbis0)
#print('nbis1:',nbis)

#Rate
PhotonRate = np.zeros(shape=(nbis))
for i in range(nbis):
    PhotonRate[i] = Rate[0, E_min_index_min+i]
    #print('PhotonRate[i]:',PhotonRate[i])
dNc_dtdE = PhotonRate


PhotonRate0 = np.zeros(shape=(nbis0))
for i in range(nbis0):
   PhotonRate0[i] = Rate[0,i]
    #print('PhotonRate[i]:',PhotonRate[i])
dNc_dtdE0 = PhotonRate0
#dNc_dtdE_model = PhotonRate/(np.diag(matrix)[:E_min_index_min])
#PhotonRate1 = np.zeros(shape=(nbis1))
#for i in range(nbis1):
#    PhotonRate1[i] = Rate[0,E_min_index_max+i]
    #print('PhotonRate[i]:',PhotonRate[i])
#dNc_dtdE1 = PhotonRate1
#print('dNc_dtdE:',dNc_dtdE) #première valeur du tableau de Rate
##plt.plot(edges[:-1], dNc_dtdE, drawstyle = 'steps-post')
#plt.yscale('log')
#plt.xscale('log')
#plt.show()

PhotonCounts = np.zeros(shape=(nbis))
for i in range(nbis):
    PhotonCounts[i] = Rate[0, E_min_index_min+i]*Time_del[E_min_index_min+i]
dNc_dtdE2 = PhotonCounts/ dE

PhotonCounts0 = np.zeros(shape=(nbis0))
for i in range(nbis0):
    PhotonCounts0[i] = Rate[0, i]*Time_del[i]
dNc_dtdE20 = PhotonCounts0/ dE0

#PhotonCounts1 = np.zeros(shape=(nbis1))
#for i in range(nbis1):
#    PhotonCounts1[i] = Rate[0,E_min_index_max+i]*Time_del[E_min_index_max+i]
#dNc_dtdE21 = PhotonCounts1/ dE1
#plt.plot(edges[:-1], PhotonCounts, drawstyle = 'steps-post')
#plt.yscale('log')
#plt.xscale('log')
#plt.show()
#print('dNc_dtdE2:',dNc_dtdE2)


PhotonFlux = np.zeros(shape=(nbis))
for i in range(nbis):
    PhotonFlux[i] = Rate[0, E_min_index_min+i]/(Area*deltaE[0]-2)
dNc_dtdE3 = PhotonFlux/ dE

PhotonFlux0 = np.zeros(shape=(nbis0))
for i in range(nbis0):
    PhotonFlux0[i] = Rate[0,i]/(Area*deltaE[0]-2)
dNc_dtdE30 = PhotonFlux0/ dE0

#PhotonFlux1 = np.zeros(shape=(nbis1))
#for i in range(nbis1):
#    PhotonFlux1[i] = Rate[0, E_min_index_max+i]/(Area*deltaE[0]-2)
#dNc_dtdE31 = PhotonFlux1/ dE1

# Photon Flux total:
PhotonFlux_bis = np.zeros(shape=(n))
for i in range(n):
    PhotonFlux_bis[i] = Rate[0, i]/(Area*deltaE[0]-2)
dNc_dtdE3_bis = PhotonFlux_bis/ dE_bis
#dNc_dtdE3_bis_model=dNc_dtdE3_bis/(np.diag(matrix))
#print('dNc_dtdE3:',dNc_dtd3)
#plt.plot(edges[:-1], dNc_dtdE3, drawstyle = 'steps-post')
#plt.yscale('log')
#plt.xscale('log')
#plt.show()

def powerlaw1D(x,amplitude, x_0, alpha):
    power_law1D_true=amplitude*(x/x_0)**(-alpha)
    return power_law1D_true


def power_law(energy, norm, index):

    return norm * np.power(energy/100.,index)

# for ease
def differential_flux(e):

    return power_law(e, .01, -2)

def differential_flux_bis(e,amplitude, x_0, alpha):
    return powerlaw1D(e,amplitude, x_0, alpha)


# integral of the differential flux
def integral(e1, e2):
    return (e2 - e1) / 6.0 * (differential_flux(e1) + 4 * differential_flux((e1 + e2) / 2.0)+ differential_flux(e2))

def integralbis(e1, e2,amplitude, x_0, alpha):
    return (e2 - e1) / 6.0 * (differential_flux_bis(e1,amplitude, x_0, alpha) + 4 *
                              differential_flux_bis((e1 + e2)/ 2.0,amplitude, x_0, alpha) +
                              differential_flux_bis(e2,amplitude, x_0, alpha))

# true photon fluxes integrated over the photon bins of the response

# This is why we need to have Energy_low and Energy_high
true_fluxes = integral(Energy_low[:], Energy_hi[:])
#print('Energy_low:',len(Energy_low))
#print(' Energy_hi:', Energy_hi)
print('true_fluxes:',true_fluxes)
print('true_fluxes[E_min_index_min:E_min_index_max+1]:',true_fluxes[E_min_index_min:E_min_index_max+1])
print('true_fluxes[:E_min_index_min]:',true_fluxes[:E_min_index_min])
#print('true_fluxes[E_min_index_max+1:]:',true_fluxes[E_min_index_max+1:])
# dNp/(dt dA)
ed_bis = np.append(Energy_low,Energy_hi[-1])

if E_min_index_max==len(E_min):
    ed=ed_bis
else:
    ed=Energy_low[E_min_index_min:E_min_index_max+2]

ed0 = Energy_low[:E_min_index_min + 1]


#ed1=np.append(Energy_low[E_min_index_max+1:],Energy_hi[-1])
print('len(ed):',len(ed))
print('len(ed_bis):',len(ed_bis))
PhE = np.diff(ed)
PhE0 = np.diff(ed0)
#PhE1 = np.diff(ed1)
PhE_bis = np.diff(ed_bis)
print('PhE:',PhE)
fig, ax = plt.subplots()
#PL = true_fluxes/PhE
print('len(true_fluxes[E_min_index_min:E_min_index_max+1]):',len(true_fluxes[E_min_index_min:E_min_index_max+1]))
print('len(PhE):',len(PhE))
if E_min_index_max==len(E_min):
    PL = true_fluxes/PhE
else :
    PL = true_fluxes[E_min_index_min: E_min_index_max + 1] / PhE

PL0 = true_fluxes[:E_min_index_min]/PhE0
#PL1 = true_fluxes[E_min_index_max+1:]/PhE1

PL_bis=true_fluxes/PhE_bis
print('PL:',PL)
# print('PL size: ', PL.shape)
# P = PL/np.transpose(matrix)
# plt.plot(edges[:-1],P,  drawstyle = 'steps-post')
# plt.yscale('log')
# plt.xscale('log')
# plt.show()
#
# ax.step(ed[:-1],PL,where='post')
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlabel('Photon Energy')
# ax.set_ylabel(r'$\frac{d N_p}{dt dA}$')
# folded_photons = np.dot(PL, np.diag(matrix))
# plt.plot(ed[:-1], folded_photons,  drawstyle = 'steps-post')
# # #plt.plot(edges[:-1], folded_counts, drawstyle = 'steps-post', label='modeled data')
# plt.yscale('log')
# plt.xscale('log')
# plt.show()
print('true_fluxes.shape:',true_fluxes.shape)
print('matrix.shape:',matrix.shape)

folded_counts = np.dot(true_fluxes, np.transpose(matrix[:,E_min_index_min:E_min_index_max+1]).T)
folded_counts0 = np.dot(true_fluxes, np.transpose(matrix[:,:E_min_index_min]).T)
#folded_counts1 = np.dot(true_fluxes, np.transpose(matrix[:,E_min_index_max+1:]).T)


folded_counts_bis = np.dot(true_fluxes, np.transpose(matrix[:,:]).T)

print('folded_counts.shape:',folded_counts.shape)
#print('folded_counts:',folded_counts)
# PhFlux = np.zeros(shape=(n))
# for i in range(n):
#     PhFlux[i] = np.sum(true_fluxes[:])
# print('PhFlux: ', PhFlux.shape)
PhFlux = np.resize(PL, len(E_min[E_min_index_min:E_min_index_max+1]))
PhFlux0 = np.resize(PL0, len(E_min[:E_min_index_min]))
#PhFlux1 = np.resize(PL1, len(E_min[E_min_index_max+1:]))
PhFlux_bis = np.resize(PL_bis, len(E_min))
#print('PhFlux:',PhFlux)
#index_10 = np.where(edges_bis == E_min_index_max)[0][0]
Photons = dNc_dtdE3*(PhFlux/folded_counts)
#Photons_optimal_part1=Photons[:index_10+1]/np.diag(matrix)[:index_10+1]
#Photons_optimal_part2=Photons[index_10+1:]
#Photons_optimal=np.concatenate((Photons_optimal_part1,Photons_optimal_part2), axis=None)

Photons_bis = dNc_dtdE3_bis*(PhFlux_bis/folded_counts_bis)
#print('Photons:',Photons)


# dNc_dtdE3 = PhotonFlux/ dE == flux REEL de photons par bande d'énergie  (real data)

#First we will define a power law differential flux function, then its integral via Simpson’s rule, and finally calculate
# the photon fluxes of the photon model for each photon energy bound in the response matrix.
# These will be the “true” fluxes of the assumed photon model :
# PL = true_fluxes/PhE

#We still have a differential area in our flux and have not accounted for dispersion.
# Thus, we must convolve these true fluxes with our matrix.
# The matrix has units of area (cm2) accounting for the effective area of the detector.
# The operation is simply a dot product between the vector of photon fluxes we just computed and the matrix. Thus:
# folded_counts = np.dot(true_fluxes, np.transpose(matrix).T)

#Photons = dNc_dtdE3*(PhFlux/folded_counts) == produit entre le flux réel et on va dire
# l'atténuation des phénomènes physiques


## For the plotting,we need to get data<100 keV :


if E_min_index_max==len(E_min):
    index_100 = np.where(edges == 100)[0][0]
    #index_50 = np.where(edges == 46)[0][0]
    index_10 = np.where(edges == 10)[0][0]
    print('index_100:', index_100)
    #print('index_50:', index_50)
    plt.plot(edges[:index_100], Photons[:index_100],  drawstyle = 'steps-post')
    plt.plot(ed[:index_100], PL[:index_100], drawstyle = 'steps-post', label='modeled data')
    plt.xscale('log')
    plt.yscale('log')
    #plt.xscale('log')
    plt.legend()
    plt.show()
    # Now we can compare our assumed spectrum with the real data TOTAL :

    plt.plot(edges[:index_100], dNc_dtdE3[:index_100], drawstyle='steps-post', label='real data')
    plt.plot(edges[:index_100], folded_counts[:index_100], drawstyle='steps-post', label='modeled data')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.show()

    # Do Fit
    amplitude = 1  # par défaut
    x_0 = 2  # par défaut
    alpha = 2 # par défaut

    ## Essai 1 Comparaison :
    # Voici ce que l'on essaie de faire
    ## Do fit in the three part of energy band:
    popt, pcov = curve_fit(powerlaw1D, edges[:index_100], dNc_dtdE3[:index_100], [amplitude, x_0, alpha], method='lm')
    print('popt:', popt)
    # popt0, pcov0=curve_fit(powerlaw1D,edges0[:-1], dNc_dtdE30, [amplitude, x_0 , alpha])#, method='lm')
    # print('popt:',popt0)
    # popt1, pcov1=curve_fit(powerlaw1D,edges1[:-1], dNc_dtdE31, [amplitude, x_0 , alpha])#, method='lm')
    # print('popt:',popt0)
    #popt1_bis, pcov1_bis = curve_fit(powerlaw1D,edges[:index_100],dNc_dtdE3[:index_100], [amplitude, x_0, alpha])  # , method='lm')

    ## Affichage des résultats pour voir ce que cela donne :

    # Super cela fonctionne pour  les bandes : 10-16



    plt.plot(edges[:index_100], dNc_dtdE3[:index_100], drawstyle='steps-post', color='blue', label="Rate")
    plt.plot(edges[:index_100], powerlaw1D(edges[:index_100], popt[0], popt[1], popt[2]), drawstyle='steps-post', color='red',
             label="PowerLaw1D")
    #plt.plot(edges[:index_100], powerlaw1D(edges[:index_100], amplitude_opt ,x_0_opt ,  alpha_opt), drawstyle='steps-post',
             #color='green',
             #label="PowerLaw1D_opt")
    # plt.plot(edges[:-1],powerlaw1D(edges[:-1],popt[0],popt[1],popt[2])/dNc_dtdE3, drawstyle='steps-post', color='yellow', label="Rapport")
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    ##  Affichage des valeurs :

    print('dNc_dtdE3[0]:', dNc_dtdE3[0])
    print('dNc_dtdE3[index_10]:', dNc_dtdE3[index_10])
    print('dNc_dtdE3[index_100]:', dNc_dtdE3[index_100-1])

    print('powerlaw1D(edges[:index_100], popt[0], popt[1], popt[2])[0]',powerlaw1D(edges[:index_100], popt[0], popt[1], popt[2])[0])
    print('powerlaw1D(edges[:index_100], popt[0], popt[1], popt[2])[index_10]',powerlaw1D(edges[:index_100], popt[0], popt[1], popt[2])[index_10])
    print('powerlaw1D(edges[:index_100], popt[0], popt[1], popt[2])[index_100]',
          powerlaw1D(edges[:index_100], popt[0], popt[1], popt[2])[index_100-1])

    ## COUNT VS PHOTONS
    amplitude_photons =0.04  # par défaut
    x_0_photons = 1 # par défaut
    alpha_photons = 0.01# par défaut

    ## Essai 1 Comparaison :
    # Voici ce que l'on essaie de faire
    ## Do fit in the three part of energy band:
    matrix_inv = np.linalg.pinv(matrix)
    matrix_inv_bis = [abs(number) for number in matrix_inv]



    popt_photons, pcov_photons = curve_fit(powerlaw1D, edges[:index_100],Photons[:index_100],[amplitude_photons, x_0_photons, alpha_photons], method='lm')
    print('popt_photons:',popt_photons)
    #print('popt_photons:',popt_photons)
    plt.plot(edges[:index_100],Photons[:index_100], drawstyle='steps-post', color='blue', label="Photons_real")
    #plt.plot(edges[:index_100],powerlaw1D(edges[:index_100], popt_photons[0], popt_photons[1], popt_photons[2]), drawstyle='steps-post',color='red',label="PowerLaw1D_Photons")
    #plt.plot(edges_bis[:index_100],(PhFlux[:index_100]/folded_counts[:index_100])*powerlaw1D(edges[:index_100], popt[0], popt[1], popt[2]), drawstyle='steps-post', color='cyan',
    #         label="Photon_real_2")
    plt.plot(edges_bis[:index_100],powerlaw1D(edges[:index_100],  popt_photons[0],  popt_photons[1], popt_photons[2]), drawstyle='steps-post',color='green',label="Photon_real_2_fit")


    #plt.plot(edges[:index_100], powerlaw1D(edges[:index_100], popt[0], popt[1], popt[2]), drawstyle='steps-post',
     #        color='orange',
      #       label="PowerLaw1D")
    #plt.plot(edges[:index_100], 0.1*powerlaw1D(edges[:index_100], popt_photons[0], popt_photons[1], popt_photons[2]),drawstyle='steps-post', color='green', label="PowerLaw1D_Photons_2")
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.show()

    print('Photons[0]:', Photons[0])
    print('Photons[index_10]:', Photons[index_10])
    print('Photons[index_100]:',Photons[index_100-1])

    print('powerlaw1D(edges[:index_100],  popt_photons[0],  popt_photons[1], popt_photons[2])[0]',
          powerlaw1D(edges[:index_100],  popt_photons[0],  popt_photons[1], popt_photons[2])[0])
    print('powerlaw1D(edges[:index_100],  popt_photons[0],  popt_photons[1], popt_photons[2])[index_10]',
          powerlaw1D(edges[:index_100],  popt_photons[0],  popt_photons[1], popt_photons[2])[index_10])
    print('powerlaw1D(edges[:index_100],  popt_photons[0],  popt_photons[1], popt_photons[2])[index_100]',
          powerlaw1D(edges[:index_100],  popt_photons[0],  popt_photons[1], popt_photons[2])[index_100-1])

else :
    index_100 = np.where(edges_bis == 100)[0][0]
    index_10 = E_min_index_max #np.where(edges == E_min_index_max)[0][0]
    plt.plot(edges[:-1], Photons, drawstyle='steps-post')
    plt.plot(ed[:-1], PL, drawstyle='steps-post', label='modeled data')
    plt.xscale('log')
    plt.yscale('log')
    # plt.xscale('log')
    plt.legend()
    plt.show()
    plt.plot(edges[:-1], dNc_dtdE3, drawstyle='steps-post', label='real data')
    plt.plot(edges[:-1], folded_counts, drawstyle='steps-post', label='modeled data')
    plt.xscale('log')
    plt.yscale('log')
    # plt.xscale('log')
    plt.legend()
    plt.show()

    # Do fit
    amplitude = 1  # par défaut
    x_0 = 2  # par défaut
    alpha = 1  # par défaut
    # Voici ce que l'on essaie de faire
    ## Do fit in the three part of energy band:
    popt, pcov = curve_fit(powerlaw1D, edges[:-1], dNc_dtdE3, [amplitude, x_0, alpha], method='lm')
    print('popt:', popt)
    popt0, pcov0=curve_fit(powerlaw1D,edges0[:-1], dNc_dtdE30, [amplitude, x_0 , alpha])#, method='lm')
    # print('popt:',popt0)
    # popt1, pcov1=curve_fit(powerlaw1D,edges1[:-1], dNc_dtdE31, [amplitude, x_0 , alpha])#, method='lm')
    # print('popt:',popt0)
    #popt1_bis, pcov1_bis = curve_fit(powerlaw1D, np.concatenate((edges, edges1[1:]))[:-1],
    #                                 np.concatenate((dNc_dtdE3, dNc_dtdE31)), [amplitude, x_0, alpha])  # , method='lm')

    #Affichage on the energy band
    plt.plot(edges[:-1], dNc_dtdE3, drawstyle='steps-post', color='blue', label="Rate")
    plt.plot(edges[:-1], powerlaw1D(edges[:-1], popt[0], popt[1], popt[2]), drawstyle='steps-post', color='red',
             label="PowerLaw1D")
    # plt.plot(edges[:-1],powerlaw1D(edges[:-1],popt[0],popt[1],popt[2])/dNc_dtdE3, drawstyle='steps-post', color='yellow', label="Rapport")
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    # Affichage on the energy band
    plt.plot(edges0[:-1], dNc_dtdE30, drawstyle='steps-post', color='blue', label="Rate")
    plt.plot(edges0[:-1], powerlaw1D(edges0[:-1], popt0[0], popt0[1], popt0[2]), drawstyle='steps-post', color='red',
             label="PowerLaw1D")
    # plt.plot(edges[:-1],powerlaw1D(edges[:-1],popt[0],popt[1],popt[2])/dNc_dtdE3, drawstyle='steps-post', color='yellow', label="Rapport")
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    # Affichage on the energy < 100 keV

    ## global Results :
    dNc_dtdE3_optimal=np.zeros_like(dNc_dtdE3_bis)
    for i in range(nbis0):
        dNc_dtdE3_optimal[i]=powerlaw1D(edges0[:-1],popt0[0],popt0[1],popt0[2])[i]
    for i in range(n-nbis0):
        dNc_dtdE3_optimal[i+nbis0] = powerlaw1D(edges_bis[nbis0:], popt[0], popt[1], popt[2])[i]
    plt.plot(edges_bis[:index_100], dNc_dtdE3_bis[:index_100], drawstyle='steps-post', color='blue', label="Rate")
    plt.plot(edges_bis[:index_100], powerlaw1D(edges_bis[:index_100], popt[0], popt[1], popt[2]), drawstyle='steps-post', color='red',
             label="PowerLaw1D")
    plt.plot(edges_bis[:index_100], dNc_dtdE3_optimal[:index_100],
             drawstyle='steps-post', color='green',
             label="PowerLaw1D_opt")
    # plt.plot(edges[:-1],powerlaw1D(edges[:-1],popt[0],popt[1],popt[2])/dNc_dtdE3, drawstyle='steps-post', color='yellow', label="Rapport")
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')

    plt.show()

    ## Affichage

    print('dNc_dtdE3_bis[0]:', dNc_dtdE3_bis[0])
    print('dNc_dtdE3_bis[E_min_index_min]:', dNc_dtdE3_bis[E_min_index_min])
    print('dNc_dtdE3_bis[index_100]:', dNc_dtdE3_bis[index_100-1])

    print('dNc_dtdE3_optimal[0]',dNc_dtdE3_optimal[0])
    print('dNc_dtdE3_optimal[E_min_index_min]',dNc_dtdE3_optimal[E_min_index_min])
    print('dNc_dtdE3_optimal[index_100]',dNc_dtdE3_optimal[index_100-1])

    ## COUNT VS PHOTONS
    amplitude_photons = 1  # par défaut
    x_0_photons = 1  # par défaut
    alpha_photons = 2  # par défaut

    ## Essai 1 Comparaison :
    # Voici ce que l'on essaie de faire
    ## Do fit in the three part of energy band:

    popt_photons, pcov_photons = curve_fit(powerlaw1D, edges[:-1], Photons, [amplitude_photons, x_0_photons, alpha_photons], method='lm')
    popt0_photons, pcov0_photons = curve_fit(powerlaw1D, edges0[:-1], dNc_dtdE30, [amplitude, x_0, alpha])  # , method='lm')
    # Affichage on the energy band
    print('popt_photons:',popt_photons)

    ## global Results :
    Photons_optimal = np.zeros_like(Photons_bis)
    for i in range(nbis0):
        Photons_optimal[i] = powerlaw1D(edges0[:-1], popt0_photons[0], popt0_photons[1], popt0_photons[2])[i]
    for i in range(n - nbis0):
        Photons_optimal[i + nbis0] = powerlaw1D(edges_bis[nbis0:], popt_photons[0], popt_photons[1], popt_photons[2])[i]

    plt.plot(edges_bis[:index_100], Photons_bis[:index_100], drawstyle='steps-post', color='blue', label="Photon_real")
    plt.plot(edges_bis[:index_100], powerlaw1D(edges_bis[:index_100], popt_photons[0], popt_photons[1], popt_photons[2]), drawstyle='steps-post', color='red', label="Photon_real_fit")
    plt.plot(edges_bis[:index_100], Photons_optimal[:index_100], drawstyle='steps-post', color='green', label="Photon_real_fit_opt")
    #plt.plot(edges_bis[:index_100], powerlaw1D(edges_bis[:index_100], popt_photons[0], popt_photons[1], popt_photons[2]), drawstyle='steps-post', color='red',
     #        label="PowerLaw1D")
    # plt.plot(edges[:-1],powerlaw1D(edges[:-1],popt[0],popt[1],popt[2])/dNc_dtdE3, drawstyle='steps-post', color='yellow', label="Rapport")
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

    ## Affichage

    print('Photons_bis[0]:', Photons_bis[0])
    print('Photons_bis[E_min_index_min]:',Photons_bis[E_min_index_min])
    print('Photons_bis[index_100]:', Photons_bis[index_100-1])

    print('Photons_optimal[0]', Photons_optimal[0])
    print('Photons_optimal[E_min_index_min]', Photons_optimal[E_min_index_min])
    print('Photons_optimal[index_100]', Photons_optimal[index_100-1])




















#d = np.dot(dNc_dtdE3*true_fluxes)
# plt.plot(edges[:-1],folded_counts/dE, drawstyle = 'steps-post', label='modeled data')
# plt.plot(edges[:-1],PL[:77], drawstyle = 'steps-post', label='real data')
# plt.yscale('log')
# plt.xscale('log')
# plt.show()

# Result : We can see how the effective area and effects of energy dispersion have altered the shape of our power law.
# This is why it can be difficult to spectrally model x-ray spectra measured by instruments
# that exhibit strong energy dispersion

# Do_fit pour voir ce que cela fait normalement mais ici cela ne fonctionne pas:
#fitter = fitting.LevMarLSQFitter()
#Essai=np.resize(true_fluxes[E_min_index_min:E_min_index_max+1], len(E_min[E_min_index_min:E_min_index_max+1]))
#result = fitter(models.PowerLaw1D(), dNc_dtdE3, true_fluxes[E_min_index_min:E_min_index_max+1])
#print('dNc_dtdE3:',dNc_dtdE3.shape)
#print('true_fluxes:',Essai.shape)




#print('amplitude_opt = %.3g' %popt[0])
#print('x_0_opt = %.3g' %popt[1])
#print('alpha_opt = %.3g' %popt[2])

#,weights=1.0 /true_fluxes) amplitude=1, x_0=3, alpha=50, fixed = {'alpha': True}
#print('popt0:',popt0)
#print('pcov0:',pcov0)
#print('popt1:',popt1)
#print('pcov1:',pcov1)

## What can we do if we have issue of convergence of the power law
# Error : RuntimeWarning: invalid value encountered in power
#   return amplitude * xx ** (-alpha)

## Affichage des résultats pour voir ce que cela donne :

# Super cela fonctionne pour  les bandes : 10-16

#plt.plot(edges[:-1],dNc_dtdE3, drawstyle='steps-post',color='blue', label="Rate")
#plt.plot(edges[:-1],powerlaw1D(edges[:-1],popt[0],popt[1],popt[2]), drawstyle='steps-post', color='red', label="PowerLaw1D")
#plt.plot(edges[:-1],powerlaw1D(edges[:-1],popt[0],popt[1],popt[2])/dNc_dtdE3, drawstyle='steps-post', color='yellow', label="Rapport")
#plt.legend()
#plt.yscale('log')
#plt.show()

#plt.xscale('log')
#plt.ylim(ymax = 100, ymin = 0.1)
#plt.plot(edges0[:-1],dNc_dtdE30, drawstyle='steps-post',color='blue', label="Rate")
#plt.plot(edges0[:-1],powerlaw1D(edges0[:-1],popt0[0],popt0[1],popt0[2]), drawstyle='steps-post', color='red', label="PowerLaw1D")
#plt.plot(edges[:-1],powerlaw1D(edges[:-1],popt[0],popt[1],popt[2])/dNc_dtdE3, drawstyle='steps-post', color='yellow', label="Rapport")
#plt.legend()
#plt.yscale('log')
#plt.show()

#plt.plot(edges1[:-1],dNc_dtdE31, drawstyle='steps-post',color='blue', label="Rate")
#plt.plot(edges1[:-1],powerlaw1D(edges1[:-1],popt1[0],popt1[1],popt1[2]), drawstyle='steps-post', color='red', label="PowerLaw1D")
#plt.plot(edges[:-1],powerlaw1D(edges[:-1],popt[0],popt[1],popt[2])/dNc_dtdE3, drawstyle='steps-post', color='yellow', label="Rapport")
#plt.legend()
#plt.yscale('log')
#plt.show()
## We can see here that in the theard energy band we don't have the convergence of the curvefit


## global Results :
#dNc_dtdE3_optimal=np.zeros_like(dNc_dtdE3_bis)
#for i in range(nbis0):
#    dNc_dtdE3_optimal[i]=powerlaw1D(edges0[:-1],popt0[0],popt0[1],popt0[2])[i]
#for i in range(nbis):
#        dNc_dtdE3_optimal[nbis0+i] = powerlaw1D(edges[:-1], popt[0], popt[1], popt[2])[i]
#for i in range(nbis1):
#        dNc_dtdE3_optimal[nbis0+nbis+i] = powerlaw1D(edges1[:-1], popt1[0], popt1[1], popt1[2])[i]









