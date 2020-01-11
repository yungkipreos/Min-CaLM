import numpy as np
from astropy.io import ascii
from pylab import *
from scipy.optimize import nnls
from matplotlib import pyplot as plt
from scipy import interpolate
from decimal import *
import pandas as pd


#--------------IMPORTING MINERAL DATA---------------#

# Routine to read mineral spectra into a Pandas DataFrame
data_dir = "Minerals"
min_data=[("Diopside", "Diopside_555.dat"),
	("Dolomite", "Dolomite_549.dat"),
	("Enstatite", "Enstatite_550.dat"),
	("Fayalite", "Fayalite_557.dat"),
	("Forsterite", "Forsterite_465.dat"),
	("Al2O3", "Al2O3.dat"),
	("Albite","Albite_438.dat"),
	("Anorthite","Anorthite_564.dat"),
	("Augite","Augite_480.dat"),
	("Bronzite","Bronzite_549.dat"),
	("Coesite","Coesite_1970.dat"),
	("Corundum","Corundum.dat"),
	("GammaAlumina", "GammaAlumina.dat"),
	("Halloysite","Halloysite.dat"),
	("Hedenbergite","Hedenbergite_600.dat"),
	("Hematite","Hematite_474.dat"),
	("Hornblende","Hornblende_488.dat"),
	("Hypersthene","Hypersthene.dat"),
	("Magnesite","Magnesite_531.dat"),
	("Magnetite","Magnetite_487.dat"),
	("Montmorillonite","montmorillonite_584.dat"),
	("Nontronite","Nontronit_587.dat"),
	("Olivine","Olivine_661.dat"),
	("Pigeonite","Pigeonite_3668.dat"),
	("Pyroxene","pyroxene.dat"),
	("Quartz","Quartz_1969.dat"),
	("Saponite","Saponite_579.dat"),
	("Siderite","Siderite_530.dat"),
	("SiO","SiOGas.dat"),
	("Talc","Talc.dat"),
	("Ti3O5","Ti3O5.dat"),
	("Tridymite","Tridymite_1967.dat"),
	("Wollastonite","Wollastonite.dat")]

# Create two arrays (names, datafiles) from the above record.
names = [x[0] for x in min_data]
datafiles = [x[1] for x in min_data]

# Now go through the list of minerals and read in their spectra
min_spectra = [] #An empty list to read the spectra into
for idx,dfile in enumerate(datafiles):
	sp = ascii.read(dfile, names=['wave','flux'])
	min_spectra.append(sp)


"""
DEFINITIONS

names = This is a list containing the names of the minerals as strings
datafiles = This is a list containing the names of the mineral FILES ie ("Dolomite_555.dat")

min_spectra = This is a two-dimensional list containing all of the mineral data. To get a specific mineral's spectra
	      You have to reference the mineral's index in the list, and then whether you want the flux or the wavelength data
	      EXAMPLE:
	    	  idx = names.index('Enstatite') 		#Find the index of the mineral by using the "names" array
                  enstatite_wave = min_spectra[idx]['wave']	#Specify if you want the flux or the wavelength
	     	  enstatite_flux = min_spectra[idx]['flux']
"""	






#--------------ASKING FOR USER INPUT---------------#
"""
The target spectrum is assumed to be in a two-column format. 
The first column is the wavelengths and the second column is 
the fluxes. Both columns must be in increasing order according
to the wavelengths column.
"""
InputSpectra = input("Please input the target's spectra (target spectra must be kept in the same folder as Min-CaLM.py): \n")
with open(InputSpectra):
	TargetSpectra_wave, TargetSpectra_flux = np.loadtxt(InputSpectra, unpack = True)






#---------------INTERPOLATING NECESSARY SPECTRA---------------#

interpolated_min_spectra = min_spectra.copy()
interp_spectra = min_spectra[19]['wave']
interp_spectra = interp_spectra.tolist()
#interp_spectra = np.flip(interp_spectra)

TargetSpectra_flux = np.interp(interp_spectra,TargetSpectra_wave, TargetSpectra_flux)
for idx,mineral in enumerate(names):
	interpolated_min_spectra_flux = np.interp(interp_spectra,min_spectra[idx]['wave'], min_spectra[idx]['flux'])
	interpolated_min_spectra[idx] = interpolated_min_spectra_flux
min_fluxes = interpolated_min_spectra







#----------------DUST BLACKBODY REMOVAL-----------------#

def planck(w,temp): 
	w = w/1.E4 
	c1 = 3.7417749E-5 
	c2 = 1.4387687 
	val = c2/w/temp 
	flux1 = c1/(w**5 * np.expm1(val) ) 
	BBflux = flux1*(w**2.0)*(5.353E-11)/np.pi 
	return BBflux 

# This calculates the Planck spectrum at a specified temperature 
hornblende_wave = min_spectra[19]['wave']
hornblende_wave = np.flip(hornblende_wave)
Cont=planck(hornblende_wave, 502)  	#ENTER DUST BB TEMP HERE
Cont=Cont/median(Cont)

# The BB Continuum is removed by division
i = 0
while (i < len(TargetSpectra_flux)):
	TargetSpectra_flux[i] = TargetSpectra_flux[i] / Cont[i]
	i = i + 1







#-------------NON-NEGATIVE LEAST SQUARE MINIMIZATION--------------#
"""
NON-NEGATIVE LEAST SQUARE MINIMIZATION COMMENTS SECTION
- In this section, the recreated target spectrum is calculated.
- The target spectrum and mineral spectra are put into a linear 
  system of equations that is then solved by using the 
  spicy.optimize.nnls() function. This function calculates the 
  relative abundance of each mineral predicted to be present in 
  the target spectrum.
- xMatrix = An array with all of the mineral data
- targetMatrix = The (reshaped) target flux array 
"""

TargetMatrix = TargetSpectra_flux.copy()
TargetMatrix.reshape((len(TargetSpectra_flux),1))

xMatrix = interpolated_min_spectra[0]
i = 1

while (i < len(min_fluxes)):
	xMatrix = np.column_stack((xMatrix,min_fluxes[i]))
	i = i + 1

wMatrix= nnls(xMatrix, TargetSpectra_flux) 
wMatrix = wMatrix[0:1]
wMatrix = wMatrix[0]







#--------MULTIPLYING THE RELATIVE ABUNDANCES WITH THE MINERAL SPECTRA----------#


i = 0
sum = []
recombinedSpectrum = np.empty([len(TargetSpectra_flux)])

while (i < len(min_fluxes)):
	sum = min_fluxes[i] * wMatrix[i]
	recombinedSpectrum = recombinedSpectrum + sum
	i = i + 1
wavelength = min_spectra[0]['wave']
wavelength = interp_spectra







#--------MULTIPLYING THE OVERALL SPECTRA BY THE BB CONTINUUM--------#

#The blackbody is recombined with each spectrum for display purposes
i = 0 
while (i < len(recombinedSpectrum)):
	recombinedSpectrum[i] = recombinedSpectrum[i] *Cont[i]
	TargetSpectra_flux[i] = TargetSpectra_flux[i] *Cont[i]
	i = i + 1

wavelength = wavelength[20:]
recombinedSpectrum = recombinedSpectrum[20:]
TargetSpectra_flux = TargetSpectra_flux[20:]

#Displaying the final spectrum:
plt.title ("BD+20 307")
plt.ylabel('Flux (Jy)')
plt.xlabel('Wavelength (\u03bcm)')
plt.plot(wavelength,TargetSpectra_flux, 'r')
plt.plot(wavelength, recombinedSpectrum, 'b')
plt.show()

i = 0
totalSum = 0
while (i< len(wMatrix)):
	totalSum = totalSum + wMatrix[i]
	i = i + 1







#-----------------DISPLAYING THE RESULTS------------------#
print("\n") #aesthetic
results = pd.DataFrame()
i = 0
while ( i < len(wMatrix)):
	if (wMatrix[i] > 0):
		percentage = round(Decimal((wMatrix[i]/totalSum)*100), 2)
		percentage = str(percentage) + "%"
		df = pd.DataFrame([[names[i], "-",percentage]], columns=['MINERAL',' ','ABUNDANCE'])
		results = results.append(df)
	i = i + 1

print(results.to_string(index=False))



