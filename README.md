# Min-CaLM
This python package can perform automated mineral compositional analysis on debris disk spectra. It will determine the minerals that are present and what their relative abundances are within the debris disk. This code is used in a paper that has been submitted to the American Journal of Undergraduate Research (AJUR), titled "An Unbiased Mineral Compositional Analysis Technique for Debris Disks" by Yung Kipreos and Dr. Inseok Song.

# Getting Started
To use Min-CaLM's code, the following python libraries must be imported: numpy, astropy.io, pylab, scipy.optimize, matplotlib, scipy, decimal.

## Debris Disk Spectrum Requirements
Min-CaLM has several assumptions regarding the format of the target debris disk spectrum. It assumes that the spectrum's data file is in a two-column format where the first column is the wavelength data and the second column is the flux data, and also that both of the columns are in increasing order (according to the wavelengths data column). The wavelength range of the debris disk spectrum should be around from ~5 to ~45 microns.

To perform mineral compositional analysis on a debris disk spectrum, there must be prominent silicate mineral features present. Below is an example of a debris disk spectrum with no silicate mineral features:

<img src="/HD192758_Debris_Disk_Spectrum.png" width = 300 >

And below is an example of a debris disk spectrum that displays prominent silicate mineral features:

<img src="/HD15407_Debris_Disk_Spectrum.png" width = 300 >

## Running the Min-CaLM Program
It is reccommended to run the Min-CaLM program through the command line. First download both the Min-CaLM.py file and the mineral spectra files into the same directory. There are some example debris disk spectra that can be downloaded as well. The Min-CaLM.py file, mineral spectra files, and example debris disk spectra files can be all found in this Git-hub repository. 

To run Min-CaLM, go to the command line and navigate to the directory that contains all of the above files. Then type "ipython" into the command line. Now that you are in ipython, type "run Min-CaLM.py" to run the program. The program will ask you to "Please input the target's spectrum (target spectrum must be kept in the same folder as Min-CaLM.py):". Type in the file name of the debris disk spectrum into the command line and press enter. The program will then display a figure of the recreated spectrum (blue) plotted over the original debris disk spectrum (red). In the command line, a table will be displayed that shows the minerals that were determined to be present in the disk and their relative abundances.

# Tutorial 1
This tutorial is a more indepth description of how to use the Min-CaLM program to perform mineral compositional analysis on debris disk spectra. This tutorial assumes that the star's photosphere blackbody contributions has aleady been removed from the spectrum. The resulting spectrum is the debris disk spectrum that contains the mineral spectra and the blackbody spectrum produced by the heated dust/debris. Here, I will use the debris disk around BD+20 307 as an example. BD+20 307 is a dusty binary star system in the Aries constellation. 

The BD+20 307 debris disk data file is in a two column format, a section of which is shown below. The left column is the wavelength data and the right column is the flux data. Notice that the data is in ascending order according to the wavelength data. 

<img src="/BD+20_307_Sample_Data.png" width = 300 >

The spectrum of this debris disk is:

<img src="/HD_69830_Spectrum.png" width = 300 >

Before the Min-CaLM program can be used on the debris disk data, the dust's blackbody spectrum must be removed















