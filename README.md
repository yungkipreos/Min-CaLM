# Min-CaLM
This python package can perform automated mineral compositional analysis on debris disk spectra. It will determine the minerals that are present and what their relative abundances are within the debris disk. This code is used in a paper that has been submitted to the American Journal of Undergraduate Research (AJUR), titled "An Unbiased Mineral Compositional Analysis Technique for Debris Disks" by Yung Kipreos and Dr. Inseok Song.

# Getting Started
To use Min-CaLM's code, the following python libraries must be imported: numpy, astropy.io, pylab, scipy.optimize, matplotlib, scipy, decimal.

## Debris Disk Spectrum Requirements
Min-CaLM assumes several things about the target debris disk spectrum. It assumes that the spectrum's data file is in atwo-column format where the first column is the wavelength data and the second column is the flux data, and also that both of the columns are in increasing order (according to the wavelengths data column).

To perform mineral compositional analysis on a debris disk spectrum, there must be prominent silicate mineral features present. Below is an example of a debris disk spectrum with no silicate mineral features:
![](Min-CaLM/HD192758_Debris_Disk_Spectrum.png)
