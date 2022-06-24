## ESSENCE: functions for evaluating spatially correlated noise in the interferometric images.

**ESSENCE** is a python package for evaluating the statistical significance of the image analysis and signal detection under correlated noise in the interferometric images (e.g., ALMA, NOEMA), namely, Evaluating Statistical Significance under Noise Correlation (ESSENCE).  
This code does the following things for you:
1. characterizing the noise statistical properties of interferometric images by noise autocorrelation function (ACF).
2. from measured noise ACF, compute the noise in the spatially integrated quantities (e.g., flux, spectrum) with a given aperture. 
3. from measured noise ACF, simulate noise maps with the same correlation property.

### Requirements (Tested version):

>  python-------------(3.7.7)  
>  astropy------------(4.3.1)  
>  spectral_cube------(0.6.0)  
>  numpy--------------(1.21.5)  
>  scipy--------------(1.7.3)  
>  multiprocess-------(0.70.13)  
>  functools 

### Installation:


Not required. Just git pull the software to a desired directory.    
>        $ git clone 
>        $ cd essence
  
Usage:


see [tutorial]() 

Contacts:


Takafumi Tsukui: tsukuitk23@gmail.com
