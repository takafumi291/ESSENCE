ESSENCE: functions for evaluating correlated noise in the interferometric images.
**ESSENCE** is a python package for evaluating the statistical significance of the image analysis and signal detection under correlated noise in the interferometric images (e.g., ALMA, NOEMA), namely, Evaluating Statistical Significance under Noise Correlation (ESSENCE).  
This code does the following things for you:
1. characterizing the noise statistical properties of interferometric images by noise autocorrelation function (ACF).
2. from measured noise ACF, compute the noise in the spatially integrated quantities (e.g., flux, spectrum) with a given aperture. 
3. from measured noise ACF, simulate noise maps with the same correlation property.
