## ESSENCE: functions for evaluating spatially correlated noise in the interferometric images.

**ESSENCE** is a Python package for evaluating the statistical significance of image analysis and signal detection under correlated noise in interferometric images (e.g., ALMA, NOEMA), namely, Evaluating Statistical Significance undEr Noise CorrElation.  
This code does the following things for you:
1. measuring noise autocorrelation function (ACF) which fully characterizes the statistical properties of spatially correlated noise in the interferometric image.  
2. computing the noise in the spatially integrated quantities (e.g., flux, spectrum) with a given aperture. 
3. simulating noise maps with the same correlation property.
4. constructing a covariance matrix from noise ACF, which can be used for a 2D image or 3D cube model fitting.

Detailed formulation of ESSENCE and its application are presented in [Tsukui et al. 2023](https://www.spiedigitallibrary.org/journals/Journal-of-Astronomical-Telescopes-Instruments-and-Systems/volume-9/issue-01/018001/Estimating-the-statistical-uncertainty-due-to-spatially-correlated-noise-in/10.1117/1.JATIS.9.1.018001.full?SSO=1).

### Requirements:
	
| Packages | Tested version |
| --------------:|---------------:|
| python | 3.7.7 |
| astropy | 4.3.1 |
| spectral_cube | 0.6.0 |
| numpy | 1.21.5 |
| scipy | 1.7.3 |
| multiprocess | 0.70.13 |
| functools | |

### Installation:

Not required. Git clone the software to a desired directory.    
>        $ git clone https://github.com/takafumi291/ESSENCE.git
>        $ cd essence

### Example data:

For running tutorial.ipynb, please download [example data](https://drive.google.com/file/d/1h0wEPHVebVSjl803r9LnQyBTxfoU2kBY/view?usp=sharing), unzip, and place it in the same directory of the ipynb file.  
The data is from Tsukui and Iguchi 2021, Sci (ADS/JAO.ALMA2017.1.00394.S PI=Gonzalez Lopez, Jorg)

### Usage:

See [tutorial](https://github.com/takafumi291/ESSENCE/blob/main/Tutorial.ipynb) for a quick example. 

### Contacts:
I am open to collaborations, e.g., any suggestion, feedback, or directly improving my codes. I am also happy to help with any difficulties you encounter using my codes. Feel free to contact me! 

Takafumi Tsukui: tsukuitk23_at_gmail.com
