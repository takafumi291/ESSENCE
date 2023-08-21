"""
    Copyright (C) 2022, Takafumi Tsukui
    E-mail: tsukuitk23@gmail.com

    Updated versions of the software are available from my web page
    https://sites.google.com/view/takafumitk/home?authuser=0

    If you have found this software useful for your research,
    I would appreciate an acknowledgement to the use of the
    "Essense (Evaluating Statistical Significance undEr CorrElation) of Tsukui (2022)"

    This software is provided as is without any warranty whatsoever.
    Permission to use, for non-commercial purposes is granted.
    Permission to modify for personal or internal use is granted,
    provided this copyright and disclaimer are included unchanged
    at the beginning of the file. All other rights are reserved.
    In particular, redistribution of the code is not allowed.

"""


import os
import sys
import matplotlib.pyplot as plt
import astropy.units as u
from spectral_cube import SpectralCube
import numpy as np
from astropy.io import fits
import time
#from multiprocessing import Pool
#https://stackoverflow.com/questions/41385708/multiprocessing-example-giving-attributeerror
import multiprocess as mp
import functools
from scipy import signal
from astropy.stats import sigma_clip
from scipy import ndimage
import scipy
from astropy.modeling import models, fitting

#################################################################################################

def mk_noisemap(image, pb, pbfrac=0.7, sigma=4):
    """
    make noisemap from data by (1) sigma clipping emission regions 
                               (2) growing the clipped region with the beam FWHM (pixels)
                               (3) clip region where primaly beam < pbfrac
    The clipped regions are filled with nan values. 
    
    image: data (2d image or cube), spectral-cube class, if image.shape==3 we assume (v, x, y)
    pb: primary beam image, spectral-cube class 
    pbfrac:  minimum primary beam fraction not to clip 
    sigma: threthold for sigma clipping
    
    return masked cube or image with nan
    """
    # first mask primary beam>pbfrac
    Pbmask=mk_pbmask(image._data,pb,pbfrac=pbfrac,plot=False)

    if len(image.shape)==3:
        Pbmasked=np.copy(image._data)
        Pbmasked[~Pbmask]=np.nan
        Sigmacliped = sigma_clip(Pbmasked, sigma=sigma, maxiters=5,masked=True, axis=(1,2))
        growpix=image.hdu.header['BMAJ']/image.hdu.header['CDELT2']
        Yy, Xx = np.indices([int(3*growpix),int(3*growpix)], dtype='float')
        CX=(int(3*growpix)-1)/2
        CY=(int(3*growpix)-1)/2
        Kernel = (((Yy-CY)**2 + (Xx-CX)**2)**0.5<growpix).astype('float')
        Finalmask=np.array([signal.fftconvolve(np.copy(i.mask.astype(float)), Kernel, mode="same")>1 for i in Sigmacliped])
        Finalmasked=np.copy(image.hdu.data)
        Finalmasked[Finalmask]=np.nan

        return Finalmasked

    if len(image.shape)==2:
        Pbmasked=np.copy(image.hdu.data)
        Pbmasked[~Pbmask]=np.nan
        Sigmacliped = sigma_clip(Pbmasked, sigma=sigma, maxiters=5,masked=True, axis=None)
        growpix=image.hdu.header['BMAJ']/image.hdu.header['CDELT2']

        Yy, Xx = np.indices([int(3*growpix),int(3*growpix)], dtype='float')
        CX=(int(3*growpix)-1)/2
        CY=(int(3*growpix)-1)/2
        Kernel = (((Yy-CY)**2 + (Xx-CX)**2)**0.5<growpix).astype('float')
        Finalmask=signal.fftconvolve(np.copy(Sigmacliped.mask.astype(float)), Kernel, mode="same")>1
        Finalmasked=np.copy(image.hdu.data)
        Finalmasked[Finalmask]=np.nan

        return Finalmasked
        
#######################################################################################################################

def mk_pbmask(image, pb, pbfrac=0.7, plot=True):
    """mask image with nan, where pb<pbfrac
    """
    if len(image.shape)==3:
        YY, XX = np.indices(image.shape[1:], dtype='float')
        Cx=image.shape[1]/2
        Cy=image.shape[2]/2
    if len(image.shape)==2:
        YY, XX = np.indices(image.shape[:], dtype='float')
        Cx=image.shape[0]/2
        Cy=image.shape[1]/2

    Radius = ((YY-Cy)**2 + (XX-Cx)**2)**0.5
    if len(image.shape)==3:
        Mask = (pb._data>pbfrac)
    if len(image.shape)==2:
        Mask = (pb._data[0]>pbfrac)
    return Mask
    
############################################################################################################################################## 
    
def zoompeak(image,Width, plot=False):
    cenx,ceny=np.where(image==np.nanmax(image))
    if plot:
        plt.imshow(image[cenx[0]-Width:cenx[0]+Width+1,ceny[0]-Width:ceny[0]+Width+1],extent=[Width-0.5,-Width+.5,-Width-0.5,Width+.5])
        plt.colorbar()
        plt.xlabel(r"$\Delta$x (pixels)")
        plt.ylabel(r"$\Delta$y (pixels)")
    return image[cenx[0]-Width:cenx[0]+Width+1,ceny[0]-Width:ceny[0]+Width+1]
    
###############################################################################################################################################    

def zoomcen(image,Width, plot=False):
    """
    zoom Width x Width pixels of image
    """
    image_s=image[int((image.shape[0]-1)/2)-Width:int((image.shape[0]-1)/2)+Width+1,int((image.shape[0]-1)/2)-Width:int((image.shape[0]-1)/2)+Width+1]
    try:
        s=np.where(image_s==np.max(image_s))==(np.array([int((image_s.shape[0]-1)/2)]),np.array([int((image_s.shape[0]-1)/2)]))
        if s:
            print("peak coincide with the image center")
        else:
            print("peak do not coicide with the image center")
    except:
        print("there seem to be multiple peaks")
    if plot:
        plt.imshow(image[int((image.shape[0]-1)/2)-Width:int((image.shape[0]-1)/2)+Width+1,int((image.shape[0]-1)/2)-Width:int((image.shape[0]-1)/2)+Width+1],extent=[Width-0.5,-Width+.5,-Width-0.5,Width+.5])
        plt.colorbar()
        plt.xlabel(r"$\Delta$x (pixels)")
        plt.ylabel(r"$\Delta$y (pixels)")
    return image[int((image.shape[0]-1)/2)-Width:int((image.shape[0]-1)/2)+Width+1,int((image.shape[0]-1)/2)-Width:int((image.shape[0]-1)/2)+Width+1]

####################################################################################################################################################

def shift_2d_replace(data, dx, dy, constant=np.nan):
    """
    Shifts the array in two dimensions while setting rolled values to constant
    :param data: The 2d numpy array to be shifted
    :param dx: The shift in x
    :param dy: The shift in y
    :param constant: The constant to replace rolled values with
    :return: The shifted array with "constant" where roll occurs
    """
    shifted_data = np.roll(data, dx, axis=1)
    if dx < 0:
        shifted_data[:, dx:] = constant
    elif dx > 0:
        shifted_data[:, 0:dx] = constant

    shifted_data = np.roll(shifted_data, dy, axis=0)
    if dy < 0:
        shifted_data[dy:, :] = constant
    elif dy > 0:
        shifted_data[0:dy, :] = constant
    return shifted_data
    
####################################################################################################################################################    
    
def acf_calc(image, dx, dy, beamarea_pix=None):
    """
    dx pixel shift for x axis
    dy pixel shift for y axis
    return (1) ACF at shift (dx, dy)
           (2) standared error on measured noise ACF if beamarea_pix is given.
    """
    acfimage=shift_2d_replace(image, dx, dy, constant=np.nan)*image
    if np.any(beamarea_pix):
        return np.nanmean(acfimage), np.nanstd(acfimage)/np.sqrt(acfimage.size/beamarea_pix)
    else:
        return np.nanmean(acfimage)  
        
####################################################################################################################################################

def mk_acf(noiseimage, pixelhwidth=30, filedir="", filename="", beamarea_pix=None, cpus2use=8):
    """
    Compute noise ACF
    pixelhwidth: required size of noise ACF (half of the width) in pixel, producing (2xpixelwidth+1) x (2xpixelwidth+1) noise ACF
                 note that we need ACF jsut include all the relative vector (in pixel) of interested aperture when computing the variance
                 associated with the integrated flux or spectrum over the aperture
    beamarea_pix: optional parameter required to compute standard error in ACF.
    output: noise ACF 
            if beamarea (in pixel) is given, it compute symmetrized standard error in ACF following Eq.(5) 
            np.array([noiseACF, symmetrized standard error in ACF])
    """
    t1=time.time()
    if isinstance(pixelhwidth, int):
        if len(noiseimage.shape)==3:
            if np.any(beamarea_pix):
                Acf_array=np.zeros([noiseimage.shape[0], int(pixelhwidth*2+1),int(pixelhwidth*2+1),2])
            else:
                Acf_array=np.zeros([noiseimage.shape[0], int(pixelhwidth*2+1),int(pixelhwidth*2+1)])

        if len(noiseimage.shape)==2:
            if np.any(beamarea_pix):
                Acf_array=np.zeros([int(pixelhwidth*2+1),int(pixelhwidth*2+1),2])
            else:
                Acf_array=np.zeros([int(pixelhwidth*2+1),int(pixelhwidth*2+1)])
    else:
        raise Exception("put integer for pixelhwidth, which is required ACF size to compute")

    print("output noise ACF array shape is ",Acf_array.shape)

    if filename and os.path.exists(filedir+"noise_acf_"+filename+"_"+str(pixelhwidth)+".npy"):
        print("file found")
        return np.load(filedir+"noise_acf_"+filename+"_"+str(pixelhwidth)+".npy")

    itest,jtest=np.meshgrid([i for i in range(0,int((Acf_array.shape[1]-1)/2)+1)],[i for i in range(-int((Acf_array.shape[1]-1)/2),int((Acf_array.shape[1]-1)/2)+1)])
    zipped_ji=np.column_stack((jtest.ravel(), itest.ravel()))

    if len(noiseimage.shape)==3:
        if cpus2use:
            #if __name__ == "__main__":
            print('multiprocess is enabled')
            p = mp.Pool(cpus2use)
            result=np.array([list(p.starmap(functools.partial(acf_calc, i, beamarea_pix=beamarea_pix), zipped_ji)) for i in noiseimage])
            #else:
            #    result=np.array([[acf_calc(i, *ij, beamarea_pix=beamarea_pix) for ij in zipped_ji] for i in noiseimage])
        else:
            result=np.array([[acf_calc(i, *ij, beamarea_pix=beamarea_pix) for ij in zipped_ji] for i in noiseimage])

        if np.any(beamarea_pix):
            Acf_array[:,int((Acf_array.shape[1]-1)/2)+itest.ravel(), int((Acf_array.shape[1]-1)/2)+jtest.ravel(),:]=result
            Acf_array[:,int((Acf_array.shape[1]-1)/2)-itest.ravel(), int((Acf_array.shape[1]-1)/2)-jtest.ravel(),:]=result
        else:
            Acf_array[:,int((Acf_array.shape[1]-1)/2)+itest.ravel(), int((Acf_array.shape[1]-1)/2)+jtest.ravel()]=result
            Acf_array[:,int((Acf_array.shape[1]-1)/2)-itest.ravel(), int((Acf_array.shape[1]-1)/2)-jtest.ravel()]=result


    if len(noiseimage.shape)==2:
        if cpus2use:
            #if __name__ == "__main__":
            print('multiprocess is enabled')
            p = mp.Pool(cpus2use)
            result=np.array(list(p.starmap(functools.partial(acf_calc, noiseimage, beamarea_pix=beamarea_pix), zipped_ji)))
            #else:
            #    result=np.array([acf_calc(noiseimage, *ij, beamarea_pix=beamarea_pix) for ij in zipped_ji])
        else:
            result=np.array([acf_calc(noiseimage, *ij, beamarea_pix=beamarea_pix) for ij in zipped_ji])

        if np.any(beamarea_pix):
            Acf_array[int((Acf_array.shape[0]-1)/2)+itest.ravel(), int((Acf_array.shape[0]-1)/2)+jtest.ravel(),:]=result
            Acf_array[int((Acf_array.shape[0]-1)/2)-itest.ravel(), int((Acf_array.shape[0]-1)/2)-jtest.ravel(),:]=result
        else:
            Acf_array[int((Acf_array.shape[0]-1)/2)+itest.ravel(), int((Acf_array.shape[0]-1)/2)+jtest.ravel()]=result
            Acf_array[int((Acf_array.shape[0]-1)/2)-itest.ravel(), int((Acf_array.shape[0]-1)/2)-jtest.ravel()]=result

    if filename and os.path.exists(filedir+"noise_acf_"+filename+"_"+str(pixelhwidth)+".npy")==False:
        np.save(filedir+"noise_acf_"+filename+"_"+str(pixelhwidth), Acf_array)

    t2=time.time()
    print("It took:",t2-t1,"sec to compute the noise ACF.")
    return Acf_array        
    
####################################################################################################################################################
    
def nodiag_view3D(a):
    #https://stackoverflow.com/questions/55588122/how-to-calculate-the-relative-vectors-from-a-list-of-points-from-one-point-to-e
    m = a.shape[0]
    p,q,r = a.strides
    return np.lib.stride_tricks.as_strided(a[:,1:], shape=(m-1,m,2), strides=(p+q,q,r))

####################################################################################################################################################

def mk_relvec(mask):
    """
    return a set of relative position vectors for pairs of True pixel in the input masked array.
    for example, given two True pixel it returns two relative position vector 1->2, 2->1
    masked array: 
    return array with size of N(N-1) of dx and dy
    """
    ytest,xtest=np.where((mask)==1)
    atest=np.column_stack([xtest, ytest])
    testd = (atest-atest[:,None,:])
    return nodiag_view3D(testd).reshape(-1,atest.shape[1])[:,0], nodiag_view3D(testd).reshape(-1,atest.shape[1])[:,1] # dx, and dy
    
####################################################################################################################################################    
    
def mk_noise_var_noiseACF(mask, acf, split=False):
    """
    compute variance (std**2) of noise in integrated flux within aperture following Eq.().
    mask 2d boolean mask, where True represent the aperture region.
    acf 2d noise auto-correlation function
    """
    dx,dy=mk_relvec(mask)                                                  # generate a set of relative position vectors from mask image
    variance=acf[int((len(acf)-1)/2),int((len(acf)-1)/2)]*(mask).sum()     #  variance term acf(0,0) position times the number of pixels
    covariance=acf[dy+int((len(acf)-1)/2),dx+int((len(acf)-1)/2)].sum()   #  covariance term acf(dx,dy) for relative position vectors 
    if split:
        return variance+covariance, variance,covariance
    else:
        return variance+covariance   
        
####################################################################################################################################################        
        
def mk_intflux(mask,img):
    """
    return integrated flux, sum of the pixels where mask==True
    """
    return img._data[mask].sum()
    
####################################################################################################################################################    
    
def mk_aperture(img, rad, plot=False):
    """
    img: 2dimage
    rad: radius of the circular aperture in pix. 
    generate circular aperture
    """
    if len(img.shape)==2:
        YY, XX = np.indices(img.shape[:], dtype='float')  
        Cx=img.shape[0]/2
        Cy=img.shape[1]/2

    Radius = ((YY-Cy)**2 + (XX-Cx)**2)**0.5  
    if len(img.shape)==2:
        Mask = (Radius <= rad)
        if plot:
            plt.imshow(Mask)
            plt.show()
            plt.title("generated aperture")
        return Mask
    else:
        raise Exception("this works only for 2d image")
        
####################################################################################################################################################        
        
def mk_sigmasqrtnumbeam(noise_image, mask, beamarea):
    """
    noise in the integrated flux estimated by sigma_N*\sqrt(N_beam)
    """
    return np.sqrt(mask.sum()/beamarea)*np.nanstd(noise_image)
    
####################################################################################################################################################    
    
def mk_noisespec(mask, acf_cube):
    """
    return underlying noise (1sigma) associated with the spectrum within the given aperture 
    mask: aperture
    acf_cube: 3d array with (velocity, y, x)
    """
    cubemask=np.array([mask for i in range(0,acf_cube[:,:,:].shape[0])])
    noisespec=np.array([np.sqrt(mk_noise_var_noiseACF(mask, acf_cube[i,:,:])) for i in range(0, acf_cube[:,:,:].shape[0])])
    return noisespec
    
####################################################################################################################################################    
    
def mk_noise_var_randaperture(noise_image, mask, plot=True):
    """
    generate noise in spatially integrated flux (2d image) or spectrum (cube) 
    by placing aperture in the noise region randomly and taking the rms across the summed values on the aperture. 
    noise_image: noise_map produced by mk_noisemap
    mask: 2d boolean aperture mask 
    """
    rng=np.random.default_rng()
    if len(noise_image.shape)==3:
        rms=np.copy(noise_image[:,0,0]*0)
        independent_data=np.copy(noise_image[:,0,0]*0)
        cp=np.copy(noise_image)
        convolved_cp=np.copy(noise_image)*0
        for i,j in enumerate(cp):
            convolved_cp[i]=ndimage.convolve(j, mask, mode="constant",cval=np.nan)
            rms[i]=np.nanstd(convolved_cp[i],ddof=1)
            independent_data[i]=np.sum(~np.isnan(convolved_cp[i]))/np.nansum(mask)
            if plot==True:
                plt.imshow(convolved_cp[i])
                plt.colorbar()
                plt.title('summed maps by the aperture')
                plt.show()
                plt.imshow(mask)
                plt.title('aperture shape')
                plt.colorbar()
                plt.show()

        return rms, rms/np.sqrt(2*independent_data-2) #standard deivation and standard error

    if len(noise_image.shape)==2:
        cp=np.copy(noise_image)
        convolved_cp=ndimage.convolve(cp, mask, mode="constant",cval=np.nan)
        rms=np.nanstd(convolved_cp,ddof=1)
        independent_data=np.sum(~np.isnan(cp))/np.nansum(mask)
        if plot==True:
            plt.imshow(convolved_cp)
            plt.title('summed maps by the aperture')
            plt.colorbar()
            plt.show()
            plt.imshow(mask)
            plt.title('aperture shape')
            plt.colorbar()
            plt.show()
            plt.hist(convolved_cp[~np.isnan(convolved_cp)].ravel(), bins='auto')
            plt.title('histgram of data points')
            plt.show()
        return rms, rms/np.sqrt(2*independent_data-2)    
        
####################################################################################################################################################

def mk_simnoise(pixwidth, acf, silent=False):
    """
    simulate noise map (pixwidth x pixwidth) from measured acf
    """
    t1=time.time()
    X,Y = np.meshgrid(np.arange(pixwidth),np.arange(pixwidth))
    # Create a vector of cells
    #XY = np.column_stack((np.ndarray.flatten(X),np.ndarray.flatten(Y)))
    X=np.ndarray.flatten(X)
    Y=np.ndarray.flatten(Y)
    # Calculate a matrix of relative distance vector between the cells
    X = X.reshape(-1,1)
    Xdist = X.T - X+int((acf.shape[0]-1)/2)
    Y = Y.reshape(-1,1)
    Ydist = Y.T - Y+int((acf.shape[0]-1)/2) #acf relative distance vectorの(0,0)
    cov =  np.copy(acf[Ydist,Xdist])
    # covariance matrix to be input for the scipy multivariate_normal
    noise = scipy.stats.multivariate_normal.rvs(mean = np.zeros(pixwidth**2),
        cov = cov)
    noisemap=noise.reshape((pixwidth,pixwidth))
    if silent==False:
        print("It took",str(time.time()-t1),'second to generate '+str(pixwidth)+'x'+str(pixwidth)+' noise map')
    return noisemap
    
####################################################################################################################################################    
    
def mk_simnoise_cube(pixwidth, acf, cpus2use=None):
    """
    simulate noise cube from measured acf
    """
    t1=time.time()
    if len(acf.shape)==3:
        if cpus2use:
            #if __name__ == "__main__":
            print("multiprocess is enabled")
            p = mp.Pool(cpus2use)
            result=np.array(list(p.map(functools.partial(mk_simnoise, pixwidth, silent=True), list(acf))))
            #else:
            #    result=np.array([mk_simnoise(pixwidth, i, silent=True) for i in list(acf)])
        else:
            result=np.array([mk_simnoise(pixwidth, i, silent=True) for i in list(acf)])
        print("It took",time.time()-t1,'second to generate '+str(acf.shape[0])+'x'+str(pixwidth)+'x'+str(pixwidth)+' (v, x, y) noise cube')      
        return result

    else:
        print("use mk_simnoise")

###################################################################################################################################################

def mk_cov(pixwidth, acf, silent=False):
    """
    compute covariance from measured acf
    the output covariance has dimension of pixwidth^2 x pixwidth^2
    """
    t1=time.time()
    X,Y = np.meshgrid(np.arange(pixwidth),np.arange(pixwidth))
    # Create a vector of cells
    #XY = np.column_stack((np.ndarray.flatten(X),np.ndarray.flatten(Y)))
    X=np.ndarray.flatten(X)
    Y=np.ndarray.flatten(Y)
    # Calculate a matrix of relative distance vector between the cells
    X = X.reshape(-1,1)
    Xdist = X.T - X+int((acf.shape[0]-1)/2)
    Y = Y.reshape(-1,1)
    Ydist = Y.T - Y+int((acf.shape[0]-1)/2) #acf relative distance vector (0,0)
    cov =  np.copy(acf[Ydist,Xdist])
    return cov
