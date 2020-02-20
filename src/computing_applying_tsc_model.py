# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 10:45:13 2019

@author: Renata Porciuncula Baptista
@email: renata.porciunculabaptista@cea.fr

Goal of this script is quantify TSC in sodium images.
This is done by a 4 point calibration, where the intensities are compensated by
T1 relaxations.

Note: this script only works if [:,:,i] be a transversal slice of tubes

f(Signal/pT1) - > Concentration mm/L/voxel


"""
# ------------------------------- IMPORTS
# Imports : system libraries
import sys
import os

# Imports : python libraries
import cv2
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from dipy.denoise.nlmeans import nlmeans
from dipy.denoise.noise_estimate import estimate_sigma
import pickle
from pylab import * 
from scipy.stats import linregress


# this is the folder where the script is install,
# to be modify if u are running in other coputing
root = 'Z:/people/Renata/Python_scripts/'
folder = 'Quantification/'
packages = ['libs']

for i in range(len(packages)):
    sys.path.insert(0, os.path.join(root, folder, packages[i]))
    
# Import : homemade libraries
from segmentation import getTube3D
from equations import pondT1

# ------------------------------ Constants, paths
# Path strings
data_root = 'Z:/people/Renata/Data/Pilot_study_1/'
sub_path = 'sub-{0:02d}/func/'
tasks_path = ['visual_sodium/', 'motor_sodium/']
filename = 'Sandro_Rec_XMRIDynamic_visual_TR30-sub{0:02d}_{1}_echo_{2}.nii'
output_filename =  filename[0:-4] + '_TSC_.nii'

#filename = 'Z:/people/Renata/Data/Test/rfunc.nii'
# Constants
sub_array = [3]
states = ['OFF','ON']
#TR = 30 #[ms]
#FA = 59

#alexa data
TR = 20 #
FA = 25

# Tubes - article LEROI 2018
b_agar_model = 52
a_agar_model = (36-52)/5

number_tubes = 2  # right left subdivide in RB RF LB LF where F:front B: back

# information on tubes _7T
tubes_concentrations_front = [105, 209]
tubes_concentrations_back = [51, 155]

# 3T 
tubes_concentrations_front = [51,105]
tubes_concentrations_back = [51,105]
Machine3T = True

agar_concentrations = np.asarray([2.0, 2.0, 2.0, 2.0]) #np.asarray changes type 
#list to np array, which is optimize to compute operation in all array at once
factor_correction_liquid_tissue = 1#0.8
factor_psf = 0.8

T1_agar_tubes = agar_concentrations*a_agar_model + b_agar_model
# Once T1_agar_tubes is known; comment the line above and replace by
T1_agar_tubes = np.asarray([36.2, 36.2, 36.2, 36.2])
#T2 = 13.4 or 17/5.2 #T2star = 4.3
            
# Pick a file to get the size
filename_ex = data_root + sub_path.format(sub_array[0]) + tasks_path[0] + filename.format(sub_array[0], states[0], 0)
filename_ex = 'Z:/people/Alexa/B0951_Na/average25.nii'
# Load imaging / to open the file
mean = False
try:
    niiDataEx = nib.load(filename_ex)
    # Dimensions
    shape_total = (niiDataEx.get_fdata()).shape
    print(shape_total)
    if len(shape_total) == 3:
        N, M, P = (niiDataEx.get_fdata()).shape
        FourDimensions = False

    else:
        N, M, P, Q = (niiDataEx.get_fdata()).shape
        FourDimensions = True
except IOError:
    print('An error occured trying to read the file:'+ filename_ex)
    print('Verify the path or the extension .nii .niigz')


# Model variables
a = 0
b = 0

# Outputs
average_intensities_front = np.zeros(number_tubes)
average_intensities_back = np.zeros(number_tubes)

concentrations_maps = np.zeros((N, M, P))

# -------------------------------- Main code
for i in range(len(sub_array)):
    
    #------------------------- Getting masks
   # filename_current = data_root + sub_path.format(sub_array[i]) + tasks_path[0] + filename.format(sub_array[i], states[0], 0)
   # output_current = data_root + sub_path.format(sub_array[i]) + tasks_path[0] + output_filename.format(sub_array[i], states[0], 0)
    filename_current = 'Z:/people/Alexa/B0951_Na/average25.nii'
    
   # Try to open the file
    try:
        niiData = nib.load(filename_current)
    except IOError:
        print('An error occured trying to read the file:'+ filename_ex)
        print('Verify the path or the extension .nii .niigz')
   
    # Retrieve the images
    image = niiData.get_fdata()
    if FourDimensions:
        image = image[:,:,:,0]
        image = np.swapaxes(image,0,1)
    affine = niiData.affine

    
    # In order to call ``non_local_means`` first you need to estimate the standard
    # deviation of the noise. We use N=4 since the Sherbrooke dataset was acquired
    # on a 1.5T Siemens scanner with a 4 array head coil.
    sigmaB = estimate_sigma(image, N=32)
    
    
    # Calling the main function ``non_local_means``
    den = nlmeans(image, sigma=sigmaB, patch_radius= 1, block_radius = 1, rician= True)

    # OTSU only uses image8
    image8 = np.array(den/np.max(np.max(den)) * 255.0, dtype = np.uint8)
    image64 = image#den * 255 # eauivalente
    plt.imshow(image64[:,20,:], cmap='gray')

    masks, idx_change = getTube3D(image8)
    
    # veryfing if its coherent the changing tubes
    if idx_change == -1:
        print("ERROR: automatic slice changing tubes FAILED\n Veryfing manually the slice where the tubes changes, not")
        if Machine3T != True:
            exit(-1)
        else:
            idx_change = P - 1
        
    for i in range(P):
        print("z",i)
        plt.imshow(np.hstack([image64[:,:,i].T,
                              masks[1,:,:,i].T/np.max(np.max(image8[:,:,i])),
                              masks[0,:,:,i].T/np.max(np.max(image8[:,:,i]))]), cmap='gray')
        plt.show()
       
    # ------------------- Keeping only center
    
    kernel = np.ones((3,3), np.uint8)       # set kernel as 3x3 matrix from numpy
    #Create erosion and dilation image from the original image
    erosion_masks = np.zeros(masks.shape)
   # erosion_masks = cv2.erode(masks, kernel, iterations=1)
    for k in range(P):
        for j in range(number_tubes):
            erosion_masks[j,:,:,k] = cv2.erode(masks[j,:,:,k], kernel, iterations=1)
        
    #erosion_masks = masks    
    # ------------------- Analyzing histogram
    
    # retrive regions where we found the tubes
    #idx_true_mask =  np.where(mask_current == 1)
    for i in range(number_tubes):
       # masks = erosion_masks
        # Saving 0 - idx change, in Z
        idx_true_mask =  np.where(erosion_masks[i,:,:,0:idx_change] > 0)
        roi_tube = image64[idx_true_mask]
        average_intensities_front[i] = np.mean(roi_tube)
        
        if Machine3T != True:
            # Saving values from idx change - Last
            tst_b = masks[i,:,:,idx_change:-1]
            idx_true_mask =  np.where(erosion_masks[i,:,:,idx_change:-1] > 0)
            roi_tube = image64[idx_true_mask]
            average_intensities_back[i] = np.mean(roi_tube)

        
        # IDENTIFYING the order
        tst_f = erosion_masks[i,:,:,0:idx_change]
        
        idx_true_mask =  np.where(tst_f > 0)
        roi_tube = image64[idx_true_mask]
        plt.hist(x = roi_tube)
        plt.title("Histogram of front, tube:" + str(i))
        plt.show()
        
        print("G max front, ",np.mean(roi_tube))
        if Machine3T != True:
            tst_b = erosion_masks[i,:,:,idx_change:-1]
            
            idx_true_mask =  np.where(tst_b > 0)
            roi_tube = image64[idx_true_mask]
            fig = plt.hist(x = roi_tube)
            plt.title("Histogram of back, tube:" +str(i))
            plt.show()
        
            print("G max back, ",np.mean(roi_tube))
          


    # -------------------------- Quantifying 
    if Machine3T != True:
        y = np.hstack([tubes_concentrations_front, tubes_concentrations_back])
        x = np.hstack([average_intensities_front, average_intensities_back])
    else:
        y = tubes_concentrations_front
        x = average_intensities_front
        T1_agar_tubes = T1_agar_tubes[0:2]
        
    factor_correction_T1_tubes  = pondT1(T1_agar_tubes, TR, FA)
    x_correct = (x*factor_correction_liquid_tissue)/(factor_correction_T1_tubes*factor_psf)

    gradient, intercept, r_value, p_value, std_err = linregress(x_correct,y)

    a = gradient
    b = intercept
    
    
    # ---------------------------- Veryfing 
    plt.plot(x_correct, y, 'yo', x_correct, a*x_correct+b, '--k') 
    plt.show() 
    print("R_value:", r_value)

    
    # ---------------------------- Save model
    model = {'a' : a, 'b' : b, 'std_error': std_err, 'TR': TR, 'FA': FA}
    
    with open('model.pickle', 'wb') as handle:
        pickle.dump(model, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open('model.pickle', 'rb') as handle:
        model_load = pickle.load(handle)
    

 #  #--------------------------- Save TSC
  #  nib.save(nib.Nifti1Image(image64*a+b, affine), output_current)

   #_____________________
   
'''
 References:
     
 [Leroi 2018] - Simultaneous multi-parametric mapping of total sodium concentration,
                 T1,T2and ADC at 7 T using a multi-contrast unbalanced SSFP
                 
 [Reetz 2012] - Increased brain tissue concentration in Huntington's Disease _
                 A sodium imaging study at 4T
'''  
    