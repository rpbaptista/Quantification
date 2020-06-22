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

version: these

"""
# ------------------------------- IMPORTS
# Imports : system librariesim
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
from datetime import datetime
startTime = datetime.now()

# this is the folder where the script is install,
# to be modify if u are running in other coputing
#root = 'Z:/people/Renata/Python_scripts/'
#root = 'C:/Users/rp258738/Documents/Codes_local/'
root = '/neurospin/ciclops/people/Renata/Codes/'
folder = 'Quantification/src/'
output_folder = 'Quantification/output/'

packages = ['libs']

for i in range(len(packages)):
    sys.path.insert(0, os.path.join(root, folder, packages[i]))
    #sys.path.append()
# Import : homemade libraries
from segmentation import getTube3D, createArtificialMask
from equations import pondT1

# ------------------------------ Constants, paths
# Path strings
#filename_list = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub01_ab160146/average_FA53.nii',
#                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub02_ab160281/average_FA55.nii',
#                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub03_ag160127/average_FA55.nii']
#filename_list = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub01_ab160146/average_FA25.nii',
#                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub02_ab160281/average_FA25.nii',
#                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub03_ag160127/average_FA25.nii']
#filename_list = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub01_ab160146/average_FA25.nii',
#                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub02_ab160281/average_FA25.nii',
#                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub03_ag160127/average_FA25.nii']
#filename_list = ['/neurospin/ciclops/people/Renata/ProcessedData/Results/sub01_ab160146/average_FA25.nii',
#                 '/neurospin/ciclops/people/Renata/ProcessedData/Results/sub02_ab160281/average_FA25.nii',
#                 '/neurospin/ciclops/people/Renata/ProcessedData/Results/sub03_ag160127/average_FA25.nii']


#filename_list = ['X:\people\Renata\ReconstructedData\Examples\M0_na.nii.gz']
filename_list = ['/neurospin/ciclops/people/Renata/ReconstructedData/23Na_patients/Examples/M0_na.nii.gz']

#                 
# Constants
#alexa data
TR = 20 #
FA = 55

# Tubes - article LEROI 2018
b_agar_model = 52
a_agar_model = (36-52)/5

number_tubes = 2  # right left subdivide in RB RF LB LF where F:front B: back

# information on tubes _7T
tubes_concentrations_front = [51, 155]
tubes_concentrations_back = [105, 209]

#tubes_concentrations_front= [105, 155]
#tubes_concentrations_back = [51, 209]


# OPTIONS ================= 3T 
Machine3T = False
ArtificialMask = True
radius = 4

if (Machine3T == True):
    tubes_concentrations_front = [51,105]
    tubes_concentrations_back = [51,105]

agar_concentrations = np.asarray([0.3, 0.9, 0.6, 1.2]) #np.asarray changes type 
#list to np array, which is optimize to compute operation in all array at once
factor_correction_liquid_tissue = 1#0.8
factor_psf = 1

T1_agar_tubes = agar_concentrations*a_agar_model + b_agar_model
# Once T1_agar_tubes is known; comment the line above and replace by
#T1_agar_tubes = np.asarray([36.2, 36.2, 36.2, 36.2])
#T2 = 13.4 or 17/5.2 #T2star = 4.3
            
# Pick a file to get the size
# Load imaging / to open the file
mean = False
try:
    niiDataEx = nib.load(filename_list[0])
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
    print('An error occured trying to read the file:'+ filename_list[0])
    print('Verify the path or the extension .nii .niigz')


# Model variables

a = 0
b = 0

# Outputs
average_intensities_front = np.zeros(number_tubes)
average_intensities_back = np.zeros(number_tubes)
median_intensities_front = np.zeros(number_tubes)
median_intensities_back = np.zeros(number_tubes)
percentil_intensities_front = np.zeros(number_tubes)
percentil_intensities_back = np.zeros(number_tubes)

concentrations_maps = np.zeros((N, M, P))


# -------------------------------- Main code
for idx_img in range(len(filename_list)):
    
    #------------------------- Getting image
    filename_current = filename_list[idx_img]
 
   # Try to open the file
    try:
        niiData = nib.load(filename_current)
    except IOError:
        print('An error occured trying to read the file:'+ filename_current)
        print('Verify the path or the extension .nii .niigz')
   
    # Retrieve the images
    image = niiData.get_fdata()
    image[np.isnan(image)] = 0
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

    masks, idx_change, list_c1, list_c2 = getTube3D(image8)
    
    # create artificial mask
    if (ArtificialMask == True):
        artificial_masks = np.zeros((2,M,N,P))
        for p in range(P):
            artificial_masks[0,:,:,p] = createArtificialMask(np.squeeze(image8[:,:,1]),list_c1[p,:],radius);
            artificial_masks[1,:,:,p] = createArtificialMask(np.squeeze(image8[:,:,1]),list_c2[p,:],radius);
         
    
    # veryfing if its coherent the changing tubes
    if Machine3T == False:
        if idx_change == -1:
            print("ERROR: automatic slice changing tubes FAILED\n Veryfing manually the slice where the tubes changes, not")
            idx_change = int(P/2 - 1)
    else:
        if idx_change == -1:
            idx_change = P - 1
               
        
  
      
    # ------------------- Keeping only center
    if (ArtificialMask == False):
        kernel = np.ones((4,4), np.uint8)       # set kernel as 3x3 matrix from numpy
        #Create erosion and dilation image from the original image
        erosion_masks = np.zeros(masks.shape)
       # erosion_masks = cv2.erode(masks, kernel, iterations=1)
        for k in range(P):
            for j in range(number_tubes):
                erosion_masks[j,:,:,k] = cv2.erode(masks[j,:,:,k], kernel, iterations=1)
    else: 
        erosion_masks = artificial_masks
        
    #erosion_masks = masks    
    initial_slice = int(np.round(P*0.12))
    secure_margin = 3
    final_slice =  int(np.round(P*0.88))

    for i in np.arange(initial_slice, final_slice, 1, dtype=np.int16):
      #  print("z",i)
        if np.sum( masks[1,:,:,i] +  masks[0,:,:,i]) > 0 and (i < idx_change-secure_margin or i > idx_change+secure_margin) :
            
            img_scaled = image64[:,:,i].T/np.max(np.max(image64[:,:,i]))
            
            mask2 = masks[1,:,:,i].T/np.max(np.max(image8[:,:,i]))
            mask1 = masks[0,:,:,i].T/np.max(np.max(image8[:,:,i]))
            
            erosion_mask2 = erosion_masks[1,:,:,i].T/np.max(np.max(image8[:,:,i]))
            erosion_mask1 = erosion_masks[0,:,:,i].T/np.max(np.max(image8[:,:,i]))
       
#            plt.imshow(np.hstack([img_scaled,
#                                  erosion_mask1*img_scaled,
#                                  erosion_mask2*img_scaled]), cmap='gray')
#            plt.show()
    # ------------------- Analyzing histogram
    
    # retrive regions where we found the tubes
    #idx_true_mask =  np.where(mask_current == 1)
    for i in range(number_tubes):
       # masks = erosion_masks
        # Saving 0 - idx change, in Z
        part_tube = np.zeros((erosion_masks[i,:,:,:]).shape)
        part_tube[:,:,initial_slice:idx_change-secure_margin] = 1
        idx_true_mask =  np.where(erosion_masks[i,:,:,:]*part_tube > 0)
        #print(size(idx_true_mask))
        roi_tube_front = image64[idx_true_mask]
       # roi_tube_front = np.ma.array(roi_tube_front, mask=roi_tube_front>0.001)
        average_intensities_front[i] = np.mean(roi_tube_front)
        median_intensities_front[i] = np.median(roi_tube_front)
        percentil_intensities_front[i] = np.sort(roi_tube_front)[int(len(roi_tube_front)*0.9)]
        plt.plot(np.sort(roi_tube_front))
        plt.show()
        if Machine3T == False:
            # Saving values from idx change - Last
            part_tube = np.zeros((erosion_masks[i,:,:,:]).shape)
            part_tube[:,:,idx_change+secure_margin:final_slice] = 1
        
            idx_true_mask =  np.where(erosion_masks[i,:,:,:]*part_tube > 0)
            roi_tube_back = image64[idx_true_mask]
        #    roi_tube_back = np.ma.array(roi_tube_back, mask=roi_tube_back>0.001)
         
            average_intensities_back[i] = np.mean(roi_tube_back)
            median_intensities_back[i] = np.median(roi_tube_back)
            percentil_intensities_back[i] = np.sort(roi_tube_back)[int(len(roi_tube_back)*0.9)]
           # plt.plot(np.sort(roi_tube_back))
           # plt.show()


      # Plotting
     #   plt.hist(x = roi_tube_front)
     #   plt.title("Histogram of front, tube:" + str(i))
     #   plt.show()       
     #   print("G mean front, ",np.mean(roi_tube_front))
   
     #   if Machine3T == False:
     #       fig = plt.hist(x = roi_tube_back)
     #       plt.title("Histogram of back, tube:" +str(i))
     #       plt.show()        
     #       print("G mean back, ",np.mean(roi_tube_back))
          


    # -------------------------- Quantifying 
    if Machine3T == False:
        y = np.hstack([tubes_concentrations_front, tubes_concentrations_back])
       # x = np.hstack([average_intensities_front, average_intensities_back])
      #  x = np.hstack([median_intensities_front, median_intensities_back])
        x = np.hstack([percentil_intensities_front, percentil_intensities_back])
    else:
        y = tubes_concentrations_front
        x = average_intensities_front
        T1_agar_tubes = T1_agar_tubes[0:2]
        
    factor_correction_T1_tubes  = pondT1(T1_agar_tubes, TR, FA)
    x_correct = (x*factor_correction_liquid_tissue)/(factor_correction_T1_tubes*factor_psf)
    x_correct = x
    gradient, intercept, r_value, p_value, std_err = linregress(x_correct,y)

    a = gradient
    b = intercept
    
    
    # ---------------------------- Veryfing 
    plt.plot(x_correct, y, 'yo', x_correct, a*x_correct+b, '--k') 
    plt.title("Calibration curve: y = {0}x{1}".format(a,b))
    plt.show() 
    print("R_value:", r_value)

    print(datetime.now() - startTime)

    # ---------------------------- Save model
    model = {'a' : a, 'b' : b, 'std_error': std_err, 'TR': TR, 'FA': FA}
    
    model_filename = 'model_{0}_FA{1}.pickle'.format(idx_img,FA)
    fullpath_output = os.path.join(root, output_folder, model_filename)

    with open( fullpath_output, 'wb') as handle:
        pickle.dump(model, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open( fullpath_output, 'rb') as handle:
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
    