# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 14:55:55 2019

@author: Renata Porciuncula Baptista
@mail: renata.porciunculabaptista@cea.fr

between this script and computing_applying_tsc_model,
The user must segment the sodium image
     
Steps: 
    1 - Segment sodium image
Goal is from CSF masks and TSC map generate the corrected version.

# WARNING: Pickle is not securized library, make sure no one has access to
#  your pickle file other than u

version: these

"""
# TODO: take into account the PSF

# ------------------------------- IMPORTS ------------------------------------
# Imports : system libraries
import sys
import os

# Imports : python libraries
import nibabel as nib
import matplotlib.pyplot as plt
import pickle
import numpy as np
import scipy.ndimage


# this is the folder where the script is install,
# to be modify if u are running in other coputing
#root = 'Z:/people/Renata/Python_scripts/'
#folder = 'Quantification/'

root = 'C:/Users/rp258738/Documents/Codes_local/'
folder = 'Quantification/src/'

packages = ['libs']

for i in range(len(packages)):
    sys.path.insert(0, os.path.join(root, folder, packages[i]))
    
# Import : homemade libraries
from equations import pondT1, fwhm2sigma

# ------------------------------ CONSTANTS, PATHS ----------------------------
T1_CSF = 70 # [ms]
T1_tissue = 35
C_CSF = 135 # [mM]
spatialResolution = 3 #[mm]
spatialResolutionVoxel = 1 # 
sigma = fwhm2sigma(spatialResolution) # give the sigma necessary to this resolution
sigmaVoxel = fwhm2sigma(spatialResolutionVoxel) # give the sigma necessary to this resolution


FA = 55
TR = 20


#data_root = 'Z:/people/Renata/Data/Pilot_study_1/'
#sub_path = 'sub-{0:02d}/func/'
#tasks_path = ['visual_sodium/', 'motor_sodium/']
#TSC_filename = 'Z:/people/Renata/Data/Pilot_study_1/sub-03/func/visual_sodium/wrSandro_Rec_XMRIDynamic_visual_TR30-sub03_OFF_echo_0_TSC_.nii'
#output_filename =  TSC_filename[0:-4] + 'correct.nii'


sodium_filename = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub01_ab160146/registered_averageFA53.nii',
                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub02_ab160281/registered_averageFA55.nii',
                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub03_ag160127/registered_averageFA55.nii']


CSF_filename = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub01_ab160146/c3anat.nii',
                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub02_ab160281/c3anat.nii',
                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub03_ag160127/c3anat.nii']
WM_filename = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub01_ab160146/c2anat.nii',
                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub02_ab160281/c2anat.nii',
                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub03_ag160127/c2anat.nii']
GM_filename = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub01_ab160146/c1anat.nii',
                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub02_ab160281/c1anat.nii',
                 'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/Results/sub03_ag160127/c1anat.nii']

Machine3T = False



# -------------------------------- Main code ---------------------------------
for idx_img in range(len(sodium_filename)): 
    output_filename = sodium_filename[idx_img][0:-4] + '_TSC.nii'
    output_filename_no_csf =  sodium_filename[idx_img][0:-4] + '_TSC_no_csf.nii'
    
   
    # Load model 
    with open( 'model_{0}_FA{1}.pickle'.format(idx_img,FA), 'rb') as handle:
        model_load = pickle.load(handle) 
    
        # f(x) = ax + b
    a = model_load['a']
    b = model_load['b']
    print("Model:", a,b)
    # Load CSF
    niiData = nib.load(CSF_filename[idx_img])   
    affine = niiData.affine
                 
    # Get dimensions and data
    N, M, P = (niiData.get_fdata()).shape 
    probability_CSF = niiData.get_fdata() 
    
    # Load WM and GM masks
    if Machine3T == False:
        niiData = nib.load(GM_filename[idx_img])                
        probability_GM = niiData.get_fdata()
        niiData = nib.load(WM_filename[idx_img])                
        probability_WM = niiData.get_fdata()
        probability_tissue = probability_GM + probability_WM
    else:
        print("TO IMPLEMENT")
      #  niiData = nib.load(brain_filename)                
      #  probability_tissue = niiData.get_fdata() /np.max(np.max(probability_tissue))
      #  probability_CSF = niiData.get_fdata() /np.max(np.max(probability_CSF)) 
    
    
    # Load native image
    niiData = nib.load(sodium_filename[idx_img])                
    sodium_image = niiData.get_fdata()
    sodium_image[np.isnan(sodium_image)] = 0

    shape_sodium_image = sodium_image.shape
    if len(shape_sodium_image) != 3:
        sodium_image = sodium_image[:,:,:,0]
    
    # Change TSC
    # reducing mask resolution
    maskCSF_T1 = sodium_image*probability_CSF/(pondT1(T1_CSF,TR,FA)) # Signal -> Signal/pondT1 
    #maskCSF_lr = scipy.ndimage.filters.gaussian_filter(probability_CSF.astype(float),
    #                                                 (sigma,sigma,sigma))
    #maskCSF_lr = scipy.ndimage.filters.gaussian_filter(probability_CSF.astype(float),
    #                                                 (sigmaVoxel,sigmaVoxel,sigmaVoxel))
    maskTissue_T1 = sodium_image*probability_tissue/(pondT1(T1_tissue,TR,FA)) # Signal -> Signal/pondT1 
    
    TSC_CSF_apparent = maskCSF_T1*a + b #signal/ponT1 -> mM
    #TSC_CSF_theorical_psf = maskCSF_lr*C_CSF
    TSC_CSF_theorical = probability_CSF*C_CSF
    TSC_tissue_apparent_no_csf =  maskTissue_T1*a - maskCSF_T1*a
    TSC_tissue_apparent =  maskTissue_T1*a + b#+ maskCSF_T1*a
    
    #TSC_tissue_apparent_no_csf[TSC_tissue_apparent < 0] = 0 # concentration negative pas de sense
    
    for i in range(N):
        diff = TSC_CSF_apparent[i,:,:] - TSC_CSF_theorical[i,:,:]
      #  if (np.max(np.max(diff)) - np.min(np.min(diff)))> 0.1:
     #   print("Slice", i, "from", M)
    
#        f, axarr = plt.subplots(1,3,figsize=(10,30))
#        axarr[0].imshow(TSC_CSF_apparent[:,i,:].T, cmap='hot', vmin=0, vmax=150)
#        axarr[1].imshow(TSC_CSF_theorical[:,i,:].T, cmap='hot', vmin=0, vmax=150)
#        axarr[2].imshow((diff.T>0)*150, cmap='hot', vmin=000, vmax=150)
#        
#        print("  MIN:",np.min(np.min(diff)))
#        print("  MAX:",np.max(np.max(diff)))
#        print("  Negative voxels:",np.sum(np.sum(diff<0)))
#        print("  Total pixels:", np.sum(np.sum(probability_CSF[:,:,i] > 0)))        
#     #   plt.colorbar()
#            
#        plt.show()
#            
    
    for i in range(N):
        if np.sum(np.sum(TSC_tissue_apparent[i,:,:].T)) != 0:
            print("Slice",i,"from", N)
            
            plt.imshow(TSC_tissue_apparent[i,:,:], cmap='hot', vmin=000,vmax=60)
            plt.colorbar()
            plt.show()
    
    
     #--------------------------- Save TSC
        nib.save(nib.Nifti1Image(TSC_tissue_apparent, affine), output_filename )
        nib.save(nib.Nifti1Image(TSC_tissue_apparent_no_csf, affine), output_filename_no_csf )