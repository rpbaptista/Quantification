# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:58:48 2020

@author: RP258738
"""

# Create syntenthic masks
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

packages = ['libs','class']

for i in range(len(packages)):
    sys.path.insert(0, os.path.join(root, folder, packages[i]))
# Import : homemade libraries
from equations import pondT1, pondT2, fwhm2sigma
from model import FourCompartmentModel

# %%    
# Data
CSF_filename = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AB160146_AB160146-1948/ab160146/1_013_t1_mp2rage_sag_iso0_75mm_UNI-DEN_20190606/c3ab160146_20190606_001_013_t1_mp2rage_sag_iso0_75mm_t1_mp2rage_sag_iso0_75mm_UNI-DEN.nii',
                'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AB160281_AB160281-2081/ab160281/1_015_t1_mp2rage_sag_iso0_75mm_UNI-DEN_20191209/c3ab160281_20191209_001_015_t1_mp2rage_sag_iso0_75mm_t1_mp2rage_sag_iso0_75mm_UNI-DEN.nii',
                'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AG160127_AG160127-2015/ag160127/1_013_t1_mp2rage_sag_iso0_75mm_UNI-DEN_20190923/c3ag160127_20190923_001_013_t1_mp2rage_sag_iso0_75mm_t1_mp2rage_sag_iso0_75mm_UNI-DEN.nii'];
WM_filename = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AB160146_AB160146-1948/ab160146/1_013_t1_mp2rage_sag_iso0_75mm_UNI-DEN_20190606/c2ab160146_20190606_001_013_t1_mp2rage_sag_iso0_75mm_t1_mp2rage_sag_iso0_75mm_UNI-DEN.nii',
                'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AB160281_AB160281-2081/ab160281/1_015_t1_mp2rage_sag_iso0_75mm_UNI-DEN_20191209/c2ab160281_20191209_001_015_t1_mp2rage_sag_iso0_75mm_t1_mp2rage_sag_iso0_75mm_UNI-DEN.nii',
                'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AG160127_AG160127-2015/ag160127/1_013_t1_mp2rage_sag_iso0_75mm_UNI-DEN_20190923/c2ag160127_20190923_001_013_t1_mp2rage_sag_iso0_75mm_t1_mp2rage_sag_iso0_75mm_UNI-DEN.nii'];
GM_filename = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AB160146_AB160146-1948/ab160146/1_013_t1_mp2rage_sag_iso0_75mm_UNI-DEN_20190606/c1ab160146_20190606_001_013_t1_mp2rage_sag_iso0_75mm_t1_mp2rage_sag_iso0_75mm_UNI-DEN.nii',
                'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AB160281_AB160281-2081/ab160281/1_015_t1_mp2rage_sag_iso0_75mm_UNI-DEN_20191209/c1ab160281_20191209_001_015_t1_mp2rage_sag_iso0_75mm_t1_mp2rage_sag_iso0_75mm_UNI-DEN.nii',
                'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AG160127_AG160127-2015/ag160127/1_013_t1_mp2rage_sag_iso0_75mm_UNI-DEN_20190923/c1ag160127_20190923_001_013_t1_mp2rage_sag_iso0_75mm_t1_mp2rage_sag_iso0_75mm_UNI-DEN.nii'];
               
               
output_filename = ['C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AB160146_AB160146-1948/synthetic_sodium.nii',
                     'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AB160281_AB160281-2081/synthetic_sodium.nii',
                'C:/Users/rp258738/Documents/VirtualMachine-Ubuntu/2020-03-24/Reconstructed/AG160127_AG160127-2015/synthetic_sodium.nii'];


TR = 20
FA = 55 
TE = 5 

# MODELED_ USED IN THE ON/OFF trial
fraction_T2short = 0.4 

# Table article Giles; Madelin 2017 :
# Multipulse sodium magnetic resonance imaging
#for multicompartment quantification: proof_of_concept
 
           # T1 T2S T2L
ISC_times = [24, 2,   12] 
ESC_times = [46, 3.5, 30] 
CSF_times = [64, 56,  56]

# Water part
wGM = 0.85
wWM = 0.7
cISC = 12 

# aplha
alphaISC = 0.70  
alphaESC = 0.30
alphaCSF = 1.00
# https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.22133
# Computing
model4C = FourCompartmentModel(TR, FA, 
                               [ISC_times[0], ESC_times[0], CSF_times[0]],
                               [alphaISC, alphaESC, alphaCSF], [cISC])

for i in range(len(output_filename)):
    affine = (nib.load(CSF_filename[i])).affine
    prob_csf = (nib.load(CSF_filename[i])).get_fdata() 
    prob_gm = (nib.load(GM_filename[i])).get_fdata() 
    prob_wm = (nib.load(WM_filename[i])).get_fdata() 
    
    imageCSF = prob_csf * model4C.get_CSF(TE,CSF_times[1],CSF_times[2],fraction_T2short)
    imageGM = wGM * prob_gm *( model4C.get_ISC(TE,ISC_times[1],ISC_times[2],fraction_T2short) +
                                     model4C.get_ESC(TE,ESC_times[1],ESC_times[2],fraction_T2short))
    imageWM = wWM * prob_wm *( model4C.get_ISC(TE,ISC_times[1],ISC_times[2],fraction_T2short) +
                                     model4C.get_ESC(TE,ESC_times[1],ESC_times[2],fraction_T2short))
    
    TSC = imageCSF + imageGM + imageWM
    nib.save(nib.Nifti1Image(TSC, affine), output_filename[i] )
    print(i)                    