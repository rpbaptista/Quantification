# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:56:38 2019

@author: RP258738
"""
import numpy as np

def pondT1(T1,TR,FA):
    FA_rad = np.radians(FA)
    num = (1 - np.exp(-TR/T1))*np.sin(FA_rad)
    den = (1 - np.cos(FA_rad)*np.exp(-TR/T1))
    return num/den

def fwhm2sigma(fwhm):
    return fwhm / np.sqrt(8 * np.log(2))