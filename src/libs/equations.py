# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:56:38 2019

@author: RP258738
"""
from math import acos, degrees, sin, cos, radians
import numpy as np
    
def pondT1(T1,TR,FA):
    FA_rad = np.radians(FA)
    num = (1 - np.exp(-TR/T1))*np.sin(FA_rad)
    den = (1 - np.cos(FA_rad)*np.exp(-TR/T1))
    return num/den
def pondT2(f, TE, T2_s, T2_l):
    return f*np.exp(-TE/T2_s)+(1-f)*np.exp(-TE/T2_l)
def fwhm2sigma(fwhm):
    return fwhm / np.sqrt(8 * np.log(2))


  

def getPondT1(FA, TR, T1):
    angleRad = radians(FA)
    return sin(angleRad)*(1 - np.exp(-TR/T1))/(1 - cos(angleRad)*np.exp(-TR/T1))


def biexpSimplified (TE,a1l,pT1):

           # T1 T2S T2L
    IC_times = [24, 2,   12]
    EC_times = [46, 3.5, 30]

    f = 0.4

    ISC_var = -pT1*(a1l*concMassConservation(a1l) - 12*0.7) * biexp (TE,IC_times[1],IC_times[2],f )
    ESC_var = -pT1*(1-a1l-0.04-0.3)*140 * biexp (TE,EC_times[1],EC_times[2],f)

    return ISC_var+ ESC_var

def KfromSNR (SNR, sigmaNoise = 1.0):
    # SNR = 20 log (SigmaSig/SigmaNoise)
    sigmaSignal = sigmaNoise * np.exp(SNR/20)
    return np.sqrt(sigmaSignal)

def concMassConservation(a1l):
    C1 = 12
    C2 = 140
    a1 = 0.7
    bloodFractionRepos = 0.03
    bloodFractionAct = 0.04
    return (( a1*C1 + (1-a1-bloodFractionRepos)*C2- (1-a1l-bloodFractionAct)*C2)/a1l)

def concentrationMassConservation(a1, a1l, C1, C2 = 140, bloodFractionRepos = 0.03, bloodFractionAct = 0.04):
    """Return the concentration(s) for a given alpha1, alpha1line; C1; C2"""
    # alpha1C1 + (1-alpha1)C2 = alpha1lineC1line + (1-alpha1line)C2
    # Thought to be used, C1 ISC, C2 ESC

    if type(a1) is float and type(a1l) is float :
        return ( a1*C1 + (1-a1-bloodFractionRepos)*C2- (1-a1l-bloodFractionAct)*C2)/a1l
    if type(a1) is float and (type(a1l) is list):
        answer = [None] * len(a1l)
        for i in range(len(a1l)):
            answer[i] = ( a1*C1 + (1-a1-bloodFractionRepos)*C2- (1-a1l[i]-bloodFractionAct)*C2)/a1l[i]
        return answer
    else:
        try:
            return ( a1*C1 + (1-a1-bloodFractionRepos)*C2- (1-a1l-bloodFractionAct)*C2)/a1l
        except:
            print("Type incorrect: a1 (float), a1l (float /list of float)")

def angleErnst(TR,T1):
    """Return the angle (degrees) Ernst for a given TR, T1 (ms)"""
    return  degrees(acos(np.exp(-TR/T1)))

def funcDiffON_OFF(alpha1, C1, alpha1l, C1l):
    return 0

def biexp (TE,T2bS,T2bL,f,K=1):
    return K*(f*np.exp(-1*TE/T2bS) + (1-f)*np.exp(-1*TE/T2bL))

def fwhm2sigma(fwhm):
    return fwhm / np.sqrt(8 * np.log(2))
