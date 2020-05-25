# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 16:20:26 2019

@author: Renata Porciuncula Baptista
mail: renata.porciunculabaptista@cea.fr

Library:
    Contains class usefull to the sodium quantification modelisation
    4Compartmentmodel is a model 
    it has attributes type compartment
    
    Model: contains attributes that are shared between compartments; i.e;
        sequence parameters (TR,FA), as dictionnary with the compartments,
        a getTSC function that allows to obtain the sum of all compartments
    
    FourCompartmental specialization:
        gives the names of the compartments, standart values for concentration,
        and physical values (T1) for each compartments, as volume fractions;
    
    Compartiments: contains spoiled, and not spoiled ponderations, and the 
    signal model for each compartment, it can be mono or bi exponential
    
    OBS: If u gonna use a model with different number of compartments, we strongly
    recommend to create another class daughter of model. It s likely, u'll instanciate
    several objects from the same class (activation, repos, different tissues),
    this will save time, and make it easier to correct mistakes
"""
#### IMPORTS
import sys
from math import exp, radians,sin,cos
from equations import biexp
import copy

#### CLASS
class Compartment:
    """This class models a compartment in a multicompartment model (Sodium MRI),
    the signal description can be mono or biexponential, spoiled or not.
    
    Its attributes:
        degree: int
            to say if is mono or biexponencial
        concentration: float/int
            determine the concentration of the subtance in the comp
        volFraction: float between 0 and 1
            proportional volume of this compart in the model
        T1: float/int, [ms]
            the relaxation time T1, of this environment
    
    Its public methods:
        function(self, TE, T2_1, T2_2 = -1, f=1): non poderate signal, non spoiled
        functionPond : ponderate signal, spoiled or not
        
         get/set for each attribute named accordingly
         get_attribute
         set_attribute
         
        """
    
    # PRIVATE METHODS
    def __init__(self, degree, volFraction, concentration, T1 = 1):
         self.deg = degree
         self.conc = concentration
         self.volFraction = volFraction
         self.T1 = T1
         
    def __getEta(self, TR):
       return 1 - exp(-TR/self.T1)
   
    def __getSpoiled(self, FA, TR):
        angleRad = radians(FA)
      
        return sin(angleRad)*(1 - exp(-TR/self.T1))/(1 - cos(angleRad)*exp(-TR/self.T1))
   
    # PUBLIC METHODS
    def function(self, TE, T2_1, T2_2, f=1):
        """  
        Return value in a given TE from a monoexponencial or biexponencial model signal
        
        '''
        
        function(self, TE, T2_1, T2_2 , f=1): non poderate signal, non spoiled
            TE = time Echo
            T2_1 = T2star short
            T2_2 = T2star long
            f = proportion of short, if f=1, monoexp
        ."""
       
        if self.deg == 1:
            return biexp (TE,T2_1,0,1)
        
        if self.deg == 2:       
            return biexp (TE,T2_1,T2_2,f)
    
    def functionPond (self,FA, TR, TE, T2_1, T2_2, f=1, spoiled = True):
        """  
        Return value in a given TE from a ponderate (alpha, concentration)
        monoexponencial or biexponencial model signal, spoiled or not
        
        '''
        
        function(self, TE, T2_1, T2_2 = -1, f=1): non poderate signal, non spoiled
            TE = time Echo
            T2_1 = T2star short
            T2_2 = T2star long
            f = proportion of short, if f=1, monoexp
            spoiled (boolean) = spoiled the signal or not

        ."""
        if spoiled == True:

            pond = self.conc * self.volFraction * self.__getSpoiled(FA,TR)
        
        else:
            pond = self.conc * self.volFraction * self.__getEta(TR)
        return pond * self.function(TE,T2_1,T2_2,f)
        
    def set_degree(self, degree):
        """Set mono/biexponencial"""
        if type(degree) is int:
            self.deg = degree
        else:
            print("Degree couldn't be set, type require is int type passed is{0}".format(type(degree)))
            sys.exit(1)        
            
    def get_degree(self):
        """Get degree mono or bi"""
        return self.deg
    
    def set_concentration (self, concentration):
        """Set concentration [mMol]"""
        if type(concentration) is float or type(concentration) is int:
            self.conc = concentration
        else:
            print("Concentration couldn't be set, type require is number (int/float) type passed is{0}".format(type(concentration)))
            sys.exit(1)
            
    def get_concentration(self):
        """Get concentration [mMol]"""
        return self.conc
    
    def set_volFraction(self, volFraction):
        """Set volume fraction, betwenn 0 and 1"""
        if type(volFraction) is float or type(volFraction) is int:
            if volFraction >= 0 and volFraction <= 1:
                self.volFraction = volFraction
            else:
                print("VolFraction needs to be between 0 and 1")
                sys.exit(1)
        else:
            print("VolFraction couldn't be set, type require is number (int/float) type passed is{0}".format(type(volFraction)))
            sys.exit(1)
            
    def get_volFraction (self):
        """Get volume fraction"""
        return self.volFraction
    
    def set_T1(self, T1):
        """Set T1 [ms]"""
        if type(T1) is float or type(T1) is int:    
            self.T1 = T1
        else:
            print("T1 couldn't be set, type require is number (int/float) type passed is{0}".format(type(T1)))
            sys.exit(1)
            
    def get_T1 (self):
        """Get T1 [ms]"""
        return self.T1
            
      
class Model:
    # Definitions
    """This class contais the shared attribute between compartments in a 
       multicompartment model (Sodium MRI),
    
    
        Its attributes:
            FA: float
                Flip angle
            TR: float/int, [ms]
                repetition time of the sequence
            __Compartments: dictionnary of compartments
                private
            
        Its public methods:
            get_TSC sum of each compartment function non weight,
            get/set for each attribute public named accordingly
            get_attribute
            set_attribute
        ."""
        
    def __init__(self, TR = 0, FA = 0):
        self.compartments = {} 
        self.TR = TR
        self.FA = FA
    
    def get_TSC(self, TE, T2_1, T2_2, f):
        """Get total sodium concentration """
        TSC = 0
        for key, value in self.compartments.items():
            TSC = TSC + value.function( TE, T2_1, T2_2, f)
        return TSC
    
    def set_TR(self, TR):
        """Set repetition time (float)"""
        if type(TR) is float or type(TR) is int:    
            self.TR = TR
        else:
            print("TR couldn't be set, type require is number (int/float) type passed is{0}".format(type(TR)))
            sys.exit(1)
            
    def get_TR(self):
        """Get repetition time"""
        return self.TR
    
    def set_FA (self, FA):
        """Set flip angle (degree)"""
        if type(FA) is float or type(FA) is int:    
            self.FA = FA
        else:
            print("FA couldn't be set, type require is number (int/float) type passed is{0}".format(type(FA)))
            sys.exit(1)

    def get_FA(self):
        """Get flip angle (degree)"""
        return self.FA
    
    
class FourCompartmentModel(Model):
    def __init__(self, TR, FA, T1_array, alpha_array, c_array):
        """This class models a with 4 multicompartment model (Sodium MRI),
        one compartment solid, empty in Sodium
    
    
        Its attributes:
            FA: float
                Flip angle
            TR: float/int
                repetition time of the sequence
            T1_array: float/int - size 3
                relaxation time T1 all compartment
            alpha_array: float/int - size 1 to 3 
                alpha all 3 compartment
            c_array: array float, size 1 to 3, standart value 140 mM
                concentration array all compartments
            __Compartments: dictionnary of compartments
                private
            
        Its public methods:
            get_TSC sum of each compartment function non weight,
            get/set for each attribute public named accordingly
            get_attribute
            set_attribute
            get_CSF signal compartment CSF
            get_ISC signal compartment ISC
            get_ESC signal compartment ESC
        ."""
        
        # Unpacking
        if len(T1_array) != 3 or len(alpha_array)!=3:
            print("T1 array and alpha array expected size 3, real sizes {0}, {1}".format(len(T1_array),
                 len(alpha_array)))
            sys.exit(1)
        if len(c_array) > 3:
            print("Concentration array expected size 3, real size {0}".format(len(c_array)))
            sys.exit(1)
        else:
            if len(c_array) == 3:
                c1 = c_array[0]
                c2 = c_array[1]
                c3 = c_array[2]
            elif len(c_array) == 2:
                c1 = c_array[0]
                c2 = c_array[1]
                c3 = 140 #CSF concentration litterature
            else:
                c1 = c_array[0]
                c2 = 140 #ESC concentration litterature
                c3 = 140 #CSF concentration litterature
                
        # setting
        self.TR = TR
        self.FA = FA
        self.compartments = {}
        
        IC = Compartment(2, alpha_array[0], c1, T1_array[0])
        EC = Compartment(2, alpha_array[1], c2, T1_array[1])
        CSF = Compartment(2, alpha_array[2], c3, T1_array[2])
    
        # adding into the dictionnary
        self.compartments.update({'IC': IC})
        self.compartments.update({'EC': EC}) 
        self.compartments.update({'CSF': CSF})
    
    def get_ISC(self, TE, T2_1, T2_2, f, pond = True):
        IC = self.compartments['IC']
        if pond == True:
            return IC.functionPond(self.FA, self.TR, TE, T2_1, T2_2, f) 
        return IC.function( TE, T2_1, T2_2, f) 
                                
    def get_ESC(self, TE, T2_1, T2_2, f, pond = True):
        EC = self.compartments['EC']
        if pond == True:
            return EC.functionPond(self.FA,self.TR, TE, T2_1, T2_2, f)     
        return EC.function(TE, T2_1, T2_2, f)     
        
    def get_CSF(self, TE, T2_1, T2_2, f, pond = True):
        CSF = self.compartments['CSF']
        if pond == True:
            return CSF.functionPond(self.FA, self.TR, TE, T2_1, T2_2, f)
        return CSF.function( TE, T2_1, T2_2, f)
  
    def set_T1All(self, T1_array):
        """Set all 3 T1 relaxation time [ms]"""
        
        if len(T1_array) != 3 :
            print("T1 array expected size 3, real sizes {0}".format(len(T1_array)))
            sys.exit(1)
        else:
            # Copy all the values from the anterior compartment
            # change only the T1
            IC_new =  copy.deepcopy(self.compartments['IC'])
            IC_new.set_T1(T1_array[0])
            
            EC_new =  copy.deepcopy(self.compartments['EC'])
            EC_new.set_T1(T1_array[1])
        
            CSF_new =  copy.deepcopy(self.compartments['CSF'])
            CSF_new.set_T1(T1_array[2])
            
            # Changing objects in the dictionnary
            self.compartments.update({'IC': IC_new})
            self.compartments.update({'EC': EC_new}) 
            self.compartments.update({'CSF': CSF_new})
    
    def set_alphaAll(self, alpha_array):        
        """Set all 3 alpha volume fraction"""
        
        if len(alpha_array) != 3 :
            print("Alpha array expected size 3, real sizes {0}".format(len(alpha_array)))
            sys.exit(1)
        else:
            # Copy all the values from the anterior compartment
            # change only the alpha
     
            IC_new =  copy.deepcopy(self.compartments['IC'])
            IC_new.set_volFraction(alpha_array[0])
            
            EC_new =  copy.deepcopy(self.compartments['EC'])
            EC_new.set_volFraction(alpha_array[1])
        
            CSF_new =  copy.deepcopy(self.compartments['CSF'])
            CSF_new.set_volFraction(alpha_array[2])
            
            # Changing objects in the dictionnary
            self.compartments.update({'IC': IC_new})
            self.compartments.update({'EC': EC_new}) 
            self.compartments.update({'CSF': CSF_new})
  