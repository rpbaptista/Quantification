# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 15:37:16 2019

@author: Renata Porciuncula Baptista
@email: renata.porciunculabaptista@cea.fr

This library contains functions that are able to extract the
quantification tube of the image most important function
getTube3D, getTube

"""

import cv2
import numpy as np
import matplotlib.pyplot as plt
import random as rng
import math
rng.seed(12345)

#TO DO; implement OTSU in 32 bits

def getConnectedComponentsInSizeRange(image,  max_size, min_size = 0):
    '''
        General function to extract connected components within a size range
        
        output : binary image, masking objects that are within a given size range
    '''
    
    #find all your connected components (white blobs in your image)
    nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(image, connectivity=8)
    
    #connectedComponentswithStats yields every seperated component with information on each of them, such as size
    #the following part is just taking out the background which is also considered a component, but most of the time we don't want that.
    sizes = stats[1:, -1]; nb_components = nb_components - 1
     
    #your answer image
    img2 = np.zeros((output.shape), dtype = np.uint8)
    
    # counter for the components withing the range
    nb_components_selected = 0
    
    #for every component in the image, you keep it only if it's above min_size
    for i in range(0, nb_components):
        if sizes[i] >= min_size and sizes[i] <= max_size:
            img2[output == i + 1] = 255
            nb_components_selected = nb_components_selected + 1
    
  #  plt.imshow(np.hstack([image, img2]), cmap='gray')
  #  plt.show()    
#        
    return img2, nb_components_selected



def getIndexRoundAsPossible(img, nb_components, number_tubes):
    '''
      Compute perimeter and area and see if the relationship is closer 
      with the circle one
      input : image 2D
      output : N images with the tubes masks
    '''
    if nb_components >= 2:
        # Find contours
        contours, hierarchy = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)[-2:]
    
        # Initialize variables 
        area = np.zeros(len(contours))
        perimeter = np.zeros(len(contours))
        diff = np.zeros(len(contours))
    
        # Compute perimeter and area
        i = 0
        for cnt in contours:
            perimeter[i] = cv2.arcLength(cnt, True)
            area[i] = cv2.contourArea(cnt)
            i = i + 1
        
        # Compute r from two metrics, closer two zero, closer to a perfect circle 
        r_perimeter = perimeter/(2*math.pi)
        r_area = np.sqrt(area/math.pi)
        r_max = np.maximum(r_area,r_perimeter)
        diff = np.abs(r_area - r_perimeter)/r_max

       # Getting idx of the objects more "circle"
        diff_sorted = np.sort(diff)
        dissymmetry_tubes = (diff_sorted[1] - diff_sorted[0])
        print(dissymmetry_tubes)
        threshold = diff_sorted[number_tubes - 1]
        idx = np.where(diff <= threshold)[0]
        masks = np.zeros((number_tubes, img.shape[0],img.shape[1]))

        idx_tube = 0
        for i in range(2):
            j = idx[i]
            drawing = np.zeros((img.shape[0],img.shape[1], 3), dtype=np.uint8)
            color = (255,255,255)
            cv2.drawContours(drawing, [contours[j]], 0, color, -1, cv2.LINE_8, hierarchy, 0)
            drawing_shape = drawing.shape
            if (dissymmetry_tubes < 0.20):
                if len(drawing_shape) == 3:
            #        print(idx_tube)
                    masks[idx_tube,:,:] = drawing[:,:,0]
            idx_tube = idx_tube + 1

        #return drawing[:,:,0]    
        return masks
    else:
        return np.zeros(img.shape)

def getTube(image, number_tubes = 2):

    # binarizing the image
    ret, th = cv2.threshold(image, 0, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    
#    plt.imshow(np.hstack([image, th]), cmap='gray')
#    plt.show()    
#    
    # minimum and maximize size of particles we want to keep (number of pixels)
    #here, it's a fixed value, but you can set it as you want, eg the mean of the sizes or whatever
    MIN_SIZE = 20
    MAX_SIZE = 300
  
    # First level selection : getting components connected within a range of value  
    output_im, nb_components = getConnectedComponentsInSizeRange(th,  MAX_SIZE, MIN_SIZE)

    # Second level selection: simmetry
    # Part 1 : retrivieng the more round N objects, N = number_tubes
    output_im = getIndexRoundAsPossible(output_im, nb_components, number_tubes)
    
    return output_im

def getTube3D(image3D):
    # Getting sizes
    N, M, P = image3D.shape
    number_tubes = 2
    
    # initialize the output image
    output = np.zeros((number_tubes, N, M, P), dtype = np.uint8)
    nb_pixels = np.full((number_tubes, P),np.inf)
    for i in range(P):
        ch1 = image3D[:,:,i]
        print("Z:",i)
        output[:,:,:,i] = getTube(ch1,  number_tubes) 
        for j in range(number_tubes):
            if np.sum(np.sum(output[j,:,:,i])) > 0:
                nb_pixels[j, i] = np.sum(np.sum(output[j,:,:,i]))/np.max(np.max(output[j,:,:,i]))
    
    # This part aim to find the transition between the tubes
    #print(nb_pixels,np.argmin(nb_pixels[0,:]),np.argmin(nb_pixels[1,:]))
    if np.argmin(nb_pixels[0,:])==np.argmin(nb_pixels[1,:]):
        return output, np.argmin(nb_pixels[0,:])
    return output, -1

#------------------- NOT USED
def getIndexSymmetric(img, nb_components, number_tubes):
    '''
        Goal of this script is two return the two most symmetric object
        based on their center
    '''
    
    if nb_components > number_tubes:
        nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(img, connectivity=8)
        
        #connectedComponentswithStats yields every seperated component with information on each of them, such as size
        #the following part is just taking out the background which is also considered a component, but most of the time we don't want that.
        centroids = centroids[1:nb_components,:]
        nb_components = nb_components - 1 
       
        # Compute the differennces
        symmetry_matrix = np.full((nb_components, nb_components,2), np.inf) 
        for i in range(nb_components):
            for j in range(i):
                diff =  centroids[i] - centroids[j]
                symmetry_matrix[i,j,0] = diff[0]
                symmetry_matrix[i,j,1] = diff[1]
                
        # Analyze the differences
        min_x = np.min(symmetry_matrix[:,:,0])
        min_y = np.min(symmetry_matrix[:,:,1])
        
        if min_x < min_y:    
            res = np.unravel_index(np.argmin(symmetry_matrix[:,:,0], axis=None), symmetry_matrix.shape)
            return res, centroids
    return None