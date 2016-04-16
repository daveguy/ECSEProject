
# coding: utf-8

# In[ ]:

# import os
import numpy as np
import SimpleITK
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from __future__ import division
from datetime import datetime as dt
get_ipython().magic(u'pylab inline')


# A method to display slices

# In[ ]:

def sitk_show(img, title=None, margin=0.0, dpi=40):
#     nda = SimpleITK.GetArrayFromImage(img)
    nda = img
    #spacing = img.GetSpacing()
    figsize = (1 + margin) * nda.shape[0] / dpi, (1 + margin) * nda.shape[1] / dpi
    #extent = (0, nda.shape[1]*spacing[1], nda.shape[0]*spacing[0], 0)
    extent = (0, nda.shape[1], nda.shape[0], 0)
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([margin, margin, 1 - 2*margin, 1 - 2*margin])
    plt.set_cmap("gray")
    ax.imshow(nda,extent=extent,interpolation=None)
    
    if title:
        plt.title(title)
    
    plt.show()


# Adds as outlier voxels 3 standard deviations away from the norm

# In[ ]:

def add_outlier(seg, img, k):
    for i in range(1,k+1):
        stdev = img[seg==i].std()
        mean = img[seg==i].mean()
        sub = np.all([seg==i, np.any([img>mean+3*stdev, img<mean-3*stdev], 0)],0)
        seg[sub] = 4
        
        


# Takes an image segmentation (gray matter, white matter and CSF) and initializes theta. outliers are more than 3 standard deviations away from norm

# In[ ]:

def init_theta(seg,img, theta, k, c):
    for i in range(0,c):
        for j in range(1,k+2):
            theta[i][2*j] = img[i][seg==(j)].mean()
            theta[i][2*j+1] = img[i][seg==(j)].var()
            
            
    
        


# Calculates the probability of a configuration of t given alpha. Equation 2 in the Menze paper            

# In[ ]:

def T_prob(Ti, alphai):
    return alphai**np.sum(Ti) * (1-alphai)**(Ti.size-np.sum(Ti))


# Calculates the probability of a configuration of y given a configuration of t, a class k and theta.
# Equation 3 in the menze paper

# In[ ]:

def Y_prob(yi, ti, ki, theta):
    mu=[theta[0][2*ki],theta[1][2*ki],theta[2][2*ki],theta[3][2*ki]]
    var=[theta[0][2*ki+1],theta[1][2*ki+1],theta[2][2*ki+1],theta[3][2*ki+1]]
    mukP1 = [theta[0][-2],theta[1][-2],theta[2][-2],theta[3][-2]]
    varkP1 = [theta[0][-1],theta[1][-1],theta[2][-1],theta[3][-1]]
    return np.prod(mlab.normpdf(yi,mu,np.sqrt(var))**(1-ti) * mlab.normpdf(yi,mukP1, np.sqrt(varkP1))**ti)


# returns the most likely class from the segmentation (1,2 or 3) for a given probability map by finding the segmentation index's with the highest corresponding mean probability. Used to match the maps for CSF, WM and GM to a corresponding segmentation

# In[ ]:

def get_class(seg, probMap):
    return np.argmax([np.mean(probMap[seg==1]), np.mean(probMap[seg==2]), np.mean(probMap[seg==3])])+1


# The function qi(ti) and wik(tik). Equation 4 in my report, Equations 5 and later unnumbered equations in Menze paper.

# In[ ]:

def Qi(yi,ti, k, theta, alphai, atlas, index):
    sum = 0
    for i in range(1,k+1): #3 classes are 1,2,3
        sum += Y_prob(yi, ti, i, theta)*T_prob(ti, alphai)*atlas[i][index]
    return sum


# In[ ]:

def Wik(yi, theta, atlas, ti, index, k):
    mu=[theta[0][2*k],theta[1][2*k],theta[2][2*k],theta[3][2*k]]
    var=[theta[0][2*k+1],theta[1][2*k+1],theta[2][2*k+1],theta[3][2*k+1]]
    return atlas[k][index]*np.prod(mlab.normpdf( yi, mu, np.sqrt(var))**(1-ti))


# Returns all possible 2^C configurations of length C for t_i

# In[ ]:

def get_all_t(c): # returnes the 2^c configurations
    if c == 1:
        return [np.array([0]), np.array([1])]
    else:
        arr = get_all_t(c-1)
        return np.array([np.append(a0,b0) for a0 in arr for b0 in np.array([0,1])])


# Performs the alpha update for a single alpha_i. Equation 5 in my report

# In[ ]:

def alpha_update(index):
    return np.sum(map(Qi, [img[:,index]]*all_t.shape[0], all_t, [K]*all_t.shape[0], [theta]*all_t.shape[0], 
        [alpha[index]]*all_t.shape[0], [atlas]*all_t.shape[0], [index]*all_t.shape[0]) * np.true_divide(np.sum(all_t,1), C))


# Functions to calculate the numerators and denominators for the theta update. Equation 6 in my report, unnumbered equations under equation 5 in Menze paper.

# In[ ]:

def mu_update_numerator(index):
    return np.sum(map(Qi, [img[:,index]]*all_t.shape[0], all_t, [K]*all_t.shape[0], [theta]*all_t.shape[0], 
        [alpha[index]]*all_t.shape[0], [atlas]*all_t.shape[0], [index]*all_t.shape[0])*
           np.array(map(Wik, [img[:,index]]*all_t.shape[0],[theta]*all_t.shape[0], [atlas]*all_t.shape[0],
                all_t, [index]*all_t.shape[0], [ki]*all_t.shape[0]))*(1-all_t[:,ci])*img[ci,index])


def update_denominator(index): # mu and stdev update denominators are the same
    return np.sum(map(Qi, [img[:,index]]*all_t.shape[0], all_t, [K]*all_t.shape[0], [theta]*all_t.shape[0], 
        [alpha[index]]*all_t.shape[0], [atlas]*all_t.shape[0], [index]*all_t.shape[0])*
           np.array(map(Wik, [img[:,index]]*all_t.shape[0],[theta]*all_t.shape[0], [atlas]*all_t.shape[0],
                all_t, [index]*all_t.shape[0], [ki]*all_t.shape[0]))*(1-all_t[:,ci]))

def stddev_update_numerator(index):
    return np.sum(map(Qi, [img[:,index]]*all_t.shape[0], all_t, [K]*all_t.shape[0], [theta]*all_t.shape[0], 
        [alpha[index]]*all_t.shape[0], [atlas]*all_t.shape[0], [index]*all_t.shape[0])*
           np.array(map(Wik, [img[:,index]]*all_t.shape[0],[theta]*all_t.shape[0], [atlas]*all_t.shape[0],
                all_t, [index]*all_t.shape[0], [ki]*all_t.shape[0]))*(1-all_t[:,ci])*(img[ci,index]-theta[ci][2*ki])**2)

def tumor_mu_update_numerator(index):
    return np.sum(map(Qi, [img[:,index]]*all_t.shape[0], all_t, [K]*all_t.shape[0], [theta]*all_t.shape[0], 
        [alpha[index]]*all_t.shape[0], [atlas]*all_t.shape[0], [index]*all_t.shape[0])*all_t[:,ci]*img[ci,index])

def tumor_update_denominator(index):
    return np.sum(map(Qi, [img[:,index]]*all_t.shape[0], all_t, [K]*all_t.shape[0], [theta]*all_t.shape[0], 
        [alpha[index]]*all_t.shape[0], [atlas]*all_t.shape[0], [index]*all_t.shape[0])*all_t[:,ci])

def tumor_stddev_update_numerator(index):
    return np.sum(map(Qi, [img[:,index]]*all_t.shape[0], all_t, [K]*all_t.shape[0], [theta]*all_t.shape[0], 
        [alpha[index]]*all_t.shape[0], [atlas]*all_t.shape[0], [index]*all_t.shape[0])*all_t[:,ci]*(img[ci,index]-theta[ci][2*ki])**2)


# Once optimaization is complete, this gets the final labels for tumors

# In[ ]:

def get_T_for_C(index):
    np.sum(map(Qi, [img[:,index]]*all_t.shape[0], all_t, [K]*all_t.shape[0], [theta]*all_t.shape[0], 
        [alpha[index]]*all_t.shape[0], [atlas]*all_t.shape[0], [index]*all_t.shape[0])*all_t[:,ci])


# Initialization of variables, and some loading that had to be done in order

# In[ ]:

print "Loading data and initializing variables"
dir_path = 'BRATS-2/Image_Data/HG/0001/'
brainT1_path = dir_path + 'VSD.Brain.XX.O.MR_T1/VSD.Brain.XX.O.MR_T1.685.mha'
brainT1c_path =dir_path + 'VSD.Brain.XX.O.MR_T1c/VSD.Brain.XX.O.MR_T1c.686.mha'
brainT2_path = dir_path + 'VSD.Brain.XX.O.MR_T2/VSD.Brain.XX.O.MR_T2.687.mha'
brainFLAIR_path = dir_path + 'VSD.Brain.XX.O.MR_Flair/VSD.Brain.XX.O.MR_Flair.684.mha'
brainT1_img = SimpleITK.ReadImage(brainT1_path)
brainT1c_img = SimpleITK.ReadImage(brainT1c_path)
brainT2_img = SimpleITK.ReadImage(brainT2_path)
brainFlair_img = SimpleITK.ReadImage(brainFLAIR_path)
initial_seg_path = 'BRATS-2/Image_Data/HG/0001/seg.nii.gz'
initial_seg = SimpleITK.ReadImage(initial_seg_path)
atlas_path = 'BRATS-2/Image_Data/HG/0001/Atlas/'


# Slices image, slice was manually selected to contain a large amount of tumor

# In[ ]:

brainT1_img = SimpleITK.GetArrayFromImage(brainT1_img)[65,:,:]
brainT1c_img = SimpleITK.GetArrayFromImage(brainT1c_img)[65,:,:]
brainT2_img = SimpleITK.GetArrayFromImage(brainT2_img)[65,:,:]
brainFlair_img =SimpleITK.GetArrayFromImage(brainFlair_img)[65,:,:]
initial_seg = SimpleITK.GetArrayFromImage(initial_seg)[65,:,:]


# In[ ]:

K = 3 #3 tissue classes
C = 4 # my data has 4 channels
N=brainT1_img.size
theta = np.empty([C, 2*(K+2)], double) # first 2 vals in a row are empty to match index's. class one is 2*1 and 2*1+1

img = np.empty([C,N])
img[0] = brainT1_img.flatten()
img[1] = brainT1c_img.flatten()
img[2] = brainT2_img.flatten()
img[3] = brainFlair_img.flatten()
initial_seg = initial_seg.flatten()

atlas = np.empty((4,N), double) # atlas[0] is empty to make index's match with segmentation.
csf = SimpleITK.GetArrayFromImage(SimpleITK.ReadImage(atlas_path + "CSF_warped.nii.gz")[:,:,65]).flatten()
gm = SimpleITK.GetArrayFromImage(SimpleITK.ReadImage(atlas_path + "GM_warped.nii.gz")[:,:,65]).flatten()
wm = SimpleITK.GetArrayFromImage(SimpleITK.ReadImage(atlas_path + "WM_warped.nii.gz")[:,:,65]).flatten()
atlas[get_class(initial_seg, csf)] = csf
atlas[get_class(initial_seg, gm)] = gm
atlas[get_class(initial_seg, wm)] = wm

add_outlier(initial_seg, img[0], K)
alpha = np.full(N, 0.3)
alpha[initial_seg==4] = 0.7

init_theta(initial_seg, img, theta, K, C)
all_t = get_all_t(C)


# #### EM, i.e. main segmentation

# Considering only non-zero( in all channels) voxels

# In[ ]:

indexs = np.where(~np.any(img,axis=0))[0] # index's of non-zero intensity vectors yi
for j in range(0,6): #for now do 10 iterations of EM
    print dt.now().time().isoformat() + " EM iteration " + str(j)
    print "\tPerforming alpha step"
    alpha[indexs] = map(alpha_update, indexs) # alpha update
    for ci in range(0,C):
        print"\t" + dt.now().time().isoformat() + " Performing theta update for normal tissues"
        for ki in range (1,K+1): #update mean and var for normal tissues
            theta[ci][2*ki] = np.true_divide(np.sum(map(mu_update_numerator, indexs)), np.sum(map(update_denominator, indexs)))
            theta[ci][2*ki+1] = np.true_divide(np.sum(map(stddev_update_numerator, indexs)), np.sum(map(update_denominator, indexs)))
        ki = 4  #update mean and var for tumor tissues
        print "\t" + dt.now().time().isoformat() + " Performing theta update for tumor tissue"
        theta[ci][2*ki] = np.true_divide(np.sum(map(tumor_mu_update_numerator, indexs)), np.sum(map(tumor_update_denominator, indexs)))
        theta[ci][2*ki+1] = np.true_divide(np.sum(map(tumor_stddev_update_numerator, indexs)), np.sum(map(tumor_update_denominator, indexs)))


# Get final segmentations

# In[ ]:

print "Computing final segmentations"
indexs = np.where(~np.any(img,axis=0))[0]
final_tumor = np.zeros((C,N), double)
for c in range(0,C):
    final_tumor[c][indexs] = map(get_T_for_C, indexs) 


# Save the segmentations

# In[ ]:

print "Saving images"
t1_seg = SimpleITK.GetImageFromArray(np.reshape(final_tumor[0],(58,72,53)))
t1c_seg = SimpleITK.GetImageFromArray(np.reshape(final_tumor[1],(58,72,53)))
t2_seg = SimpleITK.GetImageFromArray(np.reshape(final_tumor[2],(58,72,53)))
flair_seg = SimpleITK.GetImageFromArray(np.reshape(final_tumor[3],(58,72,53)))
SimpleITK.WriteImage(t1_seg, dir_path + 'my_T1_seg.nii.gz')
SimpleITK.WriteImage(t1c_seg, dir_path + 'my_T1c_seg.nii.gz')
SimpleITK.WriteImage(t2_seg, dir_path + 'my_T2_seg.nii.gz')
SimpleITK.WriteImage(flair_seg, dir_path + 'my_Flair_seg.nii.gz')


# In[ ]:

SimpleITK.WriteImage(SimpleITK.Expand(t1_seg,[3,3,3], SimpleITK.sitkWelchWindowedSinc), dir_path + 'my_T1_seg_exp.nii.gz')
SimpleITK.WriteImage(SimpleITK.Expand(t1c_seg,[3,3,3], SimpleITK.sitkWelchWindowedSinc), dir_path + 'my_T1c_seg_exp.nii.gz')
SimpleITK.WriteImage(SimpleITK.Expand(t2_seg,[3,3,3], SimpleITK.sitkWelchWindowedSinc), dir_path + 'my_T2_seg_exp.nii.gz')
SimpleITK.WriteImage(SimpleITK.Expand(flair_seg,[3,3,3], SimpleITK.sitkWelchWindowedSinc), dir_path + 'my_Flair_seg_exp.nii.gz')

