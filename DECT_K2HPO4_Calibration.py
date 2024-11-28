"""
    This code was adapted from the University of Calgary Bonelab DECT_PhantomCalibration.py (https://github.com/Bonelab/DECT_BoneAnalysis/tree/master)
    The code has been modified for use with the Mindways solid CT K2HPO4 calibration phantom with five rods labelled A to E.
"""

# Import all modules for the function
import SimpleITK as sitk
import sys
import os
import numpy as np
import pandas as pd
import argparse
from bonelab.util.echo_arguments import echo_arguments


def DECT_Calibration(filePath,lowenergy_filename,highenergy_filename,mask_fnm):
    #read in images:
    mono_low_fnm = lowenergy_filename
    mono_high_fnm = highenergy_filename
    
    mono_low = sitk.ReadImage(filePath+'/'+mono_low_fnm+'.nii', sitk.sitkFloat32)
    mono_high = sitk.ReadImage(filePath+'/'+mono_high_fnm+'.nii', sitk.sitkFloat32)
    
    img_size = mono_low.GetSize()
    img_spacing = mono_low.GetSpacing()
    img_origin = mono_low.GetOrigin()
    img_direction = mono_low.GetDirection()


    #Calculate average intensity in phantom and find line of best fit:
    RodE_label = 5
    RodD_label = 4
    RodC_label = 3
    RodB_label = 2
    RodA_label = 1

    phantom_mask = sitk.ReadImage(filePath+'/'+mask_fnm+'.nii')
    RodE_mask = phantom_mask == RodE_label
    RodD_mask = phantom_mask == RodD_label
    RodC_mask = phantom_mask == RodC_label
    RodB_mask = phantom_mask == RodB_label
    RodA_mask = phantom_mask == RodA_label

    RodE_phantom_low = sitk.Mask(mono_low,RodE_mask)
    RodD_phantom_low = sitk.Mask(mono_low,RodD_mask)
    RodC_phantom_low = sitk.Mask(mono_low,RodC_mask)
    RodB_phantom_low = sitk.Mask(mono_low,RodB_mask)
    RodA_phantom_low = sitk.Mask(mono_low,RodA_mask)

    RodE_phantom_high = sitk.Mask(mono_high,RodE_mask)
    RodD_phantom_high = sitk.Mask(mono_high,RodD_mask)
    RodC_phantom_high = sitk.Mask(mono_high,RodC_mask)
    RodB_phantom_high = sitk.Mask(mono_high,RodB_mask)
    RodA_phantom_high = sitk.Mask(mono_high,RodA_mask)


    RodE_array_low = sitk.GetArrayFromImage(RodE_phantom_low)
    RodE_array_high = sitk.GetArrayFromImage(RodE_phantom_high)
    RodE_mask_array = sitk.GetArrayFromImage(RodE_mask)
    
    RodD_array_low = sitk.GetArrayFromImage(RodD_phantom_low)
    RodD_array_high = sitk.GetArrayFromImage(RodD_phantom_high)
    RodD_mask_array = sitk.GetArrayFromImage(RodD_mask)
    
    RodC_array_low = sitk.GetArrayFromImage(RodC_phantom_low)
    RodC_array_high = sitk.GetArrayFromImage(RodC_phantom_high)
    RodC_mask_array = sitk.GetArrayFromImage(RodC_mask)
    
    RodB_array_low = sitk.GetArrayFromImage(RodB_phantom_low)
    RodB_array_high = sitk.GetArrayFromImage(RodB_phantom_high)
    RodB_mask_array = sitk.GetArrayFromImage(RodB_mask)
    
    RodA_array_low = sitk.GetArrayFromImage(RodA_phantom_low)
    RodA_array_high = sitk.GetArrayFromImage(RodA_phantom_high)
    RodA_mask_array = sitk.GetArrayFromImage(RodA_mask)

    HU_RodE_low = np.average(RodE_array_low, weights=RodE_mask_array)
    print('HU_RodE_low: '+str(HU_RodE_low))
    
    HU_RodD_low = np.average(RodD_array_low, weights=RodD_mask_array)
    print('HU_RodD_low: '+str(HU_RodD_low))
    
    HU_RodC_low = np.average(RodC_array_low, weights=RodC_mask_array)
    print('HU_RodC_low: '+str(HU_RodC_low))
    
    HU_RodB_low = np.average(RodB_array_low, weights=RodB_mask_array)
    print('HU_RodB_low: '+str(HU_RodB_low))
    
    HU_RodA_low = np.average(RodA_array_low, weights=RodA_mask_array)
    print('HU_RodA_low: '+str(HU_RodA_low))

    HU_RodE_high = np.average(RodE_array_high, weights=RodE_mask_array)
    print('HU_RodE_high: '+str(HU_RodE_high))
    
    HU_RodD_high = np.average(RodD_array_high, weights=RodD_mask_array)
    print('HU_RodD_high: '+str(HU_RodD_high))
    
    HU_RodC_high = np.average(RodC_array_high, weights=RodC_mask_array)
    print('HU_RodC_high: '+str(HU_RodC_high))
    
    HU_RodB_high = np.average(RodB_array_high, weights=RodB_mask_array)
    print('HU_RodB_high: '+str(HU_RodB_high))
    
    HU_RodA_high = np.average(RodA_array_high, weights=RodA_mask_array)
    print('HU_RodA_high: '+str(HU_RodA_high))

    y_low = [HU_RodA_low, HU_RodB_low, HU_RodC_low, HU_RodD_low, HU_RodE_low]
    x = [-51.83, -53.40, 58.88, 157.05, 375.83]
    m_low,b_low = np.polyfit(x, y_low, 1)
    print('m_low: '+str(m_low))
    print('b_low: '+str(b_low))

    y_high = [HU_RodA_high, HU_RodB_high, HU_RodC_high, HU_RodD_high, HU_RodE_high]
    x = [-51.83, -53.40, 58.88, 157.05, 375.83]
    m_high,b_high = np.polyfit(x, y_high, 1)
    print('m_high: '+str(m_high))
    print('b_high: '+str(b_high))

    #save the slope and offset values to csv file
    d = {'m_low': [m_low], 'b_low': [b_low], 'm_high': [m_high], 'b_high': [b_high]}
    df = pd.DataFrame(data=d)
    df.to_csv(filePath+'/CalibrationSlope&Offset_DECT.csv',mode='a')

    #Convert image from HU to mgK2HPO4 using best fit line from phantom.
    #Based on method from Sfeir et al., Bone 2018
    mono_le_array = sitk.GetArrayFromImage(mono_low)
    mono_he_array = sitk.GetArrayFromImage(mono_high)
    mgK2HPO4_array = ((mono_he_array - b_high) - (mono_le_array - b_low))/(m_high - m_low)
    img_K2HPO4 = sitk.GetImageFromArray(mgK2HPO4_array)
    img_K2HPO4.SetSpacing(img_spacing)
    img_K2HPO4.SetOrigin(img_origin)
    img_K2HPO4.SetDirection(img_direction)


    sitk.WriteImage(img_K2HPO4,filePath+'/'+'Calibrated_DECT.nii',True)


def main():
    # Set up description
    description='''Function to calibrate bone density images based on dual-energy CT input.

    This program will read in two simulated monoenergetic images (generated by GE), as well as a mask image
    indicating the location of calibration rods (mgK2HPO4). This script is meant to be used on images scanned with
    the Mindways solid CT calibration phantom. Phantom mask image should have the following values:
    375.83 mgK2HPO4 rod label = 5
    157.05 mgK2HPO4 rod label = 4
    58.88 mgK2HPO4 rod label = 3
    -53.40 mgK2HPO4 rod label = 2
    -51.83 mgK2HPO4 rod label = 1

    The output will be an image that's calibrated using a DECT calibration method, such that hounsfield units are converted to mgK2HPO4/ccm based on
    correlation to the Mindways solid CT calibration phantom.

    DECT calibration method based on pubications by Sfeir et al., Bone 2018, https://pubmed.ncbi.nlm.nih.gov/29704696/
    and Gluer et al., J Computer Assisted Tomography 1988, https://pubmed.ncbi.nlm.nih.gov/3351039/

'''


    # Set up argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="DECT_K2HPO4_Calibration",
        description=description
    )

    parser.add_argument("filePath",
                        type=str,
                        help="The filepath")
    parser.add_argument("--lowenergy_filename","-le",
                        type=str,
                        help="Filename for the low energy simulated monochromatic image")
    parser.add_argument("--highenergy_filename","-he",
                        type=str,
                        help="Filename for the high energy simulated monochromatic image")
    parser.add_argument("--mask_fnm","-m",
                        default='K2HPO4_mask',
                        type=str,
                        help="Filename for the K2HPO4 rod mask image")


    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('DECT_K2HPO4_Calibration', vars(args)))

    # Run program
    DECT_Calibration(**vars(args))


if __name__ == '__main__':
    main()