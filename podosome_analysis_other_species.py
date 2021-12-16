#!/usr/bin/env python

import importlib
from perform_pod_species_analysis import *

import os

import warnings
warnings.filterwarnings('ignore')

"""
Input parameters Paxillin
"""
image_loc = '../../../SH-CH/paxillin_actin/'
save_folder = 'AP_Analyzed'
save_file_name = 'AP_Analyzed/AP'
other_species = 'Paxillin'
files = ['001','002']
pix_sizes = [26.65/658,24.46/604] #um
pod_analysis_distance = 1. #um

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

perform_pod_species_analysis(image_loc,save_file_name,other_species,files,pix_sizes,pod_analysis_distance)

"""
Input parameters Talin
"""
image_loc = '../../../SH-CH/talin_actin/'
save_folder = 'AT_Analyzed'
save_file_name = 'AT_Analyzed/AT'
other_species = 'Talin'
files = ['002','004','008','010']
pix_sizes = [41.47/1024,41.47/1024,41.47/1024,28.35/700] #um
pod_analysis_distance = 1. #um

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

perform_pod_species_analysis(image_loc,save_file_name,other_species,files,pix_sizes,pod_analysis_distance)


"""
Input parameters A-Actinin
"""
image_loc = '../../../SH-CH/actinin_actin/'
save_folder = 'AA_Analyzed'
save_file_name = 'AA_Analyzed/AA'
other_species = 'Actinin'
files = ['003','010','011']
pix_sizes = [41.47/1024,41.47/1024,41.47/1024] #um
pod_analysis_distance = 1. #um

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

perform_pod_species_analysis(image_loc,save_file_name,other_species,files,pix_sizes,pod_analysis_distance)


"""
Input parameters Myosin
"""
image_loc = '../../../SH-CH/myosin_actin/'
save_folder = 'AM_Analyzed'
save_file_name = 'AM_Analyzed/AM'
other_species = 'Myosin'
files = ['003','010']
pix_sizes = [41.47/1024,41.47/1024] #um
pod_analysis_distance = 1.5 #um

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

perform_pod_species_analysis(image_loc,save_file_name,other_species,files,pix_sizes,pod_analysis_distance)

