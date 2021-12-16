#!/usr/bin/env python

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy as scpy
from scipy import stats
from scipy import ndimage
from scipy.signal import find_peaks, peak_widths
import pickle
import pandas as pd
from skimage import io

import importlib
import helper_functions
importlib.reload(helper_functions)
from helper_functions import *

import warnings
warnings.filterwarnings('ignore')

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.family'] = 'arial'
# rcParams['font.sans-serif'] = ['Tahoma']
rcParams.update({'font.size': 16})

def perform_pod_species_analysis(image_loc,save_file_name,other_species,files,pix_sizes,pod_analysis_distance):

    all_clusters = []
    all_centers = []
    for i,file in enumerate(files):
        images = io.imread(image_loc + file + '_Composite.tif')
        file_name = save_file_name+file

        actin = images[1]
        other = images[0]

        pix_size = pix_sizes[i]
        len_micron = int(1/pix_size)
        pod_filt = int(0.35*len_micron)
        site_filt = int(1.5*len_micron)
        pods, sites = find_pod_and_sites_ph(actin,pod_filt,site_filt,plot_bool=False,plot_pers=False,save_file = file_name)
        clusters,centers,radii = cluster_refine_pods_and_sites(pods,sites,actin,pix_size,plot_bool=False,save_file = file_name,upper_lim=3,lower_lim=1)
        all_clusters.append(clusters)
        all_centers.append(centers)

    # podosome_centers = np.array([list(point) for cluster in clusters for point in cluster])
    # cxs = podosome_centers.T[0]
    # cys = podosome_centers.T[1]

    all_rad_act = []
    all_rad_other = []
    for i,file in enumerate(files):
        images = io.imread(image_loc + file + '_Composite.tif')
        file_name = save_file_name+file

        actin = images[1]
        other = images[0]
        pix_size = pix_sizes[i]
        
        clusters = all_clusters[i]
        centers = all_centers[i]
        podosome_centers = np.array([list(point) for cluster in clusters for point in cluster])
        cxs = podosome_centers.T[0]
        cys = podosome_centers.T[1]

        rad_profs_act, x_ums = radial_averaging_from_pods(cxs,cys,actin,pix_size,pod_analysis_distance)
        rad_profs_other, x_ums = radial_averaging_from_pods(cxs,cys,other,pix_size,pod_analysis_distance)
        [all_rad_act.append(prof) for prof in rad_profs_act]
        [all_rad_other.append(prof) for prof in rad_profs_other]


    plt.figure(figsize=(8,6))
    mean_actin = np.mean(all_rad_act,axis=0)
    ci_actin = stats.sem(all_rad_act,axis=0)*1.96
    plt.plot(x_ums,mean_actin/np.max(mean_actin),lw=2,c='m')

    mean_other = np.mean(all_rad_other,axis=0)
    ci_other = stats.sem(all_rad_other,axis=0)*1.96
    plt.plot(x_ums,mean_other/np.max(mean_other),lw=2,c='g')

    plt.fill_between(x_ums, (mean_actin-ci_actin)/np.max(mean_actin), (mean_actin+ci_actin)/np.max(mean_actin), \
                     alpha=0.3,color='m') 
    plt.fill_between(x_ums, (mean_other-ci_other)/np.max(mean_other), (mean_other+ci_other)/np.max(mean_other), \
                     alpha=0.3,color='g') 


    plt.title('Pooled (All Cells)')
    plt.ylabel('Normalized Intensity (a.u)')
    plt.xlabel('Radial Distance From Center of Podosome (\u03BCm)')
    plt.legend(['Actin',other_species])
    plt.tight_layout()

    plt.savefig(save_file_name+'PooledRadial.pdf',dpi=300)
    plt.close()

    all_through_act = []
    all_through_other = []

    for i,file in enumerate(files):
        images = io.imread(image_loc + file + '_Composite.tif')

        actin = images[1]
        other = images[0]
        pix_size = pix_sizes[i]
        
        clusters = all_clusters[i]
        centers = all_centers[i]
        
        for i,center in enumerate(centers):
            cx,cy = center
            for px,py in clusters[i]:
                thru_act = radial_linescan(actin,cx,cy,px,py,r=4*pod_analysis_distance /pix_size)
                thru_other = radial_linescan(other,cx,cy,px,py,r=4*pod_analysis_distance /pix_size)
            all_through_act.append(thru_act)
            all_through_other.append(thru_other)

    plt.figure(figsize=(8,6))
    xs = np.linspace(0,4*pod_analysis_distance ,len(all_through_act[0]))
    plt.plot(xs,np.mean(all_through_act,axis=0)/np.max(np.mean(all_through_act,axis=0)),'m')
    plt.plot(xs,np.mean(all_through_other,axis=0)/np.max(np.mean(all_through_other,axis=0)),'g')
    plot_fill_between(plt,xs,np.mean(all_through_act,axis=0)/np.max(np.mean(all_through_act,axis=0)), \
                      1.96*stats.sem(all_through_act,axis=0)/np.max(np.mean(all_through_act,axis=0)), \
                      alpha=0.3, color='m')
    plot_fill_between(plt,xs,np.mean(all_through_other,axis=0)/np.max(np.mean(all_through_other,axis=0)), \
                      1.96*stats.sem(all_through_other,axis=0)/np.max(np.mean(all_through_other,axis=0)), \
                      alpha=0.3, color='g')

    plt.ylabel('Normalized Intensity (a.u.)')
    plt.xlabel('Radial Distance From Center of Phagocytosis Site (\u03BCm)')
    plt.title('Pooled (All Cells)')
    plt.legend(['Actin',other_species])
    plt.savefig(save_file_name+'ThroughPooled.pdf',dpi=300)

    all_along_act = []
    all_out_act = []
    all_in_act = []
    all_along_other = []
    all_out_other = []
    all_in_other = []

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i,file in enumerate(files):
        images = io.imread(image_loc + file + '_Composite.tif')
        file_name = save_file_name+file

        actin = images[1]
        other = images[0]
        pix_size = pix_sizes[i]
        
        clusters = all_clusters[i]
        centers = all_centers[i]
        
        along_actin, out_actin, in_actin = along_out_in_scans_from_pods(actin,centers,clusters,pix_size,pod_analysis_distance)
        along_other, out_other, in_other = along_out_in_scans_from_pods(other,centers,clusters,pix_size,pod_analysis_distance)
        [all_along_act.append(prof) for prof in along_actin]
        [all_out_act.append(prof) for prof in out_actin]
        [all_in_act.append(prof) for prof in in_actin]
        
        [all_along_other.append(prof) for prof in along_other]
        [all_out_other.append(prof) for prof in out_other]
        [all_in_other.append(prof) for prof in in_other]
        

    xs = np.linspace(0,pod_analysis_distance,len(along_actin[0]))

    fig, axes = plt.subplots(ncols=2,figsize=(10,4))
    out_in_xs = np.hstack([-xs[:1:-1],xs])
    out_in_act = [np.hstack([all_in_act[i][:1:-1],all_out_act[i]]) for i in range(len(all_in_act))]
    out_in_act_mean = np.mean(out_in_act,axis=0)

    out_in_other = [np.hstack([all_in_other[i][:1:-1],all_out_other[i]]) for i in range(len(all_in_other))]
    out_in_other_mean = np.mean(out_in_other,axis=0)


    axes[0].plot(out_in_xs,out_in_act_mean/np.max(out_in_act_mean),color='m')
    axes[0].plot(out_in_xs,out_in_other_mean/np.max(out_in_other_mean),color='g')
    plot_fill_between(axes[0],out_in_xs,out_in_act_mean/np.max(out_in_act_mean), \
                      1.96*stats.sem(out_in_act,axis=0)/np.max(out_in_act_mean), alpha=0.3, color='m')
    plot_fill_between(axes[0],out_in_xs,out_in_other_mean/np.max(out_in_other_mean), \
                      1.96*stats.sem(out_in_other,axis=0)/np.max(out_in_other_mean), alpha=0.3, color='g')
    axes[0].legend(['Actin',other_species])
    axes[0].set_title('Linescan')
    axes[0].set_xlabel('Distance From Podosome (\u03BCm)')
    axes[0].set_ylabel('Normalized Intensity (a.u.)')
    plt.tight_layout()

    along_act_mean = np.mean(all_along_act,axis=0)
    along_other_mean = np.mean(all_along_other,axis=0)

    axes[1].plot(xs,along_act_mean/np.max(along_act_mean),c='m')
    axes[1].plot(xs,along_other_mean/np.max(along_other_mean),c='g')

    plot_fill_between(axes[1],xs,along_act_mean/np.max(along_act_mean), \
                      1.96*stats.sem(all_along_act,axis=0)/np.max(along_act_mean), alpha=0.3, color='m')
    plot_fill_between(axes[1],xs,along_other_mean/np.max(along_other_mean), \
                      1.96*stats.sem(all_along_other,axis=0)/np.max(along_other_mean), alpha=0.3, color='g')


    axes[1].set_title('Along Arcscan')
    axes[1].set_xlabel('Arc Distance From Podosome (\u03BCm)')
    axes[1].set_ylabel('Normalized Intensity (a.u.)')
    plt.tight_layout()

    plt.savefig(save_file_name+'PooledAOI.pdf',dpi=300)

    plt.close()

    # Get individual quantifications for other species peak location, FWHM, and diameter
    all_dist2peak_in = []
    all_FWHM_in= []
    all_dist2peak_out = []
    all_FWHM_out= []
    all_diameter = []
    xs = np.linspace(0,pod_analysis_distance,len(all_in_other[0]))
    dx = xs[1]-xs[0]

    for i in range(len(all_in_other)):
        curr_in = all_in_other[i]
        curr_out = all_out_other[i]
        peaks_in, _ = find_peaks(curr_in, prominence=0.25*np.max(curr_in), height = np.max(curr_in))
        if len(peaks_in) > 0:
            results_half_in = peak_widths(curr_in/np.max(curr_in), peaks_in, rel_height=0.5)
            all_dist2peak_in.append(xs[peaks_in][0])
            all_FWHM_in.append(dx*(results_half_in[3]-results_half_in[2])[0])
            
        peaks_out, _ = find_peaks(curr_out, prominence=0.25*np.max(curr_out), height = np.max(curr_out))
        if len(peaks_out) > 0:
            results_half_out = peak_widths(curr_out/np.max(curr_out), peaks_out, rel_height=0.5)
            all_dist2peak_out.append(xs[peaks_out][0])
            all_FWHM_out.append(dx*(results_half_out[3]-results_half_out[2])[0])
        if len(peaks_in) > 0 and len(peaks_out) > 0:
            all_diameter.append(np.abs(xs[peaks_out][0]+xs[peaks_in][0]))
        

            
    all_dist2peak_in = np.array(all_dist2peak_in)
    all_FWHM_in = np.array(all_FWHM_in)
    all_dist2peak_out = np.array(all_dist2peak_out)
    all_FWHM_out = np.array(all_FWHM_out)
    all_diameter = np.array(all_diameter)

    distances = np.hstack([all_dist2peak_in,all_dist2peak_out])
    FWHM = np.hstack([all_FWHM_in,all_FWHM_out])
    directions = np.hstack([ \
        ['In Peak']* len(all_dist2peak_in) ,\
        ['Out Peak']* len(all_dist2peak_out) ,\
    ])

    other_df = pd.DataFrame(list(zip(distances,FWHM, directions)),columns = ['Distance','FWHM','Direction'])

    t, p = scpy.stats.ttest_ind(all_dist2peak_in,all_dist2peak_out)
    plt.figure(figsize=(4,6))
    sns.boxplot(x='Direction', y='Distance', data=other_df,showfliers=False,palette="Set2")
    sns.swarmplot(x='Direction', y='Distance', data=other_df,color=".2")

    plt.ylabel('Distance (\u03BCm)')
    plt.xlabel('')
    plt.title(other_species + ' Distance to Peak')

    max_dp = np.max(other_df['Distance'])
    plt.hlines(1.05*max_dp,0,1)
    plt.vlines(0,1.025*max_dp,1.05*max_dp)
    plt.vlines(1,1.025*max_dp,1.05*max_dp)
    plt.text(0.5,1.05*max_dp,'P = %.2g' %p, ha='center',va='bottom')
    plt.tight_layout()

    plt.savefig(save_file_name+'IndPeakDistances.pdf',dpi=300)
    plt.close()

    plt.figure(figsize=(3,6))

    sns.boxplot(y=all_diameter,showfliers=False,color='lightblue')
    sns.swarmplot(y=all_diameter,color=".2")

    plt.ylabel('Diameter (\u03BCm)')
    plt.xlabel('')
    plt.title(other_species)

    plt.tight_layout()
    plt.savefig(save_file_name+'IndRingDiameter.pdf',dpi=300)
    plt.close()

    t, p = scpy.stats.ttest_ind(all_FWHM_in,all_FWHM_out)

    plt.figure(figsize=(4,6))
    sns.boxplot(x='Direction', y='FWHM', data=other_df,showfliers=False,palette="Set2")
    sns.swarmplot(x='Direction', y='FWHM', data=other_df,color=".2")
    plt.title(other_species + ' FWHM')
    plt.ylabel('FWHM (\u03BCm)')
    plt.xlabel('')

    max_dp = np.max(other_df['FWHM'])
    plt.hlines(1.05*max_dp,0,1)
    plt.vlines(0,1.025*max_dp,1.05*max_dp)
    plt.vlines(1,1.025*max_dp,1.05*max_dp)
    plt.text(0.5,1.05*max_dp,'P = %.2g' %p, ha='center',va='bottom')
    plt.tight_layout()
    plt.savefig(save_file_name+'IndFWHM.pdf',dpi=300)
    plt.close()

    quant_dict = {'all_dist2peak_in':all_dist2peak_in,'all_dist2peak_out':all_dist2peak_out, \
             'all_FWHM_in':all_FWHM_in, 'all_FWHM_out':all_FWHM_out,'all_diameter':all_diameter, \
             'num_in_peaks': len(all_dist2peak_in), 'num_out_peaks': len(all_dist2peak_out), \
             'num_diameter': len(all_diameter), 'num_cells': len(files), \
             'num_pod':len(np.vstack(np.hstack(all_clusters))),'pod':all_clusters, \
             'num_sites':len(np.vstack(np.vstack(all_centers))),'sites':all_centers,
             'rad_act': all_rad_act, 'rad_other':all_rad_other,
             'along_act':all_along_act, 'along_other':all_along_other,
             'in_act':all_in_act, 'in_other':all_in_other,
             'out_act':all_out_act, 'out_other':all_out_other}

    pickle.dump(quant_dict,open(save_file_name + 'IndividualQuant.pickle','wb'))