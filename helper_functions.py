import matplotlib.pyplot as plt
import numpy as np
import scipy as scpy
from scipy import optimize
from scipy import stats
from scipy import ndimage
import csv
import pandas as pd
import dionysus as d
from ripser import ripser, lower_star_img
import persim
import pickle

from sklearn.cluster import KMeans

from mpl_toolkits.mplot3d import Axes3D


#Helps plot confidence intervals
def plot_fill_between(ph,xs, mean, ci, alpha=0.3, color='b'):
    ph.fill_between(xs, (mean-ci),(mean+ci), \
                 alpha=alpha,color=color) 
def plot_fill_between_norm(ph,xs, mean, ci, alpha=0.3, color='b'):
    ph.fill_between(xs, (mean-ci)/np.max(mean),\
                 (mean+ci)/np.max(mean), \
                 alpha=alpha,color=color) 

'''
Slightly adapted from Scipy Cookbook: https://scipy-cookbook.readthedocs.io/items/Least_Squares_Circle.html
Given x and y coordinates, returns the least squares center and radius of a best fitting circle
'''
def find_center(x,y):
    def calc_R(xc, yc):
        """ calculate the distance of each data points from the center (xc, yc) """
        return np.sqrt((x-xc)**2 + (y-yc)**2)

    def f_2b(c):
        """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R(*c)
        return Ri - Ri.mean()

    def Df_2b(c):
        """ Jacobian of f_2b
        The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
        xc, yc     = c
        df2b_dc    = scpy.empty((len(c), x.size))

        Ri = calc_R(xc, yc)
        df2b_dc[0] = (xc - x)/Ri                   # dR/dxc
        df2b_dc[1] = (yc - y)/Ri                   # dR/dyc
        df2b_dc    = df2b_dc - df2b_dc.mean(axis=1)[:, np.newaxis]

        return df2b_dc

    center_estimate = [np.mean(x),np.mean(y)]
    center_2b, ier = optimize.leastsq(f_2b, center_estimate, Dfun=Df_2b, col_deriv=True)

    xc_2b, yc_2b = center_2b
    Ri_2b        = calc_R(*center_2b)
    R_2b         = Ri_2b.mean()
    residu_2b    = sum((Ri_2b - R_2b)**2)
    
    return(center_2b,R_2b)

def dist_between_points(x0,x1,y0,y1):
    return(np.sqrt((x0-x1)**2+(y0-y1)**2))

# Determine the x,y coordinates at distance d, angle theta (rad) from x0, y0
def point_pos (x0, y0, d, theta_rad):
    return(x0 + d*np.cos(theta_rad), y0 + d*np.sin(theta_rad))

# Determine the angle between two vectors    
def vector_angle(a,b,c):
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    
    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    if a[1] < b[1]:
        angle = 2*np.pi - angle
        
    return(angle)

#scale between 0 and 1
def max_min_scaled(dist):
    return((dist-np.min(dist))/(np.max(dist)-np.min(dist)))

# #below used to help stack distribtuions that both have the same data at 0
def out_in_stack(in_dist, out_dist):
    return([np.hstack([in_dist[i][:1:-1],out_dist[i]]) for i in range(len(in_dist))])
def out_in_1d(in_dist, out_dist):
    return(np.hstack([in_dist[::-1],out_dist[1:]]))


'''
Retuns a radial averaged profile given data/image and a center point
From: https://stackoverflow.com/a/21242776
'''
def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = dist_between_points(x,center[0],y,center[1])
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return(radialprofile)

def radial_averaging_from_pods(cxs,cys,image,pix_size_um,max_r_dist = 1):
    ind_max = int(max_r_dist//pix_size_um)+1
    x_ums = np.arange(0,ind_max)*pix_size_um

    rad_profs = []
    for i in range(len(cxs)):
        rad_prof = radial_profile(image, [cxs[i],cys[i]])
        rad_profs.append(rad_prof[:ind_max])

    return(np.array(rad_profs), x_ums)

def out_in_scans_from_pods(image,centers,clusters,pix_size_um,dist_um):

    out_scans = []
    in_scans = []

    for i,center in enumerate(centers):
        cx,cy = center
        
        for px,py in clusters[i]:
            ls_out,ls_in = radial_linescans_from_pod(image,cx,cy,px,py,pix_size_um,dist_um)
            out_scans.append(ls_out)
            in_scans.append(ls_in)

    return(np.array(out_scans),np.array(in_scans))

'''
Gets radial linescans in/outward from a podosome given a center of a phago. site
'''
def radial_linescans_from_pod(z,cx,cy,px,py,pix_size_um,r_dist_um):

    r_dist = r_dist_um/pix_size_um
    r0 = dist_between_points(cx,px,cy,py)
    theta = vector_angle((px,py),(cx,cy),(cx+5,cy))
    
    rs=np.linspace(r0,r0+r_dist,100)
    xs,ys = np.array([point_pos(cx,cy,r,theta) for r in rs]).T # These are in _pixel_ coordinates!!
    # Extract the values along the line, using cubic interpolation
    ls_out = scpy.ndimage.map_coordinates(z, np.vstack((ys,xs)))

    rs=np.linspace(r0,r0-r_dist,100)
    xs,ys = np.array([point_pos(cx,cy,r,theta) for r in rs]).T # These are in _pixel_ coordinates!!
    # Extract the values along the line, using cubic interpolation
    ls_in = scpy.ndimage.map_coordinates(z, np.vstack((ys,xs)))

    return(np.array(ls_out),np.array(ls_in))


'''
Dionysus2 helper functions
'''

'''
Find podosomes and site during frustrated phagocytosis using persistent homology
'''
def find_pod_and_sites_ph(image,pod_filter,site_filter,plot_bool = True, plot_pers=False, save_file = ''):
    img_pod = ndimage.uniform_filter(image.astype(np.float64), size=pod_filter) #Smooth image
    img_pod += 0.01 * np.random.randn(*img_pod.shape) #Make each pixel's value unique
    img_sites = ndimage.uniform_filter(image.astype(np.float64), size=site_filter) #Smooth image
    img_sites += 0.01 * np.random.randn(*img_sites.shape) #Make each pixel's value unique
    
    f_upper_star_pod = d.fill_freudenthal(img_pod, reverse = True) #perform upper star filtration using fill_freudenthal
    p_pod = d.homology_persistence(f_upper_star_pod) #perform PH
    dgms_pod = d.init_diagrams(p_pod,f_upper_star_pod) #get persistence diagrams
    
    f_upper_star_sites = d.fill_freudenthal(img_sites, reverse = True) #perform upper star filtration using fill_freudenthal
    p_sites = d.homology_persistence(f_upper_star_sites) #perform PH
    dgms_sites = d.init_diagrams(p_sites,f_upper_star_sites) #get persistence diagrams
    
    pers_array_pod = dio_pers_array(dgms_pod) #get birth and death times
    pers_array_sites = dio_pers_array(dgms_sites) #get birth and death times
    
    thresh_val_pod = kmeans_pers_thresh(pers_array_pod,0,n=2,plot_bool=plot_pers,save_file=save_file) #h0 homology
    thresh_val_sites = kmeans_pers_thresh(pers_array_sites,1,n=2,plot_bool=plot_pers,save_file=save_file) #h1 homology
    
    idxs = np.arange(pers_array_pod.shape[0])
    idxs0 = idxs[np.logical_and(pers_array_pod[:,0] == 0 , np.abs(pers_array_pod[:, 2] - pers_array_pod[:, 1]) >= thresh_val_pod)]
    idxs = np.arange(pers_array_sites.shape[0])
    idxs1 = idxs[np.logical_and(pers_array_sites[:,0] == 1 , np.abs(pers_array_sites[:, 2] - pers_array_sites[:, 1]) >= thresh_val_sites)]
    
    X, Y = np.meshgrid(np.arange(img_pod.shape[1]), np.arange(img_pod.shape[0]))
    X = X.flatten()
    Y = Y.flatten()
    pods = []
    sites = []

    for idx in idxs0:
        bidx = np.argmin(np.abs(img_pod - pers_array_pod[idx, 1]))
        pods.append([X[bidx],Y[bidx]])

    for idx in idxs1:
        bidx = np.argmin(np.abs(img_sites - pers_array_sites[idx, 2]))
        sites.append([X[bidx],Y[bidx]])
        
    pods = np.asarray(pods)
    sites = np.asarray(sites)
    
    if plot_bool or save_file != '':
        plt.figure(figsize = (8,8))
        plt.imshow(image,cmap = plt.cm.bone)
        plt.axis('off')
        plt.title('Podosomes and Sites Identified')
        for pod in pods:
            plt.scatter(pod[0],pod[1],c='r',marker = 'x')
            
        for site in sites:
            plt.scatter(site[0],site[1],c='y')

        if save_file != '':
            plt.savefig(save_file+"PreRef.pdf",dpi=300)
        if plot_bool: 
            plt.show()
        else: plt.close()
    
    return(pods,sites)


#12.03.2020 swapped the order so that it merges close sites/clusters and then checks for >3 pod. per site
def cluster_refine_pods_and_sites(pods,sites,image,pix_size,plot_bool = True,save_file = '',upper_lim=3,lower_lim=1):
    #default remove podosomes more than 3 micron from estimated site center and then
    #removes podosomes less than a micron from new estimated site center

    nearest_h1 =np.array([pix_size*np.min([np.linalg.norm(pod-site) for site in sites]) for pod in pods])
    pods =  pods[nearest_h1<upper_lim*1.5]
    
    
    #k means cluster based on the number of sites
    km_clusters = KMeans(n_clusters=len(sites))
    km_clusters.fit(pods)

    # #testing initial cluster assignment
    # plt.figure()
    # plt.imshow(image,cmap=plt.cm.bone)
    # plt.scatter(pods.T[0],pods.T[1],c=km_clusters.labels_,cmap=plt.cm.rainbow)

    init_centers = []
    for i in range(np.max(km_clusters.labels_)+1):
        curr_cluster = pods[km_clusters.labels_==i]
        #find center based on a least squares fit to a circle
        try:
            center,_ = find_center(curr_cluster.T[0],curr_cluster.T[1])
        except: center =curr_cluster.T[0][0],curr_cluster.T[1][0]
        init_centers.append(center)

    from sklearn.metrics.pairwise import euclidean_distances
    cent_dist = euclidean_distances(init_centers)
    #make sure that centers aren't found to be too close (2 um) from another (causing split clusters)
    too_close = np.logical_and(cent_dist <2/pix_size , cent_dist != 0)
    # num_too_close = np.sum(np.unique(np.sum(too_close,axis=0)))

    while np.sum(too_close) > 0:
        km_clusters = KMeans(n_clusters=len(init_centers)-1)
        km_clusters.fit(pods)

        init_centers = []
        for i in range(np.max(km_clusters.labels_)+1):
            curr_cluster = pods[km_clusters.labels_==i]
            #find center based on a least squares fit to a circle
            center,_ = find_center(curr_cluster.T[0],curr_cluster.T[1])
            init_centers.append(center)

        cent_dist = euclidean_distances(init_centers)
        too_close = np.logical_and(cent_dist <1/pix_size , cent_dist != 0)

    #make sure that a cluster of podosomes has at least 3 points
    cluster_i, counts_per = np.unique(km_clusters.labels_, return_counts=True)
    num_clus_too_small = np.sum(counts_per < 3)

    if num_clus_too_small > 0:
        km_clusters = KMeans(n_clusters=len(init_centers)-num_clus_too_small)
        km_clusters.fit(pods)


    clusters = []
    centers = []
    radii = []
    dists = []
    for i in range(np.max(km_clusters.labels_)+1):
        curr_cluster = pods[km_clusters.labels_==i]
        #find center based on a least squares fit to a circle
        center,radius = find_center(curr_cluster.T[0],curr_cluster.T[1])

        dist = pix_size*np.asarray([np.linalg.norm(pod-center) for pod in curr_cluster])

        #remove points that are over 3 um from center
        remove_pods_bool = dist>upper_lim

        if np.sum(remove_pods_bool >= 1):
            curr_cluster = curr_cluster[~remove_pods_bool]
            if len(curr_cluster) < 3:
                continue
            center,radius = find_center(curr_cluster.T[0],curr_cluster.T[1])
            dist = pix_size*np.asarray([np.linalg.norm(pod-center) for pod in curr_cluster])

        #remove points that are within 1 um from center 
        remove_pods_bool = dist<lower_lim

        if np.sum(remove_pods_bool >= 1):
            curr_cluster = curr_cluster[~remove_pods_bool]
            if len(curr_cluster) < 3:
                continue
            center,radius = find_center(curr_cluster.T[0],curr_cluster.T[1])
            dist = pix_size*np.asarray([np.linalg.norm(pod-center) for pod in curr_cluster])


        clusters.append(curr_cluster)
        centers.append(center)
        radii.append(radius)
        dists.append(dist)

    clusters = np.array(clusters)
    centers = np.array(centers)
    radii = np.array(radii)

    if plot_bool or save_file != '':
        plt.figure(figsize = (8,8))
        plt.imshow(image,cmap=plt.cm.bone)
        plt.axis('off')
        plt.title('Podosomes and Sites Post-Refinement')
        plt.scatter(centers.T[0],centers.T[1],c='y')
        for clus in clusters:
            plt.scatter(clus.T[0],clus.T[1],marker = 'x',c='r')

        if save_file != '':
            plt.savefig(save_file+"PostRef.pdf",dpi=300)
        if plot_bool: 
            plt.show()
        else: plt.close()

    return(clusters,centers,radii)


#given a dionysus2 dgm obect, return an array with the homology, birth, and death times
def dio_pers_array(dgms_dio):
    pers_array = []
    for h, dgm in enumerate(dgms_dio):
        for pt in dgm:
            pers_array.append([h,pt.birth,pt.death])

    return(np.asarray(pers_array))

#Perform k-means clustering of persistence lifetimes to determine a threshold value
def kmeans_pers_thresh(pers_array,h, n=2, plot_bool=True,save_file=''):
#     from sklearn.cluster import KMeans 

    #get lifetimes
    x = np.abs(np.array(pers_array[:,2]-pers_array[:,1]))
    x = x[pers_array[:,0]==h]
    x = x[np.isfinite(x)]
    x = x.reshape(-1, 1)
    
    km = KMeans(n_clusters=n)
    km.fit(x)
    max_cluster_i = np.argmax(km.cluster_centers_)
    max_cluster = [x[i]for i in range(len(km.labels_)) if km.labels_[i] == max_cluster_i ]
    thresh_val = np.min(max_cluster)
    
    if plot_bool or save_file != '':
        bts = pers_array[pers_array[:,0]==h].T[1]
        dts = pers_array[pers_array[:,0]==h].T[2]
        bts = bts[np.isfinite(dts)]

        fig = plt.figure()
        plt.scatter(bts,x, c=km.labels_, cmap='rainbow')
        plt.title('H%i Persistence Threshold = %2.f' %(h,thresh_val))
        plt.ylabel('Lifetime (Birth-Death)')
        plt.xlabel('Birth Level')
        plt.tight_layout()
        

        if save_file != '':
            plt.savefig(save_file+"KMeansPersThreshH%i.pdf" %h,dpi=300)
        if plot_bool: 
            plt.show()
        else: plt.close()
    
    return(thresh_val)




