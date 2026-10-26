#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 10:20:28 2024

@author: saratahv
"""

import numpy as np
from skimage.measure import find_contours
from scipy.stats import circmean
from sklearn.cluster import DBSCAN
from numpy.lib.stride_tricks import sliding_window_view

# Victoria's function for extending the data grid at both ends
def wrap_data(lons_in, 
              data_in,
              hh # number of grid points to add on each side
              ):
    nlon_in=lons_in.size
    f_periodic = np.zeros([data_in.shape[0],data_in.shape[1]+hh*2])
    f_periodic[:,hh:-hh], f_periodic[:,0:hh], f_periodic[:,nlon_in+hh:nlon_in+hh*2] = data_in, data_in[:,nlon_in-hh:nlon_in], data_in[:,0:hh]
    
    lons_out = np.zeros([nlon_in+hh*2])
    lons_out[hh:-hh], lons_out[0:hh], lons_out[nlon_in+hh:nlon_in+hh*2] = lons_in, lons_in[nlon_in-hh:nlon_in]-360, lons_in[0:hh]+360
    
    return lons_out,f_periodic

# convert longitudes outside 0,360 degree range (the extensions) to be within that range
def convert_extended_longitude(lon):
    # adding 360 degrees to negative longitudes
    if lon < 0:
        lon+=360
   # subtracting 360 degrees from longitudes larger than 360
    if lon >= 360:
        lon -= 360
    return lon


# transform data indices to data grid's longitudes and latitudes
# (for output of find_contours)
def xy_to_ll(x,y,lons,lats):
    # round data coordinate to nearest int 
    # then select the lon/lat coordinate at that index
    lon = np.take(lons, np.rint(x).astype(int))
    lat = np.take(lats, np.rint(y).astype(int))
    return lon, lat


# calculate difference of longitudes ranging from 0 to 360
# (most likely also works for range -40 to 400)
def lon_diff(lon1, lon2, abs_diff = True):
    diff = lon2 - lon1
    # ensure that the result of subtraction is the smaller angle
    diff = (diff + 180) % 360 - 180
    if abs_diff:
        return abs(diff)
    else:
        return diff


# calculate the large circle distance between two coordinates
# uses the Haversine formula which isn't exact
# used in find_syopscale_contours only
def haversine_dist(lon1, lat1, lon2, lat2):
    dlon = np.deg2rad(lon_diff(lon2, lon1))
    dlat = np.deg2rad(np.abs(lat2 - lat1))
    
    lon1 = np.deg2rad(lon1)
    lon2 = np.deg2rad(lon2)
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    distance = 2 * np.arcsin(np.sqrt(a)) * 6371
    return distance

# looks for contours of given isentropes that circle the globe
def find_open_contours(tpot, isentropes, 
                         longitudes, latitudes):
    open_contours = []
    
    for i in isentropes:
        contours = find_contours(tpot,i) # returns a list of contours (as lists of data indices)
        for c in contours:
            # separate x and y coordinates of each contour
            x = np.array(c[:,1])
            y = np.array(c[:,0])
            
            # convert contour coordinates (indices) to latitudes and longitudes
            lo, la = xy_to_ll(x,y,
                              longitudes,latitudes)
            
            # exclude closed contours
            # assumption: contour open if indices of min and max longitudes in contour
            if ((0 in x) & (len(longitudes)-1 in x)):
                # extra check: the contour also needs to be at least as long
                # as the x direction/the longitude array
                if len(x) >= len(longitudes):
                    open_contours.append(np.array([lo,la]))
    return open_contours

# finds instances of longitudes repeating in open contours, these are classified as overturning contours
# returns a list of coordinates where contours turn over (2d lists of latitudes and longitudes)
# will in actuality contain duplicates since extended areas are in the data twice but with different longitudes
def find_overturning_contours(open_contours):
    overturning_contours = []
    
    # iterating over open contours
    for o_c in open_contours:
        contour_lons = o_c[0]
        contour_lats = o_c[1]
        
        ot_coords = [[],[]]
        # iterating over each longitude in the open contour
        for i,l in enumerate(contour_lons):
            # check if the current longitude occurs in the contour more than once
            lon_occurrences = np.count_nonzero(l==contour_lons)
            if lon_occurrences > 1:
                # checking that there actually is overturning, i.e. that all the instances of the repeated longitude aren't consecutive
                # in practice checking that the indices of the coordinates where longitude repeats are not in consecutive order
                # get repeating longitude indices
                lon_indices = np.where(contour_lons == l)
                # get the difference of indices: if all 1, then they're consecutive
                if np.all(np.diff(lon_indices)==1):
                    # if the indices are all consecutive, this is not a part of an overturning
                    # if no overturning found but the list containing the most
                    # recent overturning is not empty, returning it and making a new one
                    if len(ot_coords[0]) > 0:
                        overturning_contours.append(np.array(ot_coords))
                        ot_coords = [[],[]]
                    continue
                # if a longitude is in the data more than twice
                # (and the indices of the longitudes are not consecutive)
                if lon_occurrences > 2:
                    ot_coords[0].append(l)
                    ot_coords[1].append(contour_lats[i])
                # in contour turning points the longitude may be repeated only twice
                # however this can also occur when the contour is along a meridian
                # so these points are added to an overturning only when
                # an overturning has been found already (list of overturning coordinates longer than 0)
                elif (lon_occurrences == 2) and (len(ot_coords[0]) > 0):
                    ot_coords[0].append(l)
                    ot_coords[1].append(contour_lats[i])
            else:
                # if no overturning found but the list containing the most
                # recent overturning is not empty, returning it and making a new one
                if len(ot_coords[0]) > 0:
                    overturning_contours.append(np.array(ot_coords))
                    ot_coords = [[],[]]
    return overturning_contours


# this function has not been implemented and still contains some bugs
# goes over every overturning contour and removes duplicates
# based on checking whether a contour is a subset of other contours
# ideally would be a more accurate alternative for remove_repeating_center_points
# def keep_unique_overturnings(overturning_contours):
#     # list where unique overturnings will be added
#     unique_overturning_contours = []
#     # iterate over all found overturning contours
#     for oc1 in overturning_contours:
#         contour_is_unique = True
#         # # extract the coordinates of the starting point of this contour
#         # lon1 = oc1[0,0]
#         # lat1 = oc1[1,0]
#         # iterating over the current list of unique overturnings
#         for i, oc2 in enumerate(unique_overturning_contours):
#             # if contour from first level of iteration (oc1) is longer,
#             # looking for subsets of second contour inside it
#             if len(oc1[0]) > len(oc2[0]):
#                 windows = sliding_window_view(oc1, window_shape = oc2.shape)
#                 oc2_is_subarray = np.any(np.all(windows == oc2, axis = 1)) # is the axis correct? the example uses 1d arrays
#                 # if oc2 is a subarray of oc1 then it needs to be 
#                 # removed from the unique list and oc1 added instead
#                 if oc2_is_subarray: 
#                     del unique_overturning_contours[i]
#                     print("throwing out contour in unique list")
#                 # oc1's status as unique remains true and it will be added to the unique list

#             # if contour from second level of iteration/"unique" contour is longer,
#             # looking for subsets of first contour inside it
#             elif len(oc1[0]) < len(oc2[0]):
#                 windows = sliding_window_view(oc2, window_shape = oc1.shape)
#                 oc1_is_subarray = np.any(np.all(windows == oc1, axis = 1))
#                 # if oc1 is a subarray of oc2 it's a shorter duplicate
#                 # -> definitely not uat if two identical contours are slightly shifted they are still recognised as duplicates)
#                 oc1_subset = oc1[:,2:-2]
#                 windows = sliding_window_view(oc2, window_shape = oc1_subset.shape)
#                 oc1_is_subarray = np.any(np.all(windows == oc1_subset, axis = 1))
#                 # 
#                 if oc1_is_subarray: 
#                     print("throwing out contour of equal length to other contour")
#                     # contour_is_unique = False

#         if contour_is_unique:
#             unique_overturning_contours.append(oc1)
#     return unique_overturning_contoursnique, will not be added to the unique list
#                 if oc1_is_subarray: 
#                     print("Not adding a contour because it's a subset")
#                     contour_is_unique = False

#             # what if the contours are the exact same length? (presumably the most common case with duplicates)
#             elif len(oc1[0]) == len(oc2[0]):
#                 # taking a subset of oc1 and checking if that is a subset of oc2
#                 # (this makes sure that if two identical contours are slightly shifted they are still recognised as duplicates)
#                 oc1_subset = oc1[:,2:-2]
#                 windows = sliding_window_view(oc2, window_shape = oc1_subset.shape)
#                 oc1_is_subarray = np.any(np.all(windows == oc1_subset, axis = 1))
#                 # 
#                 if oc1_is_subarray: 
#                     print("throwing out contour of equal length to other contour")
#                     # contour_is_unique = False

#         if contour_is_unique:
#             unique_overturning_contours.append(oc1)
#     return unique_overturning_contours



# discards overturning contours not in the synoptic scale
# based on width (longitudinal, 5-40 degrees) and length (>1500 km) of contour
# excluding very wide contours is optional because it's biased towards high latitude contours
def find_synopscale_contours(contours, use_max_width_condition):
    kept_contours = []
    for c in contours:
        keep_contour = True
        
        lons = c[0,:]
        lats = c[1,:]
        
        # calculate width of the overturning contour
        # (distance beteen the first and last point)
        # excluding widths smaller than 5 degs lon (very small or very meridional contour)
        if lon_diff(c[0,-1], c[0,0]) < 5.0:
            keep_contour = False
        
        # excluding widths larger than 40 degrees lon (very wide)
        # this is applied only when using Bowley's parameters for the algorithm
        if use_max_width_condition:
            if lon_diff(c[0,-1], c[0,0]) > 40.0:
                keep_contour=False
        
        # checking the length of the contour
        if keep_contour:
            # calculating the length of the contour by integrating the Haversine formula over contour coordinates
            contour_length = np.sum(haversine_dist(lons[:-1], lats[:-1],
                                    lons[1:], lats[1:]))
            # keeping only contours longer than 1500 km
            if contour_length < 1500.0: keep_contour = False
        if keep_contour: kept_contours.append(c)
    
    return kept_contours


# calculating the center points for the contours
# (mean of all coordinates in contour)
def get_center_points(overturning_contours):
    # array for the output coordinates
    center_points = np.zeros((len(overturning_contours),2))
    for i, c in enumerate(overturning_contours):
        c_lon = np.deg2rad(c[0,:])
        # longitude requires using circmean
        # converts the longitude to be within [0,360] I think
        center_points[i,1] = np.rad2deg(circmean(c_lon))
        # mean latitude with regular mean
        center_points[i,0] = np.mean(c[1,:])
    return center_points


# removing duplicate overturning contours from the extensions by comparing center points
# if two center points are within 0.1 degrees in both latitude and longitude
# they're classified as duplicates and only one of the center points and corresponding contours is kept
def remove_repeating_center_points(center_points, overturnings):
    unique_center_points = []
    unique_overturnings = []
    # iterating over all center points and contours
    for p1, contour in zip(center_points, overturnings):
        point_is_unique = True
        p1_lon = convert_extended_longitude(p1[1]) # just to make sure
        p1_lat = p1[0]
        
        # iterating over unique center points
        for p2 in unique_center_points:
            p2_lon = convert_extended_longitude(p2[1])
            p2_lat = p2[0]            
            # if difference between coordinates is small enough, point is not unique
            if (lon_diff(p1_lon, p2_lon) < 0.1 and abs(p1_lat - p2_lat) < 0.1):
                point_is_unique = False

        # if no point is close enough or no points are yet in the list, the current point and overturning are unique
        if point_is_unique:
            unique_center_points.append(p1)
            unique_overturnings.append(contour)
    return np.array(unique_center_points), unique_overturnings


# clustering the center points based on distance to other points
# (at most 1000 or 1667 km) and numer of points within that distance (at least 2 or 3 points)
# values need to be in radians for DBSCAN using haversine (large circle) distance
def cluster_center_points(center_points, max_dist, min_samples):
    # converting everything to radians
    center_points_radians = np.deg2rad(center_points)
    epsilon = max_dist/6371.0
    # creating the clusterer object
    clusterer = DBSCAN(
        eps = epsilon, 
        min_samples = min_samples, 
        metric="haversine")
    # fitting data with the clusterer, extracting labels denoting clusters
    # for each center point
    labels = clusterer.fit(center_points_radians).labels_
    return labels


# separating the contours by the clusters they're in 
# (and discarding contours outside clusters)
# (needs list comprehension again because clusters can and do have 
# different numbers of contours and the contours themselves are all of different lengths)
def group_contours_by_labels(overturning_contours, labels):
    contour_groups = []
    # iterating over values in labels
    for label in range(max(labels)+1):
        group = []
        # extracting indices with labels and indexing overturning_contours to add contours to groups
        label_indices = np.where(labels==label)[0]
        for l in label_indices:
            group.append(overturning_contours[l])
        contour_groups.append(group)
    return contour_groups


# calculating RWB breaking direction for one contour
# by comparing latitudes of eastern- and westernmost points of the contour
def get_DB(contour):
    # for some reason the coordinates are in reverse order in the SH
    # fixing that here so that all coordinate arrays go from east to west and DB is calculated correctly
    # first latitude is negative -> contour in SH (not a perfect condition but we're not looking for waves crossing the equator)
    if contour[1,0] < 0: 
        contour = contour[:,::-1]
    start_point = contour[0,0], contour[1,0] # longitude,latitude
    end_point = contour[0,-1], contour[1,-1]
        
    # comparing latitudes of start and end points
    # cyclonic wave breaking when westernmost point
    # nearer to equator than easternmost point
    if abs(start_point[1]) < abs(end_point[1]):
        return 1.0 # for cyclonic wave breaking
    else:
        return -1.0 # for anticyclonic wave breaking
    
# implementing above function for a list of contours
# separates cyclonically and anticyclonically turning contours
# so that they can then be clustered separately
def calculate_DBs(list_of_contours, list_of_center_points):
    cyclonic_contours = []
    cyclonic_center_points = []
    anticyclonic_contours = []
    anticyclonic_center_points = []
    for c, p in zip(list_of_contours, list_of_center_points):
        db = get_DB(c)
        if db == 1:
            cyclonic_contours.append(c)
            cyclonic_center_points.append(p)
        else:
            anticyclonic_contours.append(c)
            anticyclonic_center_points.append(p)
    return cyclonic_contours, cyclonic_center_points, anticyclonic_contours, anticyclonic_center_points

# this function is made redundant by calculating DBs before clustering
# # Defining a DB for each group of contours
# # Special cases when all contours in group don't break in the same direction
# def calculate_DB(contour_groups, min_contours):
#     DBs = []
#     groups_to_return = []
#     # iterating over contour groups
#     for group in contour_groups:
#         DBs_in_wave = []
#         # calculating DB for each contour in the group
#         for i, c in enumerate(group):
#             DBs_in_wave.append(get_DB(c))
#         DBs_in_wave = np.array(DBs_in_wave)

#     # checking if all DBs in a group are in the same direction
#         if abs(np.sum(DBs_in_wave)) == DBs_in_wave.size:
#             # add DB and RWB with unambiguous DB to groups to return
#             DBs.append(DBs_in_wave[0])
#             groups_to_return.append(group)
            
            
#         # if group contains different DBs, 4 options:
#             # 1) remove erroneous anticyclonic turning(s) from CWB
#             # 2) remove erroneous cyclonic turning(s) from AWB
#             # 3) separate CWB and AWB (at least n=min_contours contours break in each direction)
#             # 4) discard group as it must contain fewer than min_contours contours breaking in each direction

#         # especially with unfiltered data and Bowley's parameters (large clustering distance) this sometimes results in non-ideal groups
#         else:
#             # check how many waves break in each direction
#             anticyclonic_idx = np.where(DBs_in_wave == -1)[0]
#             cyclonic_idx = np.where(DBs_in_wave==1)[0]
            
#             # adding a cyclonic wave breaking and excluding anticyclonic contours if there are not enough anticyclonic contours to consitute an RWB
#             if (anticyclonic_idx.size < min_contours and len(group) - anticyclonic_idx.size >= min_contours): 
#                 # extracting the cyclonic RWBs from the group
#                 DB_cyclonic = [group[idx] for idx in cyclonic_idx]
#                 DBs.append(1.0) # the group added is cyclonic
#                 groups_to_return.append(DB_cyclonic)
#             # adding an anticyclonic wave breaking and excluding cyclonic contours if there are not enough cyclonic contours to constitute an RWB
#             elif(cyclonic_idx.size < min_contours and len(group) - cyclonic_idx.size >= min_contours): 
#                 # extracting the anticyclonic RWBs from the group
#                 DB_anticyclonic = [group[idx] for idx in anticyclonic_idx]
#                 DBs.append(-1.0) # the group added is anticyclonic
#                 groups_to_return.append(DB_anticyclonic)
#             # if there are enough contours in both directions to have waves breaking in both directions, 
#             # separating them into two groups (quite rare but with noisy data, these may not always be within the clustering distance of each other)
#             elif(cyclonic_idx.size >= min_contours and anticyclonic_idx.size >= min_contours):
#                 DB_cyclonic = [group[idx] for idx in cyclonic_idx]
#                 DBs.append(1.0)
#                 groups_to_return.append(DB_cyclonic)
                
#                 DB_anticyclonic = [group[idx] for idx in anticyclonic_idx]
#                 DBs.append(-1.0)
#                 groups_to_return.append(DB_anticyclonic)
#             # This leaves only groups where there are not enough waves breaking in either direction to form a group
#             # they are discarded (not added to the return list)
#             else:
#                 continue
#     return groups_to_return, DBs



# main function to find breaking Rossby waves from a potential temperature at 2 PVU
def find_rwb(tpot, isentropes,
             longitudes,latitudes,
             min_contours_in_overturning, # 2 for me, 3 for Bowley
             max_clustering_distance, # 1000 (km) for me, 1667 for Bowley
             use_max_width_condition = False): # if True, keep_synopscale_contours will remove contours wider than 40 degrees longitude (Bowley's implementation)
    
    # adding 40 degrees to both ends of the TPOT field
    longitudes, tpot_extended = wrap_data(longitudes, tpot, 59) # 59 grid points = about 40 degrees

    # looking for open contours from the extended field
    open_contours_extended = find_open_contours(tpot_extended, 
                                           isentropes, 
                                           longitudes, latitudes)

    # looking for overturnings in the contours
    overturning_contours_extended = find_overturning_contours(open_contours_extended)

    # converting coordinates to be within [0,360]
    # for removing duplicate contours based on comparing coordinates
    # overturning_contours = convert_extended_coordinates(overturning_contours_extended)
    overturning_contours = overturning_contours_extended
    
    # keep only overturning contours in the synoptic scale
    # good to do before removing duplicates because that is an expensive function
    # if not using Bowley et al.'s implementation set use_max_width_condition to False (default)
    synopscale_overturnings = find_synopscale_contours(overturning_contours, use_max_width_condition)
    
    # calculate center points for each contour
    center_points = get_center_points(synopscale_overturnings)
    # remove repeating center points
    unique_center_points, unique_synopscale_overturnings = remove_repeating_center_points(center_points, synopscale_overturnings)
    
    # calculate DB for each contour so that clustering can be done separately for each DB
    cyclonic_contours, cyclonic_center_points, anticyclonic_contours, anticyclonic_center_points = calculate_DBs(unique_synopscale_overturnings, unique_center_points)
    
    # cluster center points by distance using DBSCAN
    # my original arguments for DBSCAN: 
        # 1) distance between nearest points in cluster at most 1000 km (max_clustering_distance)
        # 2) 2 points within 1000 km required to make a cluster (min_contours_in_overturning)
    # Bowley et al.'s implementation:
        # 1) distance between nearest points in cluster at most 1667 km
        # 2) 3 points within 1667 km required to make a cluster
    cyclonic_labels = cluster_center_points(cyclonic_center_points, max_clustering_distance, min_contours_in_overturning)
    anticyclonic_labels = cluster_center_points(anticyclonic_center_points, max_clustering_distance, min_contours_in_overturning)
    # group the contours themselves with the labels from clustering
    cyclonic_contour_groups = group_contours_by_labels(cyclonic_contours, cyclonic_labels)
    anticyclonic_contour_groups = group_contours_by_labels(anticyclonic_contours, anticyclonic_labels)
    
    # adding the two DBs and contour groups together because later functions want a single list
    # (DB corresponding to group still at same index so no info lost)
    DBs = [1 for i in range(len(cyclonic_contour_groups))]
    DBs.extend([-1 for j in range(len(anticyclonic_contour_groups))])
    
    contour_groups = cyclonic_contour_groups
    contour_groups.extend(anticyclonic_contour_groups)

    return contour_groups, DBs, tpot_extended, longitudes