#!/usr/bin/env python3

#####################
# Author: Lauren Grimley
# Contact: lauren.grimley@unc.edu
# Date: 8/17/21
# This script tries to match the STREAM STN name for the XS in the
# georeferenced shapefile to the data stored in the RAS geomfile
#####################

import geopandas as gpd
import pandas as pd
import os

def find_nearest_xs(df1, colname1, df2, colname2):
    sta_match = pd.DataFrame(columns=['DF1_Val', 'DF2_Ind', 'DF2_Val', 'Val_DIF'])
    for i in range(len(df1[colname1])):
        value = df1[colname1][i]
        # Find closest number to value in df2, colname 2
        nearest = df2.iloc[(df2[colname2] - value).abs().argsort()[:1]]
        # Return index, value,
        ind = nearest.index.values[0]
        value2 = nearest[colname2].values[0]
        dif = value - value2
        # Create list of matching info for output
        d = [value, ind, value2, dif]
        sta_match = sta_match.append(pd.DataFrame([d], columns=sta_match.columns), ignore_index=True)
    return sta_match


os.chdir('/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/')

# Read in HEC-RAS XS Geom Info generated from calculate_xs_area.py
xsgeo = pd.read_csv('NR_effective_depth.csv')

# XS Point shapefile where X & Y coords are the XS intersection with stream centerline
pts_shp = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/ras_geofiles/points/011003_NEW_RIVER.shp'
pts = gpd.read_file(pts_shp)
pts = pts.sort_values(by='STREAM_STN', ignore_index=True)
pts['xcoord'] = pts.geometry.centroid.x.iloc[:]
pts['ycoord'] = pts.geometry.centroid.y.iloc[:]

# Try to match the XS info read from the RAS geofile to the georeferenced XS data
# Need to find a way to join df and pts after finding the best match
df_match = find_nearest_xs(df1=xsgeo, colname1='STREAM_STN', df2=pts, colname2='STREAM_STN')
