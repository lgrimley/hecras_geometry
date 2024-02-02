#!/usr/bin/env python3

import parserasgeo as prg
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def linestring_to_points(polyline_shp):
    # This function reads in a shapefile as a geopandas dataframe and loops reads the geometry linestring as points
    # The output geodataframe has a column called 'points' with the coordinates
    gdf = gpd.read_file(polyline_shp)
    gdf['points'] = gdf.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    return gdf


def get_ras_xs_info(geo):
    # This function reads a HEC-RAS geometry file (.g01) and returns information for each cross-section
    # Empty lists where data is stored after being read from the geometry file.
    # These lists are combined into a dataframe and output
    stream_stn = list()
    chn_width = list()
    l_bnk_elev = list()
    r_bnk_elev = list()
    bed_elev = list()

    for xs in geo.get_cross_sections():
        l_bnk_sta = xs.bank_sta.left # Get the Left Bank Station
        r_bnk_sta = xs.bank_sta.right # Get the Right Bank Station
        ch_width = r_bnk_sta - l_bnk_sta # Calculate channel width
        lbe = xs.sta_elev.elevation(sta=l_bnk_sta) # Get the Left Bank elevation
        rbe = xs.sta_elev.elevation(sta=r_bnk_sta) # Get the right Bank elevation

        # Get station and elevation of the cross-section points
        xs_pts = pd.DataFrame(data=None)
        xs_pts['sta'] = [x[0] for x in xs.sta_elev.points]
        xs_pts['elev'] = [x[1] for x in xs.sta_elev.points]

        # Subset cross-section points to those within the channel
        channel_pts = xs_pts[xs_pts.sta >= l_bnk_sta]
        channel_pts = channel_pts[channel_pts.sta <= r_bnk_sta]
        min_bed = channel_pts.elev.min() # Get the minimum bed elevation in the channel

        # Append to list
        stream_stn.append(xs.header.station.value)  # Get the RAS Stream Station
        chn_width.append(ch_width)
        l_bnk_elev.append(lbe)
        r_bnk_elev.append(rbe)
        bed_elev.append(min_bed)

        # Combine into dataframe
    xsdf = pd.DataFrame()
    xsdf['STREAM_STN'] = stream_stn
    xsdf['chn_width'] = chn_width
    xsdf['l_bnk_elev'] = l_bnk_elev
    xsdf['r_bnk_elev'] = r_bnk_elev
    xsdf['bed_elev'] = bed_elev

    # Sort the dataframe based on Stream Station
    xsdf_sort = xsdf.sort_values(by=['STREAM_STN'], ignore_index=True)

    return xsdf_sort


def convert_to_metric(df):
    conversion = 0.3048
    df['chn_width_m'] = df['chn_width'] * conversion
    df['l_bnk_elev_m'] = df['l_bnk_elev'] * conversion
    df['r_bnk_elev_m'] = df['r_bnk_elev'] * conversion
    df['bed_elev_m'] = df['bed_elev'] * conversion
    return df


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


def convert_list_to_tuple(list):
    return tuple(x for x in list)


# Read in HEC-RAS Geometry File using Parserasgeo package
geo = prg.ParseRASGeo(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/ras_geofiles/011003_NEW_RIVER.g02')
df = get_ras_xs_info(geo)
df = convert_to_metric(df)


#####################################
# Code to match RAS stream station to shapefile stream station (not good)
pts_shp = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/Research/NewRiver/bathy/shapefile/NewRiver_upper_Inter.shp'
pts = gpd.read_file(pts_shp)
pts = pts.sort_values(by='STREAM_STN', ignore_index=True)
pts['xcoord'] = pts.geometry.centroid.x.iloc[:]
pts['ycoord'] = pts.geometry.centroid.y.iloc[:]

df_match = find_nearest_xs(df1=df, colname1='STREAM_STN', df2=pts, colname2='STREAM_STN') # Need to find a way to join df and pts after finding the best match

help1 = df.set_index('STREAM_STN', drop=True, inplace=False)
help2 = pts.set_index('STREAM_STN', drop=True, inplace=False)

#####################################
# Code to read RAS stream centerline shapefile and write cutline to HEC-RAS geometry file
stream_file = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/Research/NewRiver/bathy/shapefile/S_STREAMCNTRLINE_NR.shp'
st = linestring_to_points(stream_file)
st = st[st['HYDRAID'] == '3710402014011232']  # '3710402014011230'

reach_cutline = pd.DataFrame(data=None, columns=['Schematic_X', 'Schematic_Y'])
for reach in geo.get_reaches():
    reach_name = reach.header.reach_name
    river_name = reach.header.river_name
    geom = st['points'][:].tolist()[0]

    # [(x1,y1),(x2,y2),(x3,y3),...] Values are currently stored as strings
    geom_list = []
    for ii in range(len(geom)):
        reach_cutline = reach_cutline.append({'Schematic_X': geom[ii][0], 'Schematic_Y': geom[ii][1]}, ignore_index=True)
        #coord = list(geom[ii])
        #coord[0] = str(round(coord[0], 2))
        #coord[1] = str(round(coord[1], 2))
        #coord_tuple = convert_list_to_tuple(coord)
        #geom_list.append(coord_tuple)

    #reach.geo.points = geom_list

reach_cutline.to_csv('/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/Research/NewRiver/bathy/011003_NEW_RIVER.g01.reach.csv')

#####################################
# Code to read RAS xs shapefile and write cutline to HEC-RAS geometry file
shp_file = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/Research/NewRiver/bathy/shapefile/S_HYDRACROSSSECTION_NR.shp'
gdf_base = linestring_to_points(polyline_shp=shp_file)
gdf_base['STREAM_STN2'] = round(gdf_base.STREAM_STN)
gdf = gdf_base[gdf_base['HYDRAID'] == '3710402014011230']  # '3710402014011232'

xs_cutline = pd.DataFrame(data=None, columns=['Schematic_X', 'Schematic_Y'])

missing = []
for xs in geo.get_cross_sections():
    # Stream Station
    stream_stn = round(np.float64(xs.header.station.value))
    xs_geom = gdf[gdf['STREAM_STN2'] == stream_stn]['points'][:]
    #if xs_geom.empty is True:
        #print('No match for STREAM_STN', stream_stn)
        #missing.append(np.float64(xs.header.station.value))
    #else:
        xs_geom = list(xs_geom.items())
        points = xs_geom[0][1]

        # [(x1,y1),(x2,y2),(x3,y3),...] Values are currently stored as strings
        geom_list = []
        for ii in range(len(points)):
            coord = list(points[ii])
            coord[0] = str(round(coord[0], 2))
            coord[1] = str(round(coord[1], 2))
            if ii == 0:
                coord[0] = '10000'
            coord_tuple = convert_list_to_tuple(coord)
            geom_list.append(coord_tuple)

        #xs.cutline.points = geom_list
        break

geo.write('/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/Research/NewRiver/bathy/test.g01')

