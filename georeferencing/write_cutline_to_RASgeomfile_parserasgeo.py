#!/usr/bin/env python3

#####################
# Author: Lauren Grimley
# Contact: lauren.grimley@unc.edu
# Date: 8/17/21
# This script is a starting point for georeferencing HEC-RAS models using toolbox parserasgeo
#####################

import parserasgeo as prg
import geopandas as gpd
import pandas as pd
import numpy as np
import os


def linestring_to_points(polyline_shp):
    # This function reads in a shapefile as a geopandas dataframe and loops reads the geometry linestring as points
    # The output geodataframe has a column called 'points' with the coordinates
    gdf = gpd.read_file(polyline_shp)
    gdf['points'] = gdf.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    return gdf


def convert_list_to_tuple(list):
    return tuple(x for x in list)


# Set working directory
os.chdir()

# Load HEC-RAS geometry file
geo = prg.ParseRASGeo(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/ras_geofiles/011003_NEW_RIVER.g02')

# Load shapefile line of the stream network
st = linestring_to_points('S_STREAMCNTRLINE_NR.shp')
st = st[st['HYDRAID'] == '3710402014011232'] # Subset the shapefile by HYDRAID (NC state naming convention for HEC-RAS models)

# Loop through the geometry file reaches and write the X & Y cutline
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

# Write the cuteline to CSV to review
reach_cutline.to_csv('011003_NEW_RIVER.g01.reach.csv')

#--------------------------------#
# Loop through the geometry file XS and write the X & Y cutline
gdf_base = linestring_to_points(polyline_shp='S_HYDRACROSSSECTION_NR.shp')
gdf_base['STREAM_STN2'] = round(gdf_base.STREAM_STN)
gdf = gdf_base[gdf_base['HYDRAID'] == '3710402014011230']

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

geo.write('test.g01')

