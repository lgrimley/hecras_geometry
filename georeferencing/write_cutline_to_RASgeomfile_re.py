#!/usr/bin/env python3

#####################
# Author: Lauren Grimley
# Contact: lauren.grimley@unc.edu
# Date: 8/17/21
# This script is a starting point for georeferencing HEC-RAS models using Regular Expression
#####################

import re
import geopandas as gpd
import pandas as pd


# Functions

def linestring_to_points(polyline_shp):
    gdf = gpd.read_file(polyline_shp)
    gdf['points'] = gdf.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    return gdf


def get_cross_sections(geom_file):
    xs_stations = pd.DataFrame(data=None, columns=['STREAM_STN'])
    with open(geom_file, 'r') as geo:
        for line in geo:
            if "Type RM Length" in line:
                station = re.findall(r"[a-zA-Z\s]+=\s{0,100}\d+\.?\d{0,100}?\s{0,100},\s{0,100}?(\d+)", line)
                stream_sta = float(station[0])
                xs_stations = xs_stations.append({'STREAM_STN': stream_sta}, ignore_index=True)
    return xs_stations


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


# HECRAS Geometry File
geom_in = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/NewRiver/bathy/ras_geom/011003_NEW_RIVER.g02'
geom_out = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/NewRiver/bathy/ras_geom/011003_NEW_RIVER_georef.g02'

# Centerline
stream_file = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/NewRiver/bathy/shapefile/S_STREAMCNTRLINE_NR.shp'
st = linestring_to_points(stream_file)
st = st[st['HYDRAID'] == '3710402014011232']
st_num_pts = len(st.points.values[0])

# Cross-section Shapefile
xs_shpfile = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/NewRiver/bathy/shapefile/S_HYDRACROSSSECTION_NR.shp'
xs_gdf = linestring_to_points(polyline_shp=xs_shpfile)
xs_sub = xs_gdf[xs_gdf['HYDRAID'] == '3710402014011232']  # '3710402014011230'

xs_stream_sta = get_cross_sections(geom_file=geom_in)
xs_match = find_nearest_xs(xs_stream_sta, 'STREAM_STN', xs_sub, 'STREAM_STN')
xs_match.set_index('DF1_Val', drop=True, inplace=True)
xs_match.to_csv(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/NewRiver/bathy/match.csv')
skip_please = False
with open(geom_in, 'r') as geo:
    with open(geom_out, 'w') as w:
        for line in geo:
            if "Reach XY" in line:
                w.write('Reach XY=' + str(st_num_pts) + '\n')

                coords = st['points'].values[0]
                num_coords = st_num_pts
                for ind, val in enumerate(coords):
                    x = format(round((val[0]), 5), '.2f')
                    y = format(round((val[1]), 5), '.2f')
                    newline2 = "      " + str(x) + "       " + str(y)
                    if ind == num_coords:
                        w.write(newline2 + '\n')
                    else:
                        if ind % 2 == 0:
                            w.write(newline2)  # + "     ")
                        if ind % 2 == 1:
                            w.write(newline2 + '\n')
                w.write('\n')
                w.write('Rch Text X Y=' + '\n')
                w.write('Reverse River Text=0' + '\n')
                skip_please = True

            if "Type RM Length" in line:
                skip_please = False
                w.write('\n')
                station = re.findall(r"[a-zA-Z\s]+=\s{0,100}\d+\.?\d{0,100}?\s{0,100},\s{0,100}?(\d+)", line)
                station = float(station[0])
                join = xs_match[xs_match.index == station]['DF2_Ind']
                coords = xs_sub['points'][join.values[0]]

            if 'END DESCRIPTION' in line:
                w.write(line)
                newline1 = 'XS GIS Cut Line' + '=' + str(len(coords)) + "\n"
                w.write(newline1)

                num_coords = len(coords)
                for ind, val in enumerate(coords):
                    x = format(round((val[0]), 5), '.2f')
                    y = format(round((val[1]), 5), '.2f')
                    newline2 = "      " + str(x) + "       " + str(y)
                    if ind == num_coords:
                        w.write(newline2 + '\n')
                    else:
                        if ind % 2 == 0:
                            w.write(newline2)  # + "     ")
                        if ind % 2 == 1:
                            w.write(newline2 + '\n')
                w.write('\n')
            if skip_please is False:
                w.write(line)
