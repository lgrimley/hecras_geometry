#!/usr/bin/env python3

import os
import parserasgeo as prg
import pandas as pd
import geopandas as gpd


def get_filenames(dir):
    f = list()
    for filename in os.listdir(dir):
        if filename.endswith('.g01'):
            f.append(filename)
    return f


def get_hydraid(filename):
    hydraid = filename.split('_')[0]
    return hydraid


def get_ras_xs_info(geo):
    stream_stn = list()
    chn_width = list()
    l_bnk_elev = list()
    r_bnk_elev = list()
    bed_elev = list()

    for xs in geo.get_cross_sections():
        # Stream Station
        stream_stn.append(xs.header.station.value)

        # Left and Right Bank Station
        l_bnk_sta = xs.bank_sta.left
        r_bnk_sta = xs.bank_sta.right

        # Calculate channel width
        chn_width.append(r_bnk_sta - l_bnk_sta)

        # Left and Right Bank Elevation
        l_bnk_elev.append(xs.sta_elev.elevation(sta=l_bnk_sta))
        r_bnk_elev.append(xs.sta_elev.elevation(sta=r_bnk_sta))

        # Get station and elevation of the cross-section points
        xs_pts = pd.DataFrame(data=None)
        xs_pts['sta'] = [x[0] for x in xs.sta_elev.points]
        xs_pts['elev'] = [x[1] for x in xs.sta_elev.points]

        # Subset cross-section points to those within the channel
        channel_pts = xs_pts[xs_pts.sta >= l_bnk_sta]
        channel_pts = channel_pts[channel_pts.sta <= r_bnk_sta]

        # Get the minimum bed elevation in the channel
        bed_elev.append(channel_pts.elev.min())

        # Combine into dataframe
        xsdf = pd.DataFrame()
        xsdf['STREAM_STN'] = stream_stn
        xsdf['chn_width'] = chn_width
        xsdf['l_bnk_elev'] = l_bnk_elev
        xsdf['r_bnk_elev'] = r_bnk_elev
        xsdf['bed_elev'] = bed_elev

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



dir = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/Research/NewRiver/bathy/ras_geom'
filenames = get_filenames(dir)
f = filenames[0]

# Read the HECRAS Geometry File
hydraid = get_hydraid(f)
geo = prg.ParseRASGeo(os.path.join(dir, f))
df = get_ras_xs_info(geo)
df = convert_to_metric(df)

# Read in Stream Shapefile
stream_file = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/Research/NewRiver/bathy/shapefile/S_STREAMCNTRLINE_NR.shp'
stream = gpd.read_file(stream_file)
stream_subset = stream[stream['HYDRAID'] == hydraid]

# Read in XS Shapefile
xs_file = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/Research/NewRiver/bathy/shapefile/S_HYDRACROSSSECTION_NR.shp'
xs = gpd.read_file(xs_file)
xs_subset = xs[xs['HYDRAID'] == hydraid]

points = pd.DataFrame(data=None, columns=['STREAM_STN', 'X', 'Y'])
for stream in stream_subset.itertuples():
    for xs in xs_subset.itertuples():
        xs_geom = xs.geometry
        intersect = xs_geom.intersection(stream.geometry)
        if len(intersect.coords) != 0:
            points = points.append({'STREAM_STN': xs.STREAM_STN, 'X': intersect.x, 'Y': intersect.y}, ignore_index=True)

match = find_nearest_xs(df1=df, colname1='STREAM_STN', df2=xs_subset, colname2='STREAM_STN')
df['stn_match'] = match['DF2_Val']

points.set_index('STREAM_STN', inplace=True)
df.set_index('stn_match', inplace=True)
data = points.join(df)

fout = f.replace('.g01', '.csv')
data.to_csv(os.path.join(dir, fout))
