#!/usr/bin/env python3
import parserasgeo as prg
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def get_xs_info(xs):
    l_bnk_sta = xs.bank_sta.left  # Get the Left Bank Station
    r_bnk_sta = xs.bank_sta.right  # Get the Right Bank Station
    ch_width = r_bnk_sta - l_bnk_sta  # Calculate channel width
    lbe = xs.sta_elev.elevation(sta=l_bnk_sta)  # Get the Left Bank elevation
    rbe = xs.sta_elev.elevation(sta=r_bnk_sta)  # Get the right Bank elevation
    stream_stn = xs.header.station.value  # Get the RAS Stream Station
    xs_info = pd.DataFrame(data=np.array([[stream_stn, ch_width, rbe, lbe, r_bnk_sta, l_bnk_sta]]),
                           columns=['STREAM_STN', 'CH_WIDTH', 'R_BNK_ELEV', 'L_BNK_ELEV', 'R_BNK_STA', 'L_BNK_STA'])
    xs_info = round(xs_info, 2)
    return xs_info


def get_channel_data(xs):
    l_bnk_sta = xs.bank_sta.left  # Get the Left Bank Station
    r_bnk_sta = xs.bank_sta.right  # Get the Right Bank Station

    # Get cross-section station and elevation
    xs_pts = pd.DataFrame(data=None)
    xs_pts['sta'] = [x[0] for x in xs.sta_elev.points]
    xs_pts['elev'] = [x[1] for x in xs.sta_elev.points]

    # Subset cross-section points to those within the channel
    channel_pts = xs_pts[xs_pts.sta >= l_bnk_sta]
    channel_pts = channel_pts[channel_pts.sta <= r_bnk_sta]
    return channel_pts


def calculate_trapezoid_area(x1, y1, x2, y2, s):
    dx = abs(x2 - x1)
    if s < y1 and s < y2:
        area = 0
    elif s == y1 or s == y2:
        area = (dx * abs(y1 - y2)) / 2
    elif s > y1 and s > y2:
        area = (s - max(y1, y2)) * dx + ((s - min(y1, y2)) * dx) / 2
    else:
        area = ((s - min(y1, y2)) * dx) / 2
    return area


def calculate_wetted_channel_area(channel_pts):
    channel_pts.elev = min(channel_pts.elev) + channel_pts.elev
    area = list()
    i = 0
    while i < (len(channel_pts.sta) - 1):
        bankelev = max(channel_pts.elev.iloc[0], channel_pts.elev.iloc[-1])
        a = calculate_trapezoid_area(x1=channel_pts.sta.iloc[i], y1=channel_pts.elev.iloc[i],
                                     x2=channel_pts.sta.iloc[i + 1], y2=channel_pts.elev.iloc[i + 1],
                                     s=bankelev)
        area.append(a)
        i += 1
    area = sum(area)
    return area


geo = prg.ParseRASGeo(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/ras_geofiles/011003_NEW_RIVER.g02')

xsdf = pd.DataFrame()
plt_frequency = None
i = 0
for xs in geo.get_cross_sections():
    xs_info = get_xs_info(xs)
    channel_pts = get_channel_data(xs)
    station = channel_pts.sta.to_numpy()
    elevation = channel_pts.elev.to_numpy()
    xs_info['BED_ELEV'] = min(channel_pts.elev)
    area = calculate_wetted_channel_area(channel_pts)
    xs_info['CH_AREA'] = area
    eff_depth = area / xs_info.CH_WIDTH.iloc[0]
    xs_info['EFF_DEPTH'] = eff_depth
    max_bank_elev = max(xs_info.R_BNK_ELEV.iloc[0], xs_info.L_BNK_ELEV.iloc[0])
    depth = abs(max_bank_elev-elevation)
    xs_info['MAX_DEPTH'] = max(depth)
    xs_info['MIN_DEPTH'] = min(depth)
    xs_info['MAX_DEPTH_RECT_AREA'] = xs_info['MAX_DEPTH'] * xs_info['CH_WIDTH']
    xs_info['AVG_DEPTH'] = depth.mean()


    if i == plt_frequency:
        xxx = [station[0], station[0], station[-1], station[-1]]
        max_elev = max(elevation[0], elevation[-1])
        min_elev = min(elevation)
        rect = [max_elev, min_elev, min_elev, max_elev]
        avg_bnk = (elevation[0] + elevation[-1]) / 2
        bed_eff = avg_bnk - xs_info.EFF_DEPTH.iloc[0]
        rect_eff = [avg_bnk, bed_eff, bed_eff, avg_bnk]

        plt.figure()
        plt.plot(xxx, rect, label='Rectangle', linewidth=3)
        plt.plot(xxx, rect_eff, label='Effective Rectangle', linewidth=3)
        plt.plot(station, elevation, label='Survey', linewidth=3)
        plt.xlabel('Station')
        plt.ylabel('Elevation (ft NAVD88)')
        plt.grid()
        plt.legend()
        plt.title('Cross-section ID:' + str(xs_info.STREAM_STN.iloc[0]))
        plt.show()
        outname = '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/NR_xs' + str(
            xs_info.STREAM_STN.iloc[0]) + '.png'
        plt.savefig(outname)
        plt.close()
        i = 0
    else:
        i += 1
    xsdf = xsdf.append(xs_info)


xsdf.to_csv('/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/scripts/NR_effective_depth.csv')
# Plot Channel Width
plt.figure()
plt.plot(xsdf.STREAM_STN, xsdf.CH_WIDTH, linestyle='dashed',
         marker='o', markerfacecolor='yellow', markersize=5)
plt.xlim(max(xsdf.STREAM_STN), min(xsdf.STREAM_STN))
plt.xlabel('Stream Station')
plt.ylabel('Channel Width (ft)')
plt.grid()
plt.show()
plt.savefig('/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/NR_channelwidth.png')
plt.close()

# Plot channel bank and bed elevations
plt.figure()
plt.plot(xsdf.STREAM_STN, xsdf.BED_ELEV, label='Bed Elevation', linewidth=3)
plt.plot(xsdf.STREAM_STN, xsdf.R_BNK_ELEV, label='Right Bank', linewidth=3)
plt.plot(xsdf.STREAM_STN, xsdf.L_BNK_ELEV, label='Left Bank', linewidth=3)
plt.xlim(max(xsdf.STREAM_STN), min(xsdf.STREAM_STN))
plt.xlabel('Stream Station')
plt.ylabel('Elevation (ft NAVD88)')
plt.grid()
plt.legend()
plt.show()
plt.savefig(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/NR_channelelevations.png')
plt.close()

# Plot Channel Area Scatter
plt.figure()
plt.plot(xsdf.CH_AREA, xsdf.CH_AREA, linewidth=1, color='k')
plt.scatter(x=xsdf.MAX_DEPTH_RECT_AREA, y=xsdf.CH_AREA)
plt.xlim(0, 7000)
plt.ylim(0, 7000)
plt.ylabel('Trapezoidal Area (sq.ft.)')
plt.xlabel('Rectangular Area (sq.ft.)')
plt.grid()
plt.show()
plt.savefig(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/NR_channelarea_scatter.png')
plt.close()

# Plot Channel Area
plt.figure()
plt.plot(xsdf.STREAM_STN, xsdf.CH_AREA, label='Actual XS Area', linewidth=3)
plt.plot(xsdf.STREAM_STN, xsdf.MAX_DEPTH_RECT_AREA, label='Max Depth Rectangular Area',
         linewidth=2, linestyle=':')
plt.xlim(max(xsdf.STREAM_STN), min(xsdf.STREAM_STN))
plt.xlabel('Stream Station')
plt.ylabel('Area (sq.ft.)')
plt.grid()
plt.legend()
plt.show()
plt.savefig(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/NR_channelarea.png')
plt.close()

# Plot Depth
plt.figure()
plt.plot(xsdf.STREAM_STN, xsdf.Assumed_DEPTH, label='Max Depth', linewidth=3)
plt.plot(xsdf.STREAM_STN, xsdf.EFF_DEPTH, label='Effective Rect Depth',
         linewidth=2, linestyle=':')
plt.xlim(max(xsdf.STREAM_STN), min(xsdf.STREAM_STN))
plt.xlabel('Stream Station')
plt.ylabel('Depth (feet NAVD88)')
plt.grid()
plt.legend()
plt.show()
plt.savefig(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/NR_depth.png')
plt.close()

# Plot depth scatter
plt.figure()
plt.plot(xsdf.EFF_DEPTH, xsdf.EFF_DEPTH, label='1:1', linewidth=1, color='k')
plt.scatter(x=xsdf.MAX_DEPTH, y=xsdf.EFF_DEPTH, label='Max Depth')
plt.scatter(x=xsdf.AVG_DEPTH, y=xsdf.EFF_DEPTH, label='Avg Depth')
plt.scatter(x=xsdf.MIN_DEPTH, y=xsdf.EFF_DEPTH, label='Min Depth')
plt.xlim(0, 30)
plt.ylim(0, 30)
plt.ylabel('Observed Depth (ft)')
plt.xlabel('Effective Depth (ft)')
plt.grid()
plt.legend()
plt.show()
plt.savefig(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/NR_depth_scatter.png')
plt.close()

# Plot Area Diff vs Stream Station
plt.figure()
plt.plot(xsdf.STREAM_STN, (xsdf.CH_AREA - xsdf.MAX_DEPTH_RECT_AREA), linewidth=3)
plt.xlim(max(xsdf.STREAM_STN), min(xsdf.STREAM_STN))
plt.xlabel('Stream Station')
plt.ylabel('Trapezoidal area minus Rectangular area (sq. ft.)')
plt.grid()
plt.show()
plt.savefig(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/NR_area_diff.png')
plt.close()

# Plot Area Diff vs Ch Width
plt.figure()
plt.scatter(xsdf.CH_WIDTH, (xsdf.CH_AREA - xsdf.MAX_DEPTH_RECT_AREA))
plt.xlabel('Channel Width (ft)')
plt.ylabel('Trapezoidal area minus Rectangular area (sq. ft.)')
plt.show()
plt.savefig(
    '/Users/laurengrimley/OneDrive - University of North Carolina at Chapel Hill/topobathy/NR_area_diff_vs_chnwidth.png')
plt.close()