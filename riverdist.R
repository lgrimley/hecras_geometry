library(riverdist)

setwd('Z:\\users\\lelise\\data\\topobathy\\hecras_bathy\\Bachelor Creek')
load_network = TRUE
river_layer = 'stream'
point_layer = 'points_bed'

# -------------- Functions --------------
loess_smoothing <- function(x, y, alpha, interpolate){
  df = cbind(x, y)
  vars = colnames(df)
  
  colnames(df)= c('x','y')
  smoothed=predict(loess(y ~ x,
                         data=df, 
                         span=alpha,
                         control=loess.control(surface="direct")), interpolate)
  
  png('interpolated_bed_elevations.png')
  plot(df$y, x=df$x, type="p",col='gray', 
       xlab = vars[1], ylab=vars[2])
  lines(smoothed,x=interpolate, col="red", lwd=1)
  dev.off()
  
  return(smoothed)
}

calculate_slope <- function(distance, elevation, window){
    slope=c()
    i = 1
    total_distance=0
    while (total_distance <= window){
      d = distance[i+1]-distance[i]
      s = abs((elevation[i+1]-elevation[i]))/d
      slope=c(slope,s)
      i = i+1
      total_distance=total_distance+d
    }
    
    avg_slope = mean(slope) 
    bed_elev_start = (distance[1]*avg_slope)+elevation[1]
    return(bed_elev_start)
    }

# --------------  River Network --------------
proj4 = '+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '

if (load_network == FALSE){
  
  # Import stream Shapefile as a network
  network = line2network(path = '.', layer=river_layer, reproject=proj4)
  
  # Cleanup network
  streams = cleanup(network) 
  
  # Check network connectivity, clean up network before saving
  topologydots(rivers=streams)
  removeduplicates(rivers=streams)
  removeunconnected(rivers=streams)
  buildsegroutes(rivers=streams, lookup=TRUE, verbose=TRUE)
  sequenceverts(rivers=streams)
  
  save(streams, file=paste0(river_layer,'_network.Rdata'))
  
} else{
  
  load(file=paste0(river_layer,'_network.Rdata'))
}

# --------------  Add Point Data to Network --------------

# Check network connectivity
topologydots(rivers=streams)


# Import and transform point data using xy2segvert() or ptshp2segvert()
stations = pointshp2segvert(path='.',layer=point_layer, rivers=streams)

# Plot network specify segment

# Snapping distance
hist(stations$snapdist, main='Snapping Distance to Stream (m)')


# Plot network and points
zoomtoseg(seg=1, rivers=streams)
riverpoints(seg=stations$seg, vert=stations$vert, rivers=streams,
            pch=15, col='blue', cex=0.6)

# Tool to detect network route, specify start and end segments
detectroute(start=1, end=1, rivers=streams, verbose=TRUE)
showends(seg=1, rivers=streams)

startseg = 1
startvert = 1
endseg = 1
endvert = 2416 # Use plot generated above to identify end vertex

# Calculate river distances along route
riverdistance(startseg=startseg, startvert=startvert, endseg=endseg, endvert=endvert, rivers=streams, map=TRUE)


# Computing network distance
dist=c()
i=1
while (i <= length(stations$vert)){
  dist = c(dist, riverdistance(startseg=startseg, 
                               startvert=endseg,
                               endseg=stations$seg[i], 
                               endvert = stations$vert[i],
                               rivers=streams,
                               map=TRUE))
  i = i+1
}
stations['distance'] = dist
stations = stations[order(stations$distance),]  # Order points upstream to downstream


# -------------- Smooth point data --------------
stations['BED_ELEV_m'] = stations['BED_ELEV']*0.3048 # Convert to meters


# Working dataframe
df = cbind(stations['distance'], stations['BED_ELEV_m'])
df[df==0] = NA

# Load the points that you want to add to the stream network to calculate the bed elevation
interp = pointshp2segvert(path='.',layer='points_bed_30m_spacing', rivers=streams)
dist=c()
i=1
# Loop through the points and calculate the distance along the route
while (i <= length(interp$vert)){
  dist = c(dist, riverdistance(startseg=startseg, 
                               startvert=startvert,
                               endseg=interp$seg[i], 
                               endvert = interp$vert[i],
                               rivers=streams,
                               map=TRUE))
  i = i+1
}
interp['distance'] = dist
interp = interp[order(interp$distance),]  # Order points upstream to downstream

# Use loess smoothing to interpolate the bed elevation along the stream
interp['BED_ELEV_m'] = loess_smoothing(x=df['distance'],
                                       y=df['BED_ELEV_m'],
                                       alpha=0.2,
                                       interpolate=interp$distance)

# Output the data
dataout = data.frame(interp['x_coord'],
                     interp['y_coord'],
                     interp['distance'],
                     interp['BED_ELEV_m'])

dataout = na.omit(dataout)
write.csv(dataout,'bachelor_creek_interpolated_bed_elevations.csv', row.names = FALSE)

dev.off()










