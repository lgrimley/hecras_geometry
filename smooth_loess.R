#Smooth channel bottom elevations

setwd("C:/Users/lelise/Desktop/Bathy_Smoothing/")

NR_rec <- read.csv("NR_rec.csv", header=TRUE)
data = read.csv("NR_BE_clean_NewRiver.csv", header=TRUE)

# Sort data by stream station
NR_rec = NR_rec[order(NR_rec$STREAM_STN),]
data = data[order(data$STREAM_STN),]


# Plot channel bed
plot(data$BED_ELEV_m, x=data$STREAM_STN, type="l",lwd=3,col='black',
     xlab="Station", ylab="Bed Elevation_m")
points(NR_rec$BED_ELEV_m,x=NR_rec$STREAM_STN, col='red')

plot(NR_rec$BED_ELEV_m, x=NR_rec$STREAM_STN, type="l",lwd=3,col='black',
     xlab="Station", ylab="Bed Elevation_m")



# BED smoothed output
smoothed10 <- predict(loess(bed_m ~ FID, data=NR_rec, span=0.10)) 
smoothed15 <- predict(loess(bed_m ~ FID, data=NR_rec, span=0.15)) 
smoothed20 <- predict(loess(bed_m ~ FID, data=NR_rec, span=0.20)) 
smoothed25 <- predict(loess(bed_m ~ FID, data=NR_rec, span=0.25)) 

plot(NR_rec$bed_m, x=NR_rec$STREAM_STN, type="l",lwd=3,col='darkgray',
     #xlim=c(0,1500), ylim=c(-4.5,6.5),
     main="Loess Smoothing and Prediction",
     xlab="Distance from Upstream", ylab="Bed Elevation_m")
lines(NR2$BED_ELEV_m, x=NR2$STREAM_STN, col='red', lwd=2)

lines(smoothed10, x=NR_rec$FID, col="red", lwd=2)
lines(smoothed15, x=NR_rec$FID, col="purple", lwd=2)
lines(smoothed20, x=NR_rec$FID, col="green", lwd=2)
lines(smoothed25, x=NR_rec$FID, col="blue", lwd=2)


#BANK smoothed output
smoothed10 <- predict(loess(bank_m ~ FID, data=NR_rec, span=0.10)) 
smoothed15 <- predict(loess(bank_m ~ FID, data=NR_rec, span=0.15)) 
smoothed20 <- predict(loess(bank_m ~ FID, data=NR_rec, span=0.20)) 
smoothed25 <- predict(loess(bank_m ~ FID, data=NR_rec, span=0.25)) 

plot(NR_rec$bank_m, x=NR_rec$FID, type="l", main="Loess Smoothing and Prediction", xlab="Distance from Upstream", ylab="Bank Elevation (m)")
abline(v=616, col="black", lty=2)
abline(v=1418, col="black", lty=2)
lines(smoothed10, x=NR_rec$FID, col="red")
lines(smoothed15, x=NR_rec$FID, col="purple")
lines(smoothed20, x=NR_rec$FID, col="green")
lines(smoothed25, x=NR_rec$FID, col="blue")

#WIDTH smoothed output
smoothed10 <- predict(loess(width_m ~ FID, data=NR_rec, span=0.10)) 
smoothed15 <- predict(loess(width_m ~ FID, data=NR_rec, span=0.15)) 
smoothed20 <- predict(loess(width_m ~ FID, data=NR_rec, span=0.20)) 
smoothed25 <- predict(loess(width_m ~ FID, data=NR_rec, span=0.25)) 

plot(NR_rec$width_m, x=NR_rec$FID, type="l", main="Loess Smoothing and Prediction", xlab="Distance from Upstream", ylab="Bank Elevation (m)")
abline(v=616, col="black", lty=2)
abline(v=1418, col="black", lty=2)
lines(smoothed10, x=NR_rec$FID, col="red")
lines(smoothed15, x=NR_rec$FID, col="purple")
lines(smoothed20, x=NR_rec$FID, col="green")
lines(smoothed25, x=NR_rec$FID, col="blue")

#ADD SMOOTHED DATA TO TABLE
NR_rec$bed_m_smooth <- smoothed15 <- predict(loess(bed_m ~ FID, data=NR_rec, span=0.15))
NR_rec$bank_m_smooth <- smoothed15 <- predict(loess(bank_m ~ FID, data=NR_rec, span=0.15))
NR_rec$width_m_smooth <- smoothed15 <- predict(loess(width_m ~ FID, data=NR_rec, span=0.15))

plot(NR_rec$width_m_smooth)

NR_rec2 <- NR_rec[c("X","Y","bed_m","bank_m","width_m","bed_m_smooth","bank_m_smooth","width_m_smooth")]
write.table(NR_rec, "NR_rec2.csv", append = FALSE, sep = ",")
