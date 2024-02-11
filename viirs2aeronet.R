#to read viirs and aeronet for comparison

library("ncdf4")
library("ggplot2")
library("ggpointdensity")
library("plyr")
#library("dplyr")

date_start <- "2022010100"
date_end <- "2022060101"

type <- "_solar" 
#type <- "_lunar"
threshold <- 0.75

season <- "JFMAM"

sim <- "radius_25km"

wave0="550"
#obtype <- paste("AOD15",type,sep="")
obtype <- paste("AOD20",sep="")
simob <- paste(sim,"_",obtype,sep="")

indir <- paste("/work/noaa/gsd-fv3-dev/pagowski/R/viirs2aeronet_",
		 simob,sep="")
infile <- paste(indir,"/","viirs2aeronet",type,"_",season,".nc",sep="")


minaod_wave0 <-  1.e-3
maxaod_wave0 <-  3.0

small <- 1.e-5

cycle_interval <- 1

year_start <- substr(date_start,1,4)
month_start <- substr(date_start,5,6)
day_start <- substr(date_start,7,8)
hour_start <- substr(date_start,9,10)
minute_start <- "00"
second_start <- "00"

yyyymmdd_start <- paste(year_start,month_start,day_start,sep="-")
time_start <- paste(hour_start,minute_start,second_start,sep=":")

date_start <- as.POSIXlt(paste(yyyymmdd_start,time_start),"UTC")

year_end <- substr(date_end,1,4)
month_end <- substr(date_end,5,6)
day_end <- substr(date_end,7,8)
hour_end <- substr(date_end,9,10)
minute_end <- "00"
second_end <- "00"

yyyymmdd_end <- paste(year_end,month_end,day_end,sep="-")
time_end <- paste(hour_end,minute_end,second_end,sep=":")
date_end <- as.POSIXlt(paste(yyyymmdd_end,time_end),"UTC")

new_date <- date_start

ident <- format(new_date,"%Y%m%d%H")


nc <- nc_open(infile,readunlim=FALSE, write=FALSE )
nlocs <- nc$dim[["nlocs"]]$len
lats <- ncvar_get(varid="latitude",nc)
lons <- ncvar_get(varid="longitude",nc)
aod_aeronet <- ncvar_get(varid="aod_aeronet",nc)
aod_viirs_nn <- ncvar_get(varid="aod_viirs_nn",nc)
aod_viirs_arithmetic <- ncvar_get(varid="aod_viirs_arithmetic",nc)
aod_viirs_geometric <- ncvar_get(varid="aod_viirs_geometric",nc)
aod_viirs_random <- ncvar_get(varid="aod_viirs_random",nc)
coverage_ratio <- ncvar_get(varid="coverage_ratio",nc)
nc_close(nc)

i_obs_valid <-  which( aod_aeronet > 0. & aod_viirs_nn > 0. 
	    	       & coverage_ratio >= threshold )

lats <- lats[i_obs_valid]
lons <- lons[i_obs_valid]
aod_aeronet <- aod_aeronet[i_obs_valid]
aod_viirs_nn <- aod_viirs_nn[i_obs_valid]
aod_viirs_arithmetic <- aod_viirs_arithmetic[i_obs_valid]
aod_viirs_geometric <- aod_viirs_geometric[i_obs_valid]
aod_viirs_random <- aod_viirs_random[i_obs_valid]

aggregation <- c("NN","AR","GM","RM")

for (i in seq(1,4)) {
    agg <- aggregation[i]
    if (i==1) aod_viirs <- aod_viirs_nn	
    if (i==2) aod_viirs <- aod_viirs_arithmetic	
    if (i==3) aod_viirs <- aod_viirs_geometric	
    if (i==4) aod_viirs <- aod_viirs_random	

df <- data.frame(aod_aeronet=aod_aeronet,aod_viirs=aod_viirs)

regression <- lm(aod_viirs ~ aod_aeronet, data=df)
rsq <- as.numeric(summary(regression)$r.squared)
bias <- mean(aod_viirs-aod_aeronet)
string <- paste("nobs = ",format(length(i_obs_valid),width=5),
                " aggregation: ",agg,
       	  	" bias = ",sprintf("%+01.3f",round(bias,3)),
		" r^2 = ",format(round(rsq,3),nsmall = 3,width=5),sep="")

print(string)
coeffs <- as.numeric(coefficients(regression))
#print(coeffs)

df <- data.frame(x = aod_aeronet,y = aod_viirs,
                 d = densCols(aod_aeronet,aod_viirs,
		 nbin=4096,
             colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

library(MASS)
library(scales)

xmin <- minaod_wave0
xmax <- maxaod_wave0
ymin <- xmin
ymax <- xmax

picdir <- "./pics_viirs2aeronet/"
picname <- paste(picdir,"v2a","_",simob,"_",agg,"_",season,".png",sep="")
png(picname,width = 600, height = 600,bg="white")

p <- ggplot(df) + geom_point(aes(x, y, col = d), size = 0.5)  +
     scale_color_identity()  +

     scale_x_continuous(
                trans = 'log10',
                breaks = trans_breaks('log10', function(x) 10^x, n=3),
		labels = trans_format('log10', math_format(10^.x)),
		limits = c(xmin,xmax),expand=c(0,0)) +

     scale_y_continuous(
                trans = 'log10',
                breaks = trans_breaks('log10', function(x) 10^x, n=3),
		labels = trans_format('log10', math_format(10^.x)),
		limits = c(xmin,xmax),expand=c(0,0)) +

     theme_bw() +
     theme(plot.margin = margin(0.75, 0.75, 0.25, 0.25, "cm"),
           axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
	   axis.text.x = element_text(size=16), 
	   axis.text.y = element_text(size=16), 
	   axis.title.x = element_text(size=16), 
	   axis.title.y = element_text(size=16)) +
     geom_abline(intercept=0,slope=1,color="black",
		 linetype="solid",size=1) +
     xlab(paste("AERONET AOD ",wave0," nm",sep="")) +
     ylab(paste("VIIRS ",agg, " AOD ",wave0," nm",sep=""))

print(p)

dev.off()


}