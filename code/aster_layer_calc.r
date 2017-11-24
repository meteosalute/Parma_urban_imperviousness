#!/usr/bin/env Rscript

############################################################################################################
# Convert ASTER L1T HDF-EOS VNIR/TIR datasets from Radiance 
# and  exports  GeoTIFF files.
# Author: Alfonso Crisci 
# Contact: a.crisci@ibimet.cnr.it
# Organization: Istituto di Biometeorologia
# Date last modified: 24-11-2017

# DESCRIPTION:
# This script is used for batch processing of ASTER Images
# The script uses an ASTER L1T HDF-EOS file (.hdf) as the input.
# It based on the work of Cole Krehbiel
# https://git.earthdata.nasa.gov/projects/LPDUR/repos/aster-l1t

# USAGE: Rscript aster_layer_calc.r path/filename name_city city type


##############################################################################################################

in_dir <- args[1]
name <- args[2]
city  <- args[3]
type  <- args[4]

emis=0.95

message(paste(in_dir,name,city,type))

if (length(args) !=4) {
                    stop("At least four argument must be supplied (working dir, input file, name of output directory )\n", call.=FALSE)
} 



calc_lst=function(rbright,emis=0.92,c2micro=14388,Lwave=10.16,kelvin=F) {
                  temp=rbright / ( 1 + ( Lwave * (rbright / c2micro)) * log(emis))
                  if ( kelvin==F) { temp=temp-273.15}
                  return(temp)
} 


calc_tbright=function(r,band=13) {
                                  k=band-9
                                  ucc <- c(0.006822, 0.00678, 0.00659, 0.005693, 0.005225)
                                  k1 <- c(3024, 2460, 1909, 890, 646.4)
                                  k2 <- c(1733, 166, 1581, 1357, 1273)
                                  tir_rad <- ((r - 1) * ucc[k])
                                  sat_brit <- k2[k]/log((k1[k]/tir_rad) + 1)
                                  return(sat_brit)
                                  }

retrieve_aster_hdf=function(file,user='',password="") 
                           {download.file(url=paste0("http://e4ftl01.cr.usgs.gov/ASTT/AST_L1T.003/",
                           substr(file,16,19),".",
                           substr(file,12,13),".",
                           substr(file,14,15),"/",file,'.hdf'),
                           destfile=paste0(file,'.hdf'),
                           method="wget",
                           extra=paste0("--http-user=",user," --http-password=",password))
                           }

aster_date_hdf=function(file)  {aster_date=as.Date(paste(substr(file,16,19),substr(file,12,13),substr(file,14,15),sep="-"))
                               return(aster_date)
                               }

#############################################################################################################
# Convert ASTER L1T HDF-EOS VNIR/TIR datasets from Radiance 
# and  exports as GeoTIFF files.
#-------------------------------------------------------------------------------
# Author: Cole Krehbiel
# Contact: LPDAAC@usgs.gov  
# Organization: Land Processes Distributed Active Archive Center
# Date last modified: 03-06-2017
#-------------------------------------------------------------------------------
# DESCRIPTION:
# This script demonstrates how to convert ASTER L1T data from Digital Number(DN)
# to radiance (w/m2/sr/�m), and from radiance into Top of Atmosphere (TOA) 
# reflectance. 

# The script uses an ASTER L1T HDF-EOS file (.hdf) as the input and outputs 
# georeferenced tagged image file format (GeoTIFF) files for each of the VNIR 
# and SWIR science datasets contained in the original ASTER L1T file.
##############################################################################################################

require(rgdal)
require(raster)
require(gdalUtils)
require(tiff)
require(rgdal)


args = commandArgs(trailingOnly=TRUE)



#############################################################################################################
# Set up calculations
# 1. DN to Radiance (Abrams, 1999)
# Radiance = (DN-1)* Unit Conversion Coefficient

# 2. Radiance to TOA Reflectance
# Reflectance_TOA = (pi*Lrad*d2)/(esuni*COS(z))

# Define the following:
# Unit Conversion Coefficient = ucc
# pi = pi
# Radiance,Lrad  = rad
# esuni = esun 
# z = solare

# Order for ucc (Abrams, 1999) is: Band 1 high, normal, low; Band 2 h, n, l; 
# b3 h, n, l (3N & 3B the same) 
# Construct a dataframe for the UCC values:

bands <- c('1', '2', '3N', '3B', '4', '5', '6', '7', '8', '9')
gain_names <- c('Band', 'High Gain', 'Normal', 'Low Gain 1', 'Low Gain 2')
ucc_vals <- matrix( c(0.676, 1.688, 2.25, 0, 0.708, 1.415, 1.89, 0, 0.423, 
                      0.862, 1.15, 0, 0.423, 0.862, 1.15, 0, 0.1087, 0.2174, 
                      0.2900, 0.2900, 0.0348, 0.0696, 0.0925, 0.4090, 0.0313,
                      0.0625, 0.0830, 0.3900, 0.0299, 0.0597, 0.0795, 0.3320,
                      0.0209,0.0417, 0.0556, 0.2450, 0.0159, 0.0318, 0.0424, 
                      0.2650), nrow = 10, ncol = 4, byrow = TRUE)
ucc <- data.frame( bands, ucc_vals )

names(ucc) <- gain_names

# Remove unneccessary variables
rm(bands,gain_names,ucc_vals)

# Thome et al (B) is used, which uses spectral irradiance values from MODTRAN
# Ordered b1, b2, b3N, b4, b5...b9

irradiance <- c(1848,1549,1114,225.4,86.63,81.85,74.85,66.49,59.85)

#############################################################################################################
# Next, define functions for calculations
# Write a function to convert degrees to radians

calc_radians <- function(x) {(x * pi) / (180)}
 
# Write a function to calculate the Radiance from DN values
calc_radiance <- function(x){(x - 1) * ucc1}

# Write a function to calculate the TOA Reflectance from Radiance

calc_reflectance <- function(x){
                               (pi * x * (earth_sun_dist^2)) / (irradiance1 * sin(pi * sza / 180))
                               }







old_dir=getwd()

#############################################################################################################

setwd(in_dir)

# Maintains the original filename
file_name <- paste0(name,".hdf")
 
# if (!file.exists(file_list)) {retrieve_aster_hdf(file_list)}}

# Create and set output directory


out_dir <- paste(in_dir,'/',city,"_",as.character(aster_date_hdf(file_name)),'/', sep='')

suppressWarnings(dir.create(out_dir))

###################################################################
# grab DOY from the filename and convert to day of year

  month <- substr(file_name, 12, 13)
  day <- substr(file_name, 14, 15)
  year <- substr(file_name, 16, 19)
  date <- paste(year, month, day, sep = '-')  
  doy <- as.numeric(strftime(date, format = '%j'))
  
# Remove unneccessary variables

  rm(month, day, year, date)
  
  # need SZA--calculate by grabbing solar elevation info 

  md  <- gdalinfo(file_name)
  sza <- md[grep('SOLARDIRECTION=', md)]
  clip3 <- regexpr(', ', sza) 
  sza <- as.numeric((substr(sza, (clip3 + 2), 10000)))
  
  # Need the gain designation for each band
  gain_01 <- gsub(' ', '', strsplit(md[grep('GAIN.*=01', md)], ',')[[1]][[2]])
  gain_02 <- gsub(' ', '', strsplit(md[grep('GAIN.*=02', md)], ',')[[1]][[2]])
  gain_04 <- gsub(' ', '', strsplit(md[grep('GAIN.*=04', md)], ',')[[1]][[2]])
  gain_05 <- gsub(' ', '', strsplit(md[grep('GAIN.*=05', md)], ',')[[1]][[2]])
  gain_06 <- gsub(' ', '', strsplit(md[grep('GAIN.*=06', md)], ',')[[1]][[2]])
  gain_07 <- gsub(' ', '', strsplit(md[grep('GAIN.*=07', md)], ',')[[1]][[2]])
  gain_08 <- gsub(' ', '', strsplit(md[grep('GAIN.*=08', md)], ',')[[1]][[2]])
  gain_09 <- gsub(' ', '', strsplit(md[grep('GAIN.*=09', md)], ',')[[1]][[2]])
  gain_03b <- gsub(' ', '', strsplit(md[grep('GAIN.*=3B', md)], ',')[[1]][[2]])
  gain_03n <- gsub(' ', '', strsplit(md[grep('GAIN.*=3N', md)], ',')[[1]][[2]])
  
  # Calculate Earth Sun Distance (Achard and D'Souza 1994; Eva and Lambin, 1998)
  earth_sun_dist <- (1 - 0.01672 * cos(calc_radians(0.9856 * (doy - 4))))
  #-----------------------------------------------------------------------------
  # Define CRS
  # Define Upper left and lower right--need for x, y min/max
  # For offset (pixel size / 2), needs to be defined for ASTER pixel 
  # resolutions (15, 30)
  
  # Grab LR and UL values
  lr <- substr(md[grep('LOWERRIGHTM', md)], 15, 50)
  ul <- substr(md[grep('UPPERLEFTM', md)], 14, 50)
  clip4 <- regexpr(', ' , ul) 
  clip5 <- regexpr(', ', lr) 
  
  # Define LR and UL x and y values for 15m VNIR Data
  ul_y <- as.numeric((substr(ul, 1, (clip4 - 1)))) + 7.5
  ul_x <- as.numeric((substr(ul, (clip4 + 2), 10000))) - 7.5
  lr_y <- as.numeric((substr(lr, 1, (clip5 - 1)))) - 7.5
  lr_x <- as.numeric((substr(lr, (clip5 + 2) , 10000))) + 7.5
  
  # Define LR and UL x and y values for 30m SWIR Data
  ul_y_30m <- as.numeric((substr(ul, 1, (clip4 - 1)))) + 15
  ul_x_30m <- as.numeric((substr(ul, (clip4 + 2), 10000))) - 15
  lr_y_30m <- as.numeric((substr(lr, 1, (clip5 - 1)))) - 15
  lr_x_30m <- as.numeric((substr(lr, (clip5 + 2) , 10000))) + 15
  
  # Define LR and UL x and y values for 90m TIR Data
  ul_y_90m <- as.numeric((substr(ul, 1, (clip4 - 1)))) + 45
  ul_x_90m <- as.numeric((substr(ul, (clip4 + 2), 10000))) - 45
  lr_y_90m <- as.numeric((substr(lr, 1, (clip5 - 1)))) - 45
  lr_x_90m <- as.numeric((substr(lr, (clip5 + 2) , 10000))) + 45

  # Define UTM zone
  utm_row <- grep('UTMZONECODE', md)
  utm_zone <- substr(md[utm_row[1]], 1, 50)
  clip6 <- regexpr('=', utm_zone) 
  utm_zone <- substr(utm_zone, clip6 + 1, 50)
  
  # Configure extent properties (15m VNIR)
  y_min <- min(ul_y, lr_y)
  y_max <- max(ul_y, lr_y)
  x_max <- max(ul_x, lr_x)
  x_min <- min(ul_x, lr_x)
  
  # Configure extent properties (30m SWIR)

  y_min_30m <- min(ul_y_30m, lr_y_30m)
  y_max_30m <- max(ul_y_30m, lr_y_30m)
  x_max_30m <- max(ul_x_30m, lr_x_30m)
  x_min_30m <- min(ul_x_30m, lr_x_30m)

  # Configure extent properties (90m TIR)

  y_min_90m <- min(ul_y_90m, lr_y_90m); 
  y_max_90m <- max(ul_y_90m, lr_y_90m)
  x_max_90m <- max(ul_x_90m, lr_x_90m); 
  x_min_90m <- min(ul_x_90m, lr_x_90m)
  
 
  raster_dims_15m <- extent(x_min, x_max, y_min, y_max)
  raster_dims_30m <- extent(x_min_30m, x_max_30m, y_min_30m, y_max_30m)
  raster_dims_90m <- extent(x_min_90m, x_max_90m, y_min_90m, y_max_90m)
  
  # Compile Cordinate Reference System string to attach projection information
  crs_string <- paste('+proj=utm +zone=', utm_zone, ' +datum=WGS84 +units=m 
                      +no_defs +ellps=WGS84 +towgs84=0,0,0', sep = '')
  
  # Remove unneccessary variables
  rm(clip4, clip5, clip6, lr, lr_x, lr_y, md, ul, ul_x, ul_y, utm_zone,
     x_min, x_max, y_min, y_max, utm_row)
  
  # Get a list of sds names
  sds <- get_subdatasets(file_name)
  
  # Limit loop to SDS that contain VNIR/SWIR data (9 max)
 
  match_vnir <- grep('VNIR_Swath', sds)
  match_tir <- grep('TIR_Swath', sds)   
 
###################################################################################################################à

  if (length(match_vnir) > 1){
    for (k in min(match_vnir):max(match_vnir)) {
      
     # Isolate the name of the first sds
      sub_dataset <- sds[k]
      
      # Get the name of the specific SDS
      clip2 <- max(unlist((gregexpr(':', sub_dataset)))) 
      
      # Generate output name for tif
      new_file_name <- strsplit(file_name, '.hdf')
      tif_name <- paste(out_dir, new_file_name, '_', substr(sub_dataset, 
                       (clip2 + 1), 10000),'.tif', sep='')
      sd_name <- paste(new_file_name, substr(sub_dataset, (clip2 + 1), 10000), 
                       sep = '_')
      sub_name <- paste(new_file_name, 'ImageData', sep = '_')
      ast_band_name <- gsub(sub_name, '', sd_name)
    
      # Extract specified SDS and export as Geotiff
      
      gdal_translate(file_name, tif_name, sd_index=k)
      
      aster_file <- raster(tif_name, crs = crs_string)
      
      if (ast_band_name == '1'){
        # Need to know which gain value you use
        if (gain_01 == 'HGH'){
          ucc1 <- ucc[1, 2] 
        }else if(gain_01 == 'NOR'){
          ucc1 <- ucc[1, 3] 
        }else{
          ucc1 <- ucc[1, 4] 
        }
        irradiance1 <- irradiance[1]
        
        # Define Extent
        extent(aster_file) <- raster_dims_15m
    
      }else if (ast_band_name == '2'){
        # Need to know which gain value you use
        if (gain_02 == 'HGH'){
          ucc1 <- ucc[2, 2] 
        }else if(gain_02 == 'NOR'){
          ucc1 <- ucc[2, 3] 
        }else{
          ucc1 <- ucc[2, 4] 
        }
        irradiance1 <- irradiance[2]
        
        # Define Extent
        extent(aster_file) <- raster_dims_15m
    
      } else if (ast_band_name == '3N'){
        # Need to know which gain value you use
        if (gain_03n == 'HGH'){
          ucc1 <- ucc[3, 2] 
        } else if(gain_03n == 'NOR'){
          ucc1 <- ucc[3, 3] 
        } else{
          ucc1 <- ucc[3, 4] 
        }
        irradiance1 <- irradiance[3]
        
        # Define Extent
        extent(aster_file) <- raster_dims_15m

      }


      # Set up output file names
      ref_out_name <- gsub(paste(ast_band_name, '.tif', sep = ''), 
                           paste(ast_band_name, '_reflectance.tif', sep = ''), 
                           tif_name)
      rad_out_name <- gsub(paste(ast_band_name, '.tif', sep = ''), 
                           paste(ast_band_name, '_radiance.tif', sep = ''),
                           tif_name)
      # Convert DN to large raster layer
      aster_file <- calc(aster_file, fun =function(x){x})
      
      # Export the DN raster layer file (Geotiff format) to the output directory
      writeRaster(aster_file, filename = tif_name,  options = 'INTERLEAVE=BAND',
                  NAflag = 0, format = 'GTiff', datatype = 'INT1U', 
                  overwrite = TRUE)
      
      # Convert from DN to Radiance
      rad <- calc(aster_file, calc_radiance)
      rad[rad == calc_radiance(0)] <- 0
      rm(aster_file)
      # export the raster layer file (Geotiff format) to the output directory
      writeRaster(rad, filename = rad_out_name,  options = 'INTERLEAVE=BAND', 
                  NAflag = 0, format = 'GTiff', datatype = 'FLT8S', 
                  overwrite = TRUE)
      
      # Convert from Radiance to TOA Reflectance
      ref <- calc(rad, calc_reflectance)
      rm(rad)
      # export the raster layer file (Geotiff format) to the output directory
      
      writeRaster(ref, filename = ref_out_name,  options = 'INTERLEAVE=BAND', 
                  NAflag = 0, format = 'GTiff', datatype = 'FLT8S', 
                  overwrite = TRUE)
      
       # Remove unneccessary variables
 
       rm(ucc1, irradiance1, ref_out_name, rad_out_name, ref, sub_dataset,sd_name, sub_name, tif_name, new_file_name)
    
    }
  }
  
###################################################################################################################à
  
   if (length(match_tir) > 0) {
   
   for (k in min(match_tir):max(match_tir)){
   
   # Isolate the name of the first sds
    sub_dataset<- sds[k]
    
    # Get the name of the specific SDS
    clip2 <- max(unlist((gregexpr(':', sub_dataset)))) 
    
    # Generate output name for tif
    new_file_name <- strsplit(file_name, '.hdf')
    tif_name <- paste(out_dir, new_file_name, '_', substr(sub_dataset, 
                     (clip2 + 1), 10000),'.tif', sep='')
    sd_name <- paste(new_file_name, substr(sub_dataset, (clip2 + 1), 10000),
                     sep = '_')
    sub_name <- paste(new_file_name, 'ImageData', sep = '_')
    ast_band_name <- gsub(sub_name, '', sd_name)
    
    # Extract specified SDS and export as Geotiff
    
    gdal_translate(file_name, tif_name, sd_index=k, output_Raster = FALSE)
    
    # Open geotiff and add projection (CRS)
    
    aster_file <- suppressWarnings(raster(readGDAL(tif_name,as.is=T)))
    proj4string(aster_file)=crs_string
    extent(aster_file) <- raster_dims_90m

    # Convert to large raster layer
    aster_file <- calc(aster_file, fun =function(x){x} )
    
    # Export the raster layer file (Geotiff format) to the output directory
    writeRaster(aster_file, filename = tif_name,  options = 'INTERLEAVE=BAND',
                format = 'GTiff', datatype = 'INT2U', overwrite = TRUE, 
                NAflag = 0)
    # Remove unneccessary variables
    
    
    rm(aster_file, sub_dataset, sd_name, sub_name, tif_name, new_file_name)
    }
   } 
  
###################################################################################################################à
# Remove unneccessary variables
  
  
  rm( earth_sun_dist, doy,gain_01,gain_02, gain_03b, 
    gain_03n, gain_04, gain_05, gain_06, gain_07, gain_08, gain_09, sza,
    ast_band_name, clip2, crs_string, file_name,sds, lr_x_30m, lr_x_90m, 
     lr_y_30m, lr_y_90m,   raster_dims_15m, 
     raster_dims_30m, raster_dims_90m, ul_x_30m, ul_x_90m, ul_y_30m, ul_y_90m,
     x_max_30m, x_max_90m, x_min_30m, x_min_90m, y_max_30m, y_max_90m,y_min_30m,
     y_min_90m)

###################################################################################################################à
res=list()

if (length(match_vnir) > 0 ) {
    blue_radiance=raster(paste0(out_dir,name,"_","ImageData1_radiance.tif"))
    blue_reflectance=raster(paste0(out_dir,name,"_","ImageData1_reflectance.tif"))
    green_radiance=raster(paste0(out_dir,name,"_","ImageData2_radiance.tif"))
    green_reflectance=raster(paste0(out_dir,name,"_","ImageData2_reflectance.tif"))
    red_radiance=raster(paste0(out_dir,name,"_","ImageData3N_radiance.tif"))
    red_reflectance=raster(paste0(out_dir,name,"_","ImageData3N_reflectance.tif"))
    ASTERVNIR=stack(blue_radiance,blue_reflectance, green_radiance,green_reflectance,red_radiance,red_reflectance)
   # writeRaster(ASTERVNIR, filename=paste0(out_dir,name,"_stackVNIR.tif"), overwrite=TRUE)
}

if (length(match_tir) > 0 ) {
  
  res=list()
  
  res$b13_tbright=calc_tbright(suppressWarnings(raster(readGDAL(paste0(out_dir,name,"_","ImageData13.tif"),as.is=T))),band=13)
  res$b14_tbright=calc_tbright(suppressWarnings(raster(readGDAL(paste0(out_dir,name,"_","ImageData14.tif"),as.is=T))),band=14)
  
  res$b13_LST=calc_lst(res$b13_tbright,emis=emis,Lwave=10.6)
  res$b14_LST=calc_lst(res$b14_tbright,emis=emis,Lwave=11.3)
  
  writeRaster(res$b13_LST, filename=paste0(out_dir,name,"_",type,"_b13_LST.tif"), overwrite=TRUE)
  writeRaster(res$b14_LST, filename=paste0(out_dir,name,"_",type,"_b14_LST.tif"), overwrite=TRUE)
  writeRaster(res$b13_tbright, filename=paste0(out_dir,name,"_",type,"_b13_TB.tif"), overwrite=TRUE)
  writeRaster(res$b14_tbright, filename=paste0(out_dir,name,"_",type,"_b14_TB.tif"), overwrite=TRUE)
}


setwd(out_dir)
olds=list.files(pattern=name, full.names = F)
news=gsub(name,paste0(city,"_",as.character(aster_date_hdf(name)),"_",type),olds)
news=gsub(paste0("_",type,"_",type),news)
file.rename(olds,news)
saveRDS(res,file = paste0(city,"_",as.character(aster_date_hdf(file_name)),"_",type,'.rds'))

setwd(old_dir)


###################################################################################################################à
# References
# ABRAMS, M., HOOK, S., and RAMACHANDRAN, B., 1999, Aster user handbook,
# Version 2, NASA/Jet Propulsion Laboratory, Pasadena, CA, at 
# https://asterweb.jpl.nasa.gov/content/03_data/04_Documents/
# aster_user_guide_v2.pdf

# ARCHARD, F., AND D'SOUZA, G., 1994, Collection and pre-processing of 
# NOAA-AVHRR 1km resolution data for tropical forest resource assessment. 
# Report EUR 16055, European Commission, Luxembourg, at 
# http://bookshop.europa.eu/en/collection-and-pre-processing-of-noaa-avhrr-1-km
# -resolution-data-for-tropical-forest-resource-assessment-pbCLNA16055/

# EVA, H., AND LAMBIN, E.F., 1998, Burnt area mapping in Central Africa using 
# ATSR data, International Journal of Remote Sensing, v. 19, no. 18, 3473-3497, 
# at http://dx.doi.org/10.1080/014311698213768

# Thome, K.J., Biggar, S.F., and Slater, P.N., 2001, Effects of assumed solar
# spectral irradiance on intercomparisons of earth-observing sensors. In 
# International Symposium on Remote Sensing, International Society for Optics
# and Photonics, pp. 260-269, at http://dx.doi.org/10.1117/12.450668.
###################################################################################################################à
