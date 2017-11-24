require(gdalUtils)
require(raster)

#------------------------------------------------------------------------------
# Next, define functions for calculations
# Write a function to convert degrees to radians
calc_radians <- function(x) {(x * pi) / (180)}

# Write a function to calculate the Radiance from DN values
calc_radiance <- function(x,ucc){(x - 1) * ucc1}

# Write a function to calculate the TOA Reflectance from Radiance
calc_reflectance <- function(x,earth_sun_dist,irradiance,sza){
  (pi * x * (earth_sun_dist^2)) / (irradiance1 * sin(pi * sza / 180))
}

read_metadata_ASTER=function(file_name) {
  
  meta_data <- gdalinfo(file_name)
  
  lr_row <- grep('LOWERRIGHTM', meta_data)
  ul_row <- grep('UPPERLEFTM', meta_data)
  lr <- substr(meta_data[lr_row[1]], 15, 50)
  ul <- substr(meta_data[ul_row[1]], 14, 50)
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
  utm_row <- grep('UTMZONECODE', meta_data)
  utm_zone <- substr(meta_data[utm_row[1]], 1, 50)
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
  y_min_90m <- min(ul_y_90m, lr_y_90m)
  y_max_90m <- max(ul_y_90m, lr_y_90m)
  x_max_90m <- max(ul_x_90m, lr_x_90m)
  x_min_90m <- min(ul_x_90m, lr_x_90m)
  
  raster_dims_15m <- extent(x_min, x_max, y_min, y_max)
  raster_dims_30m <- extent(x_min_30m, x_max_30m, y_min_30m, y_max_30m)
  raster_dims_90m <- extent(x_min_90m, x_max_90m, y_min_90m, y_max_90m)
  
  
 
  # grab DOY from the filename and convert to day of year
  month <- substr(file_name, 12, 13)
  day <- substr(file_name, 14, 15)
  year <- substr(file_name, 16, 19)
  datep <- paste(year, month, day, sep = '-')  
  doy <- as.numeric(strftime(datep, format = '%j'))
  
  # Remove unneccessary variables
  rm(month, day, year)
  
  # need SZA--calculate by grabbing solar elevation info 
  
  sza <- meta_data[grep('SOLARDIRECTION=', meta_data)]
  clip3 <- regexpr(', ', sza) 
  sza <- as.numeric((substr(sza, (clip3 + 2), 10000)))
  
  # Need the gain designation for each band
  gain_01 <- substr(meta_data[grep('GAIN=01', meta_data)], 12, 15)
  gain_02 <- substr(meta_data[grep('GAIN=02', meta_data)], 12, 15)
  gain_04 <- substr(meta_data[grep('GAIN=04', meta_data)], 12, 15)
  gain_05 <- substr(meta_data[grep('GAIN=05', meta_data)], 12, 15)
  gain_06 <- substr(meta_data[grep('GAIN=06', meta_data)], 12, 15)
  gain_07 <- substr(meta_data[grep('GAIN=07', meta_data)], 12, 15)
  gain_08 <- substr(meta_data[grep('GAIN=08', meta_data)], 12, 15) 
  gain_09 <- substr(meta_data[grep('GAIN=09', meta_data)], 12, 15)
  gain_03b <- substr(meta_data[grep('GAIN=3B', meta_data)], 12, 15)
  gain_03n <- substr(meta_data[grep('GAIN=3N', meta_data)], 12, 15)
  
  # Calculate Earth Sun Distance (Achard and D'Souza 1994; Eva and Lambin, 1998)
  earth_sun_dist <- (1 - 0.01672 * cos(calc_radians(0.9856 * (doy - 4))))
  return(list(utm_zone=utm_zone, 
              raster_dims_15m = raster_dims_15m,
              raster_dims_30m = raster_dims_30m,
              raster_dims_90m = raster_dims_90m,
              earth_sun_dist=earth_sun_dist,
              date=as.Date(datep),
              doy=doy,
              gain_01 = gain_01,
              gain_02  = gain_02,
              gain_04  = gain_04,
              gain_05  = gain_05,
              gain_06  = gain_06,
              gain_07  = gain_07,
              gain_08  = gain_08,
              gain_09  = gain_09,
              gain_03b  = gain_03b,
              gain_03b  = gain_03b 
              ))
  
}







retrieve_aster_hdf=function(file,user='alf_crisci_1',password="la_querce_69M") 
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

readband_ASTER=function(file_name,index_k=4,nameband="b13",utm_zone,geoextent,type_data='TIR_Swath',remove=F) {
                        require(gdalUtils)
                        crs_string <- paste('+proj=utm +zone=', utm_zone, ' +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0', sep = '')
                        sds <- get_subdatasets(file_name)
                        match_tir <- grep(type_data, sds)
                        k=match_tir[index_k]
                        sub_dataset<-sds[k]
                        new_file_name <- strsplit(file_name,'.hdf')
                        tif_name <- paste(new_file_name,"_",nameband,'.tif', sep='')
                        gdal_translate(file_name, tif_name, sd_index=k, output_Raster = T)
                        aster_file <- suppressWarnings(raster(readGDAL(tif_name)))
                        proj4string(aster_file)=CRS(crs_string)
                        extent(aster_file) <- geoextent
                        aster_file[aster_file==0]=NA
                        if (remove==T) {file.remove(tif_name)}
                        return(aster_file)
}

###########################################################################################
# retrieve lst ASTER
# K1 = c1 /(lambda^5) and K2 = c2/lambda
# where c1 and c2 are first and second radiation constants from CODATA Recommended Values of the Fundamental Physical Constants: 2014

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

# Thome et al (B) is used, which uses spectral irradiance values from MODTRAN
# Ordered b1, b2, b3N, b4, b5...b9

irradiance <- c(1848,1549,1114,225.4,86.63,81.85,74.85,66.49,59.85)

klist_ASTER=list()
klist_ASTER$B10$k1 = 3.0236879000387607831e3
klist_ASTER$B10$k2 = 1.73346669879518072289e3
klist_ASTER$B11$k1 = 2.45950033786752380082e3
klist_ASTER$B11$k2 = 1.66332642774566473988e3
klist_ASTER$B12$k1 = 1.9086243591011446439e3
klist_ASTER$B12$k2 = 1.58107402197802197802e3
klist_ASTER$B13$k1 = 8.90016580863773201879e2
klist_ASTER$B13$k2 = 1.35733713207547169811e3
klist_ASTER$B14$k1 = 6.4645039694287387514e2
klist_ASTER$B14$k2 = 1.27325430088495575221e3

klist_ASTER$B10$wavelength=28.3
klist_ASTER$B11$wavelength=8.65
klist_ASTER$B12$wavelength=9.1
klist_ASTER$B13$wavelength=10.6
klist_ASTER$B14$wavelength=11.3

klist_ASTER$B10$unitconversion=6.882e-3
klist_ASTER$B11$unitconversion=6.780e-3
klist_ASTER$B12$unitconversion=6.590e-3
klist_ASTER$B13$unitconversion=5.693e-3
klist_ASTER$B14$unitconversion=5.225e-3
klist_ASTER$c2micro=14388

emiss_ASTER=list()
emiss_ASTER$mattoni_rossi$emis = 0.94
emiss_ASTER$rame_lucido$emis = 0.15
emiss_ASTER$marmo_bianco$emis = 0.95
emiss_ASTER$vernice_bianca$emis = 0.89
emiss_ASTER$vernice_nera$emis = 0.95
emiss_ASTER$acciaio_ossidato$emis = 0.94

#saveRDS(emiss_ASTER,"emiss_ASTER.rds")

table_ASTER=list(ucc=ucc,
                 irradiance=irradiance,
                 emiss_ASTER=emiss_ASTER,
                 klist_ASTER=klist_ASTER)

# Remove unneccessary variables
rm(bands,gain_names,ucc_vals,ucc,irradiance)

#saveRDS(table_ASTER,"tables_ASTER.rds")

#############################################################################################################################



retrieve_aster_hdf=function(file,user='alf_crisci_1',password="la_querce_69M") 
{download.file(url=paste0("http://e4ftl01.cr.usgs.gov/ASTT/AST_L1T.003/",
                           substr(file,16,19),".",
                           substr(file,12,13),".",
                           substr(file,14,15),
                           "/",file,'.hdf'),
                destfile=paste0(file,'.hdf'),
                method="wget",
                extra=paste0("--http-user=",user," --http-password=",password))

}

aster_date_hdf=function(file) 
{aster_date=as.Date(paste(substr(file,16,19),substr(file,12,13),substr(file,14,15),sep="-"))
            return(aster_date)
}



