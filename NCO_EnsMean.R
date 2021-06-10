# #NCO code to average netcdfs across ESMs


variables <- c("ild_0.5C","oxygen_bottom","sst","temp_bottom","zoo_200m","zoo_50m")
for (variable in variables){
  had <- paste0('~/Dropbox/"WRAP Location^3"/2d_fields/monthly/',variable,'_monthly_roms_gfdl_1980_2100.nc')
  gfdl <- paste0('~/Dropbox/"WRAP Location^3"/2d_fields/monthly/',variable,'_monthly_roms_had_1980_2100.nc')
  ipsl <- paste0('~/Dropbox/"WRAP Location^3"/2d_fields/monthly/',variable,'_monthly_roms_ipsl_1980_2100.nc')
  file_out <- glue('~/Dropbox/"WRAP Location^3"/2d_fields/monthly/{variable}_monthly_roms_ensmean_1980_2100.nc')
  system2('nces', args=c(had,gfdl,ipsl, file_out))
}