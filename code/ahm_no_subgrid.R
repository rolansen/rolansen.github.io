library(terra)
library(sf)

#TODO: write logs with timestamps
my_options <- options(digits.secs = 3)  

INTERP_GRID_DISTANCE_FROM_SAMPLE_TOLERANCE <- 1 #in units of crs
TERRA_DESCRIBE_DATA_TYPE_REGEX <- '(?<=Type=).+?(?=[ ,])' #Perl-compatible. we use a "positive lookbehind" to find all text after 'Type=', then .+? to lazily grab any character . at least once, then the "positive lookahead" (?=[ ,]) to look for anything which occurs before a ' ' or a ','
GDAL_DATA_TYPE_TO_L_MATRIX <- matrix(c('Byte',2^8, 
                                       'UInt16',2^16,
                                       'Int16',2^16,
                                       'UInt32',2^32,
                                       'Int32',2^32,
                                       'Float32',2^32,
                                       'Float64',2^64), 
                                     dimnames = list(c(), c("Data type", "Depth")),
                                     byrow=TRUE, ncol=2) #1st col is GDAL data type, 2nd # possible values. complex types not listed
CDF_GRID_CENTER_BUFFER_SIZE <- 0.1 #should be small

get_raster_icdf <- function(in_raster, added_noise=0) {
  #largely ripped off from RSToolbox: https://github.com/bleutner/RStoolbox/blob/master/R/histMatch.R
  img_cdf <- get_raster_cdf(in_raster, added_noise) #ecdf(as.vector(in_raster))
  img_icdf_kn <- knots(img_cdf) #unique cell values from reference image used to make CDF, ie x axis on CDF
  img_icdf_y <- img_cdf(img_icdf_kn) #probability for given cell value from reference image, ie y axis on CDF
  img_icdf_limits <- c(quantile(img_cdf, 0), quantile(img_cdf, 1)) #min/max of ref cdf
  img_icdf <- approxfun(img_icdf_y, img_icdf_kn, method = "linear",
                        yleft = img_icdf_limits[1] , yright = img_icdf_limits[2], ties = "ordered")
  return(img_icdf)
}

get_raster_cdf <- function(in_raster, added_noise=0) {
  in_raster_as_vector <- as.vector(in_raster)
  if (added_noise != 0) {
    noise_vect <- runif(length(in_raster_as_vector), min=-abs(added_noise), max=-abs(added_noise))
    in_raster_as_vector <- in_raster_as_vector + noise_vect
  }
  return(ecdf(in_raster_as_vector))
}

get_data_type_from_desc <- function(in_desc_text) {
  #input is output from terra::describe
  return(regmatches(in_desc_text, regexpr(TERRA_DESCRIBE_DATA_TYPE_REGEX, in_desc_text, perl=TRUE))[1])
}

get_L_from_image_path <- function(in_raster_path) {
  #L= bit depth. we assume all bands have the same bit depth as the first band
  img_desc <- describe(in_raster_path)
  img_dt <- get_data_type_from_desc(img_desc)
  img_L <- as.numeric(GDAL_DATA_TYPE_TO_L_MATRIX[GDAL_DATA_TYPE_TO_L_MATRIX[,"Data type"]==img_dt, "Depth"])
  return(img_L)
}

ahm <- function(aoi_poly, aerial_path, sat_paths, 
                c_region_size, grid_lindim_length, out_name,
                added_aerial_noise=0, added_sat_noise=0) {
  #aoi_poly is some sf object with one polygon record. should be within aerial and sat
  #aerial_path is image with first 3 bands being red, green, blue
  #sat_paths is e.g. c(red_path, green_path, blue_path)
  #TODO: what if aerials are multiple files? sat in 1?
  
  print('Making directory for outputs')
  scratch_dir_path <- file.path(dirname(out_name), 'scratch') 
  dir.create(scratch_dir_path, showWarnings = FALSE, recursive = TRUE)
  
  #read images
  print('Reading images')
  aerial_not_cropped_to_aoi <- rast(aerial_path)
  sat_not_cropped_to_aoi <- rast(sat_paths) 
  
  #reproject aoi_poly to aerial & sat crs's
  aoi_poly_aerial <- st_transform(aoi_poly, crs(aerial_not_cropped_to_aoi))
  aoi_poly_sat <- st_transform(aoi_poly, crs(sat_not_cropped_to_aoi))
  
  #crop images
  print('Cropping images to AOI')
  aerial <- crop(aerial_not_cropped_to_aoi, vect(aoi_poly_aerial))
  sat <- crop(sat_not_cropped_to_aoi, vect(aoi_poly_sat))
  
  #ideally we'd add noise to images here if requested
  #I tried this on 9/18/22 tho and crashed my machine
  #the method used seemed to basically be recommended by Hijmans on terra's r-spatial site,
  #although his image was much smaller
  #so instead we do this when finding CDF's/iCDF's, for now
  
  #get data types for aerial and sat  
  #note we assume all sat files have same bit depth as the first file
  aerial_L <- get_L_from_image_path(aerial_path)
  sat_L <- get_L_from_image_path(sat_paths[1])
  sat_dt <- get_data_type_from_desc(describe(sat_paths[1])) #used for writing output images
  
  #define cdf grid
  print('Creating CDF/iCDF grid')
  cdf_grid <- st_as_sf(st_make_grid(aoi_poly_aerial, cellsize=grid_lindim_length))
  cdf_grid <- cdf_grid[aoi_poly_aerial,]
  cdf_grid$cdf_id <- 1:nrow(cdf_grid)
  cdf_grid_centers <- st_centroid(cdf_grid)
  
  #half grid for border/corner points
  print('Creating interpolation grid')
  half_grid_lindim_length <- grid_lindim_length/2
  half_grid <- st_as_sf(st_make_grid(cdf_grid, cellsize=half_grid_lindim_length))
  half_grid <- st_intersection(half_grid, st_union(cdf_grid))
  half_grid <- half_grid[st_geometry_type(half_grid)=='POLYGON',]
  half_grid$n_neighbors <- lengths(st_intersects(half_grid)) - 1
  half_grid$interp_type <- 'ordinary'
  half_grid[half_grid$n_neighbors == 7,'interp_type'] <- 'hinge' #TODO: check if this is right
  half_grid[half_grid$n_neighbors < 7,'interp_type'] <- 'border' #need >= 4 points for bilinear interpolation
  half_grid[half_grid$n_neighbors <= 3,'interp_type'] <- 'corner'
  
  #interp grid for ordinary cells
  ordinary_cells_from_half_grid <- half_grid[half_grid$interp_type=='ordinary',]
  interp_grid_ordinary <- st_as_sf(st_make_grid(ordinary_cells_from_half_grid, 
                                                cellsize=grid_lindim_length))
  interp_grid_ordinary <- st_intersection(interp_grid_ordinary, st_union(ordinary_cells_from_half_grid))
  interp_grid_ordinary <- interp_grid_ordinary[st_geometry_type(interp_grid_ordinary)=='POLYGON',]
  interp_grid_ordinary$cdf_neighbors <- st_is_within_distance(interp_grid_ordinary,
                                                              cdf_grid_centers,
                                                              dist=INTERP_GRID_DISTANCE_FROM_SAMPLE_TOLERANCE)
  interp_grid_ordinary$interp_type <- 'ordinary'
  
  #interp grid for corner & border cells
  #TODO?: not any practical compute time cost, but we are a little inefficient when
    #we find cdf_neighbors again for agg cells. should we do something about this?
  unusual_cells_from_half_grid <- half_grid[half_grid$interp_type!='ordinary','interp_type']
  unusual_cells_from_half_grid_cdf_neighbors <- st_is_within_distance(
    unusual_cells_from_half_grid, cdf_grid_centers, 
    dist= half_grid_lindim_length + INTERP_GRID_DISTANCE_FROM_SAMPLE_TOLERANCE)
  unusual_cdf_neighbors_as_character <- as.character(
    unusual_cells_from_half_grid_cdf_neighbors[1:nrow(unusual_cells_from_half_grid)])
  unusual_cells_from_half_grid_agg <- aggregate(unusual_cells_from_half_grid, 
                                            by = list(unusual_cdf_neighbors_as_character,
                                                      unusual_cells_from_half_grid$interp_type),
                                            FUN=min)[,c('interp_type')]
  unusual_cells_from_half_grid_agg$cdf_neighbors <- st_is_within_distance(
    unusual_cells_from_half_grid_agg, cdf_grid_centers, 
    dist= half_grid_lindim_length + INTERP_GRID_DISTANCE_FROM_SAMPLE_TOLERANCE)
  st_geometry(unusual_cells_from_half_grid_agg) <- 'x'
  
  #combine interp grids for corner and border cells
  interp_grid <- rbind(interp_grid_ordinary, unusual_cells_from_half_grid_agg)
  interp_grid$interp_id <- 1:nrow(interp_grid)
  
  #define matrix which will hold objects for hist matching:
  #sat iCDF's, aerial CDF's, and n pixels in each cdf_grid cell
  print('Finding aerial CDF\'s and satellite iCDF\'s')
  hist_match_params_mat <- matrix(ncol=9, dimnames = list(c(), 
                                                     c("cdf_id",  
                                                       "red_sat_icdf", "green_sat_icdf", "blue_sat_icdf",
                                                       "ncell_sat",
                                                       "red_aerial_cdf", "green_aerial_cdf", "blue_aerial_cdf",
                                                       "ncell_aerial")))[c(-1),]
  
  #find sat iCDF's, aerial CDF's, and number of cells of each image for each cell of cdf_grid
  for (row in 1:nrow(cdf_grid)) {
    print(paste('CDF/iCDF', row))
    cdf_cell <- cdf_grid[row,]
    
    #find computational region
    if (c_region_size != grid_lindim_length) {
      half_c_region_size <- c_region_size/2
      extent_center <- cdf_grid_centers[cdf_grid_centers$cdf_id==cdf_cell$cdf_id,]
      center_coords <- st_coordinates(extent_center)
      center_x <- center_coords[1]
      center_y <- center_coords[2]
      crop_xmin <- center_x - half_c_region_size
      crop_xmax <- center_x + half_c_region_size
      crop_ymin <- center_y - half_c_region_size
      crop_ymax <- center_y + half_c_region_size
      crop_polygon <- st_polygon(list(matrix(c(crop_xmin,crop_ymin, crop_xmax,crop_ymin,
                                               crop_xmax,crop_ymax, crop_xmin,crop_ymax,
                                               crop_xmin,crop_ymin), byrow=TRUE, ncol=2)))
      crop_polygon_sat <- st_transform(crop_polygon, crs(sat))
      
    } else {
      crop_polygon <- cdf_cell
      crop_polygon_sat <- st_transform(cdf_cell, crs(sat))
    }
    
    #crop images to computational region
    aerial_cdf_cropped <- crop(aerial, crop_polygon)
    sat_cdf_cropped <- crop(sat, crop_polygon_sat)
    
    #find LS iCDF's & ncells
    sat_red_icdf <- get_raster_icdf(sat_cdf_cropped[[1]], added_sat_noise)
    sat_green_icdf <- get_raster_icdf(sat_cdf_cropped[[2]], added_sat_noise)
    sat_blue_icdf <- get_raster_icdf(sat_cdf_cropped[[3]], added_sat_noise)
    ncell_sat <- ncell(sat_cdf_cropped)
    
    #find aerial cdf's & ncells
    aerial_red_cdf <- get_raster_cdf(aerial_cdf_cropped[[1]], added_aerial_noise) #ecdf(as.vector(aerial_cdf_cropped[[1]]))
    aerial_green_cdf <- get_raster_cdf(aerial_cdf_cropped[[2]], added_aerial_noise) #ecdf(as.vector(aerial_cdf_cropped[[2]]))
    aerial_blue_cdf <- get_raster_cdf(aerial_cdf_cropped[[3]], added_aerial_noise) #ecdf(as.vector(aerial_cdf_cropped[[3]]))
    ncell_aerial <- ncell(aerial_cdf_cropped)
    
    #append to matrix
    masked_list <- list(cdf_cell$cdf_id, sat_red_icdf, sat_green_icdf, sat_blue_icdf, ncell_sat,
                        aerial_red_cdf, aerial_green_cdf, aerial_blue_cdf, ncell_aerial)
    hist_match_params_mat <- rbind(hist_match_params_mat, masked_list)
  }
  
  #write cdf/icdf matrix as RDS + grid layers to geopackage
  saveRDS(hist_match_params_mat, paste(out_name, '.rds', sep=''))
  write_sf(cdf_grid, paste(out_name, '.gpkg', sep=''), layer='cdf_grid', 
           fid_column_name='cdf_id', layer_options='FID=cdf_id') 
  interp_grid_for_write <- interp_grid
  interp_grid_for_write$cdf_neighbors <- as.character(interp_grid_for_write$cdf_neighbors[1:nrow(interp_grid_for_write)])
  write_sf(interp_grid_for_write, paste(out_name, '.gpkg', sep=''), layer='interp_grid', 
           fid_column_name='interp_id', layer_options='FID=interp_id') 
  
  ##find buffered cdf_grid centers. need this because terra::distance doesn't work with points
  #cdf_grid_center_polys <- st_buffer(cdf_grid_centers, CDF_GRID_CENTER_BUFFER_SIZE)
  
  #rasterize center points so we can use terra::distance (might be bugged for vector data)
  
  
  #hist match and interpolate in each interp_grid cell
  print('Matching aerial to satellite')
  for (row in 1:nrow(interp_grid)) {
    print(paste(Sys.time(), 'Begin matching', row))
    interp_cell <- interp_grid[row,]
    interp_type <- interp_cell$interp_type
    inv_distance_rast_set <- list() #will be multi-layer raster
    matched_aerial_list <- list()
    
    #label neighbor points in aerial
    #print(paste(Sys.time(), 'Labeling neighbor points', row))
    cdf_neighbor_points <- cdf_grid_centers[cdf_grid_centers$cdf_id %in% interp_cell$cdf_neighbors[[1]],]
    
    #sort cdf neighbors by geometries: ll, lr, ul, ur
    #this seems to already be the case, but this ensures that
    #print(paste(Sys.time(), 'Getting cdf points geom info', row))
    neighbor_geoms <- geom(vect(cdf_neighbor_points))
    ordered_cdf_points <- cdf_neighbor_points[order(neighbor_geoms[,'y'], neighbor_geoms[,'x']),]
    
    #crop aerial to interp_cell
    #print(paste(Sys.time(), 'Cropping aerial to cell', row))
    aerial_interp_cropped <- crop(aerial, vect(st_geometry(interp_cell))) #terra doesn't like neighbors field
    
    #make relative distance rasters for interpolation, if this cell isn't a corner cell 
    if (interp_type == 'ordinary') {
      
      #find pixel coordinates
      #print(paste(Sys.time(), 'Making pixel coordinate matrix', row))
      aerial_coords <- crds(aerial_interp_cropped)
      #print(paste(Sys.time(), 'Extracting stats from pixel coordinate matrix', row))
      aerial_x_coords <- aerial_coords[,1]
      aerial_y_coords <- aerial_coords[,2]
      min_x_image <- min(aerial_x_coords)
      min_y_image <- min(aerial_y_coords)
      max_x_image <- max(aerial_x_coords)
      max_y_image <- max(aerial_y_coords)
      
      #make x/y coordinate rasters
      #print(paste(Sys.time(), 'Making x and y coordinate rasters', row))
      x_coords_rast <- rast(nrows=nrow(aerial_interp_cropped), ncols=ncol(aerial_interp_cropped),
                            extent=ext(aerial_interp_cropped), crs=crs(aerial_interp_cropped))
      x_coords_rast[] <- aerial_x_coords
      y_coords_rast <- rast(nrows=nrow(aerial_interp_cropped), ncols=ncol(aerial_interp_cropped),
                            extent=ext(aerial_interp_cropped), crs=crs(aerial_interp_cropped))
      y_coords_rast[] <- aerial_y_coords
      
      #make relative distance rasters
      #print(paste(Sys.time(), 'Making relative distance from minimum rasters', row))
      scaled_x_coords_rast <- (x_coords_rast - min_x_image) / (max_x_image - min_x_image)
      scaled_y_coords_rast <- (y_coords_rast - min_y_image) / (max_y_image - min_y_image)
    } else if (interp_type == 'border') {
      y_range <- max(neighbor_geoms[,'y']) - min(neighbor_geoms[,'y'])
      x_range <- max(neighbor_geoms[,'x']) - min(neighbor_geoms[,'x'])
      interp_along_y <- y_range > x_range
      
      aerial_coords <- crds(aerial_interp_cropped)
      if (interp_along_y) {
        aerial_y_coords <- aerial_coords[,2]
        min_y_image <- min(aerial_y_coords)
        max_y_image <- max(aerial_y_coords)
        y_coords_rast <- rast(nrows=nrow(aerial_interp_cropped), ncols=ncol(aerial_interp_cropped),
                              extent=ext(aerial_interp_cropped), crs=crs(aerial_interp_cropped))
        y_coords_rast[] <- aerial_y_coords
        scaled_y_coords_rast <- (y_coords_rast - min_y_image) / (max_y_image - min_y_image)
      } else {
        aerial_x_coords <- aerial_coords[,1]
        min_x_image <- min(aerial_x_coords)
        max_x_image <- max(aerial_x_coords)
        x_coords_rast <- rast(nrows=nrow(aerial_interp_cropped), ncols=ncol(aerial_interp_cropped),
                              extent=ext(aerial_interp_cropped), crs=crs(aerial_interp_cropped))
        x_coords_rast[] <- aerial_x_coords
        scaled_x_coords_rast <- (x_coords_rast - min_x_image) / (max_x_image - min_x_image)
      }
    }
    
    #find hist matched aerials for each neighbor
    #print(paste(Sys.time(), 'Begin neighbor loop', row))
    for (neighbor_id in ordered_cdf_points$cdf_id) {

      #match aerials
      #print(paste(Sys.time(), 'Getting hist match params', row, neighbor_id))
      hist_match_params <- hist_match_params_mat[hist_match_params_mat[,'cdf_id']==neighbor_id,]
      #sat_n <- hist_match_params$ncell_sat
      #aerial_n <- hist_match_params$ncell_aerial
      #cdf_factor <- sat_n/(sat_L-1) * (aerial_L-1)/aerial_n #cdf's and icdf's are normalized already, this was a mistake
      
      ##red
      #print(paste(Sys.time(), 'Matching red', row, neighbor_id))
      sat_red_icdf <- hist_match_params$red_sat_icdf
      aerial_red_cdf <- hist_match_params$red_aerial_cdf
      matched_aerial_red <- app(app(aerial_interp_cropped[[1]], fun=aerial_red_cdf),
                                fun=sat_red_icdf)
      
      ##green
      #print(paste(Sys.time(), 'Matching green', row, neighbor_id))
      sat_green_icdf <- hist_match_params$green_sat_icdf
      aerial_green_cdf <- hist_match_params$green_aerial_cdf
      matched_aerial_green <- app(app(aerial_interp_cropped[[2]], fun=aerial_green_cdf),
                                  fun=sat_green_icdf)
      
      ##blue
      #print(paste(Sys.time(), 'Matching blue', row, neighbor_id))
      sat_blue_icdf <- hist_match_params$blue_sat_icdf
      aerial_blue_cdf <- hist_match_params$blue_aerial_cdf
      matched_aerial_blue <- app(app(aerial_interp_cropped[[3]], fun=aerial_blue_cdf),
                                 fun=sat_blue_icdf)
      
      #stack each matched band and append to list
      #print(paste(Sys.time(), 'Appending  matched results to list', row, neighbor_id))
      matched_aerial <- c(matched_aerial_red, matched_aerial_green, matched_aerial_blue)
      matched_aerial_list <- append(matched_aerial_list, matched_aerial)
    }

    #interpolate images. recall matched_aerial_list is sorted according to ordered_cdf_points
    print(paste('Interp type', interp_type, ', ',nlyr(matched_aerial_list), 'layers')) #TODO: drop
    if (interp_type == 'ordinary') {
      #print(paste(Sys.time(), 'Interpolating on x axis for ymin match', row))
      min_y_interpolated_matched_image <- matched_aerial_list[[1:3]] +
        (matched_aerial_list[[4:6]] - matched_aerial_list[[1:3]]) * scaled_x_coords_rast
      #print(paste(Sys.time(), 'Interpolating on x axis for ymax match', row))
      max_y_interpolated_matched_image <- matched_aerial_list[[7:9]] +
        (matched_aerial_list[[10:12]] - matched_aerial_list[[7:9]]) * scaled_x_coords_rast
      #print(paste(Sys.time(), 'Interpolating on y axis', row))
      interpolated_matched_image <- min_y_interpolated_matched_image + 
        (max_y_interpolated_matched_image - min_y_interpolated_matched_image) * scaled_y_coords_rast
    } else if (interp_type == 'corner' && nlyr(matched_aerial_list)==3) {
      interpolated_matched_image <- matched_aerial_list[[1:3]] #condition asserts there's only 1 image in list 
    } else if (interp_type == 'border') {
      if (interp_along_y) {
        interpolated_matched_image <- matched_aerial_list[[1:3]] + 
          (matched_aerial_list[[4:6]] - matched_aerial_list[[1:3]]) * scaled_y_coords_rast
      } else {
        interpolated_matched_image <- matched_aerial_list[[1:3]] + 
          (matched_aerial_list[[4:6]] - matched_aerial_list[[1:3]]) * scaled_x_coords_rast
      }
    } else {
      next #TODO: figure out what to do for hinges
    }
    
    ##append interpolated/matched image to mosaic
    ##print(paste(Sys.time(), 'Merging results with mosaic', row))
    #if (!exists("interpolated_matched_image_mosaic")) {
    #  interpolated_matched_image_mosaic <- interpolated_matched_image
    #} else {
    #  interpolated_matched_image_mosaic <- merge(interpolated_matched_image_mosaic, 
    #                                              interpolated_matched_image)
    #}
    
    #if (row==200) {remove(interpolated_matched_image_mosaic);break} #TODO: drop this
    
    #mosaicing is making things take too long. let's build write a tif for each file
    #then make a vrt 
    out_path_for_cell <- file.path(scratch_dir_path, paste('match', interp_cell$interp_id, '.tif', sep=''))
    writeRaster(interpolated_matched_image, out_path_for_cell, overwrite=TRUE)
  }
  
  #write merged output raster
  print('Making vrt of scratch tiles')
  output_tile_list <- list.files(path = scratch_dir_path, 
                                 pattern = paste('match[0-9]+.tif$', sep=''),
                                 full.names = TRUE)
  matched_aerial_vrt <- vrt(output_tile_list)
  
  print('Writing vrt to GeoTiff')
  writeRaster(matched_aerial_vrt, paste(out_name, '.tif', sep=''), overwrite=TRUE) 
  
  #remove results from scratch dir 
  print('Deleting scratch tiles')
  unlink(scratch_dir_path, recursive = TRUE)
}

