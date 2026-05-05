library(sp); library(terra); library(sf); library(geodata)
library(viridis); library(tidyverse); library(pheatmap)


coho <- read.csv("data/coho_site_info.csv", row.names = 1)[,c(2,5:6)] %>% 
  filter(!site %in% c("CapilanoEarly", "Seton", "Smith"))
chin <- read.csv("data/chin_site_info.csv")

getbioclim <- \(sp) {
  
  # Establish a blank list for populating.
  c.list <- list()
  
  sites <- sp
  
  # For each site, extract bioclimatic data.
  for (i in 1:nrow(sites)) {
    
    # Isolate site-specific lat and long.
    lon <- sites[i, "Longitude"]; lat <- sites[i, "Latitude"]
    
    # Download data from WorldClim at 5 arcsmin.
    # Store downloaded *.tiff files locally.
    # First, download current (1970 - 2000) conditions.
    c.clim <- geodata::worldclim_tile(var  = "bio",
                                      lon  = lon, 
                                      lat  = lat, 
                                      res  = 5,
                                      path = "./data/wcdata",
                                      download = FALSE)
    
    print(paste("Population", i, "out of", nrow(sites)))
    
    # Download projected bioclimatic variables.
    # Following Tigano et al. 2023 [Evo. App.]
    # Climate data at RCP2.6 (~best case).
    # Set download = TRUE for first time.
    f.clim26 <- cmip6_tile(ssp  = "126",
                           lon  = lon,
                           lat  = lat,
                           res  = 5,
                           path = "./data/fclim26",
                           time = "2041-2060",
                           var  = "bioc",
                           model = "UKESM1-0-LL",
                           download = FALSE)
    
    f.clim45 <- cmip6_tile(ssp  = "245",
                           lon  = lon,
                           lat  = lat,
                           res  = 5,
                           path = "./data/fclim45",
                           time = "2041-2060",
                           var  = "bioc",
                           model = "UKESM1-0-LL",
                           download = FALSE)
    
    
    # Similar to above, but RCP8.5 (~worst case).
    f.clim85 <- cmip6_tile(ssp  = "585",
                           lon  = lon,
                           lat  = lat,
                           res  = 5,
                           path = "./data/fclim85",
                           time = "2041-2060",
                           var  = "bioc",
                           model = "UKESM1-0-LL",
                           download = FALSE)
    
    # Extract site coordinates and project.
    points <- vect(sites[,c("Longitude", "Latitude")], crs = "EPSG:4326",
                   geom = c("Longitude", "Latitude"))
    
    # Set up a function to extract data from each of the 
    # three climate datasets isolated above.
    # Also renames columns consistently and adds two qualifiers.
    f <- function(clim, ssp, time) cbind(sites[i, 1],
                                         terra::extract(clim, y = points)[i, 2:20]) %>% 
      `colnames<-`(., c("Site", paste0("bio", seq(1, 19, 1)))) %>% 
      mutate(period = time, ssp = ssp) 
    
    future_data <- rbind(
      f(f.clim26, time = "2041-2060", ssp = "26"),
      f(f.clim45, time = "2041-2060", ssp = "45"),
      f(f.clim85, time = "2041-2060", ssp = "85"))
    
    # For each site, bind together and specify time period and ssp.
    # NOTE: First (1970s-2000) data set is in a different order!!!
    # bio1, bio10, bio11,...bio9. Unlike the other two. 
    # Arrange separately otherwise there is a severe mismatch.
    contemporary_data <- cbind(sites[i, 1], terra::extract(c.clim, points)[i, 2:20]) %>%
      `colnames<-`(., c("Site", paste0("bio", c(1:19)))) %>% 
      mutate(period = "1970-2000", ssp = NA)
    
    c.list[[i]] <- rbind(future_data, contemporary_data,
                         make.row.names = FALSE)
    
  }
  
  # Put the bioclim data into a single dataframe. 
  biovar <- bind_rows(c.list); head(biovar)
  return(biovar)
  
}


co_bio <- getbioclim(coho)
write.csv(co_bio, "data/coho_bioclim.csv")

ch_bio <- getbioclim(chin)
write.csv(ch_bio, "data/chinook_bioclim.csv")

j <- read.csv("data/chinook_bioclim_n106.csv") %>% 
  filter(is.na(ssp)) %>% 
  dplyr::select(-c(period, ssp)) %>% 
  column_to_rownames("site")
j2 <- prcomp(j, scale. = T)
screeplot(j2, type = "barplot", main = NULL)
var_explained <- j2$sdev^2 / sum(j2$sdev^2)

scatterplot3d(j[1, 1:3])
