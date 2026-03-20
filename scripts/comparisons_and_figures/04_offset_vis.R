library(tidyverse); library(sf); library(geodata); library(ggnewscale)


# Assemble datasets ------------------------------------------------------------

chin_sites <- read.csv("data/chin_site_info.csv")
chin_off   <- read.csv("data/GOs/chinook_offset_19bio3SD_afsGT_2pca_n819.csv")
chinook <- merge(chin_sites, chin_off)

coho_sites <- read.csv("data/coho_site_info.csv")
coho_off   <- read.csv("data/GOs/coho_offset_19bio3SD_afsGTs_2pca_n650.csv")
coho <- merge(coho_sites, coho_off)


# Spatial info -----------------------------------------------------------------

# Load in shape files. Already downloaded in the map folder.
USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "../map/"))
CAN <- sf::st_as_sf(geodata::gadm(country = "CAN", level = 0, path = "../map/"))

pcd <- data.frame(lon = c(-110, -190), lat = c(30, 80)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_bbox() %>% st_as_sfc() %>% st_cast("POLYGON")

sf_use_s2(FALSE) # Required for st_difference on class MULTIPOYLGON(s).
# Remove Canadian landmass from the above polygon.
oceanp1 <- st_difference(x = pcd, y = CAN, dimension = "polygon")
# Remove the American landmass from the above polygon.
ocean   <- st_difference(x = oceanp1, y = USA, dimension = "polygon")
# Borders become incomplete if done in a single step.
rm(oceanp1)

bound <- terra::ext(c(-180, -110, 30, 70))

pres_bio <- worldclim_global(var = "bio", res = 5,
                             path = "data/wcdata")

f26_bio <- cmip6_world(model = "UKESM1-0-LL", ssp = "126", time = "2041-2060",
                       path = "data/fclim26/", res = 5, var = "bio")

f45_bio <- cmip6_world(model = "UKESM1-0-LL", ssp = "245", time = "2041-2060",
                       path = "data/fclim45/", res = 5, var = "bio")

f85_bio <- cmip6_world(model = "UKESM1-0-LL", ssp = "585", time = "2041-2060",
                       path = "data/fclim85/", res = 5, var = "bio")

rivers <- read_sf("data/keep_rivers_2.shp")

# Overlay function -------------------------------------------------------------

offset_plot <- \(df, variable, cgr_rev, plot_title, abs_cols, brks, rast, LP, rb) {

  if(rb == TRUE) bound <- terra::ext(c(min(df$Longitude)-5, max(df$Longitude)+5,
                                       min(df$Latitude)-2,  max(df$Latitude)+2)) else
                 bound <- terra::ext(c(-180, -110, 30, 70))
                                       
  if(is.null(rast) == FALSE) range <- c(round(terra::minmax(terra::crop(f26_bio[[names(rast)]], bound))[1]), # High mitigation minimum.
                                        round(terra::minmax(terra::crop(f85_bio[[names(rast)]], bound))[2])) # Low mitigation maximum.
                                         
  # Same as above, but for a shared/independent middle value for the legend.
  if(abs_cols == TRUE) mp <- median(c(df$go26, df$go85)) else
    mp <- mean(as.numeric(df[,deparse(substitute(variable))]))

  if(is.null(rast) == TRUE) cols <- "gray70" else cols <- c("#35A052", "#FCFEBB", "red")

  # Plot land and populations on it.
  p1 <- ggplot(data = df) + theme_bw() +
    {if(is.null(rast) == FALSE)  tidyterra::geom_spatraster(data = terra::crop(rast, bound)) } +
      scale_fill_gradientn(colours = cols, 
                           limits = range,
                           oob = scales::squish, name = names(rast),
                           guide = guide_colorbar(frame.colour = 'black',
                                                  ticks.colour = 'black')) +
    new_scale_color() + new_scale_fill() +
    geom_sf(data = ocean, fill = "white", linewidth = 1/10, colour = 'black') +
    labs(x = NULL, y = NULL) +
    geom_sf(data = rivers, linewidth = 1/10, colour = "gray70", fill = "gray70") +
    geom_point(aes(x = Longitude,
                   y = Latitude,
                   fill  = {{variable}},   # Colour is metric of interest.
                   size  = {{variable}}),  # Point size is also proportional to value.
               shape = 21, stroke = 1/2) + 
    ggtitle(plot_title) +
    theme(legend.background = element_rect(colour = 'black', fill = 'white'),
          legend.position = "right",
          legend.justification = "top",
          panel.grid = element_blank(),
          legend.key = element_blank(),
          panel.background = element_rect(fill = "gray90")) +
    coord_sf(xlim = c(max(df$Longitude), min(df$Longitude)),
             ylim = c(min(df$Latitude),  max(df$Latitude))) +
    scale_fill_gradient2(low  = if(cgr_rev == FALSE) "skyblue" else "red2",
                         high = if(cgr_rev == FALSE) "red2" else "skyblue",
                         midpoint = mp, 
                         name = if(is.null(rast) == FALSE) paste0("Genomic", "\n", "offset") else "Heterozygosity",
                         breaks   = scales::breaks_pretty(n = brks),
                         limits   = if(abs_cols == TRUE) c(min(df[,"go26"]), max(df[,"go85"])) else
                           c(min(df[,deparse(substitute(variable))]),
                             plyr::round_any(max(df[,deparse(substitute(variable))]),
                                             accuracy = 0.01, f = ceiling))) +
    scale_size_continuous(breaks = scales::breaks_pretty(n = brks),
                          name = if(is.null(rast) == FALSE) paste0("Genomic", "\n", "offset") else "Heterozygosity",
                          range = c(1, 4),
                          limits   = if(abs_cols == TRUE) c(min(df[,"go26"]), max(df[,"go85"])) else
                            c(min(df[,deparse(substitute(variable))]),
                              plyr::round_any(max(df[,deparse(substitute(variable))]),
                                              accuracy = 0.01, f = ceiling))) +
    guides(fill = guide_legend(order = 1),
           size = guide_legend(order = 1)) +
    theme(legend.position = "inside",
          legend.position.inside = LP,
          legend.justification = "left",
          legend.box.just = "left",
          legend.spacing.y = unit(2, "pt"))
  # Above ensures size and fill legends are combined.
  
  p1

}


# Coho offsets -----------------------------------------------------------------

(coho85 <- offset_plot(coho, variable = go85, rast = f85_bio$bio05,
                      cgr_rev = FALSE, abs_cols = T, brks = 4, 
                      plot_title = "Genomic offset (ssp85)",
                      LP = c(0.01, 2/3), rb = FALSE))
ggsave("plots/coho_GO85_19bio3RDA3SD2PCA_afsGT_n650.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

(coho45 <- offset_plot(coho, variable = go45, rast = f45_bio$bio05,
                       cgr_rev = FALSE, abs_cols = T, brks = 4,
                       plot_title = "Genomic offset (ssp45)",
                       LP = c(0.01, 2/3), rb = FALSE))
ggsave("plots/coho_GO45_19bio3RDA3SD2PCA_afsGT_n650.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

(coho26 <- offset_plot(coho, variable = go26, rast = f26_bio$bio05,
                       cgr_rev = FALSE, abs_cols = T, brks = 4,
                       plot_title = "Genomic offset (ssp26)",
                       LP = c(0.01, 2/3), rb = FALSE))
ggsave("plots/coho_GO26_19bio3RDA3SD2PCA_afsGT_n650.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)


# Chinook offsets --------------------------------------------------------------


(chin85 <- offset_plot(chinook, variable = go85, rast = f85_bio$bio05,
                       cgr_rev = FALSE, abs_cols = T, brks = 4,
                       plot_title = "Genomic offset (ssp85)",
                       LP = c(0.01, 0.3), rb = FALSE))
ggsave("plots/chinook_GO85_19bio3RDA3SD2PCA_afsGT_n819.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

(chin45 <- offset_plot(chinook, variable = go45, rast = f45_bio$bio05,
                       cgr_rev = FALSE, abs_cols = T, brks = 4,
                       plot_title = "Genomic offset (ssp45)",
                       LP = c(0.01, 0.3), rb = FALSE))
ggsave("plots/chinook_GO45_19bio3RDA3SD2PCA_afsGT_n819.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

(chin26 <- offset_plot(chinook, variable = go26, rast = f26_bio$bio05,
                       cgr_rev = FALSE, abs_cols = T, brks = 4,
                       plot_title = "Genomic offset (ssp26)",
                       LP = c(0.01, 0.3), rb = FALSE))
ggsave("plots/chinook_GO26_19bio3RDA3SD2PCA_afsGT_n819.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)


# # Coho heterozygosity ----------------------------------------------------------



coho_hets <- read.csv("data/heterozygosity_files/coho_lcwgs_hets_n650.csv",
                      row.names = 1) %>%
  merge(., coho_sites, by.x = "pop", by.y = "site")

(cohoHet <- offset_plot(coho_hets, variable = het, rast = NULL,
                        cgr_rev = TRUE, abs_cols = F, brks = 8,
                        plot_title = NULL,
                        LP = c(0.04, 0.15),
                        rb = TRUE))
ggsave("plots/coho_lcwgs_het.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)



chin_hets <- read.csv("data/heterozygosity_files/chinook_lcwgs_hets_n819.csv",
                      row.names = 1) %>% 
  merge(., chin_sites, by.x = "pop", by.y = "Site")

(chinHet <- offset_plot(chin_hets, variable = het, rast = NULL,
                        cgr_rev = TRUE, abs_cols = F, brks = 5,
                        plot_title = NULL,
                        LP = c(0.04, 0.15),
                        rb = TRUE))
ggsave("plots/chin_lcwgs_het.tiff", bg = 'white',
       dpi = 300, width = 8, height = 6)

sites <- readxl::read_excel("data/sample_sites.xlsx")
expand_same <- \(x) x + 
  coord_sf(xlim = c(max(sites$longitude), min(sites$longitude)),
           ylim = c(min(sites$latitude),  max(sites$latitude))) +
  theme(legend.position = "top",
        legend.background = element_rect(color = NA)) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE),
         size = guide_legend(nrow = 1, byrow = TRUE))

(sp_het <- cowplot::plot_grid(plotlist = list(
  expand_same(chinHet),
  expand_same(cohoHet)),
  ncol = 2, nrow = 1
))
ggsave("plots/lcwgs_het.tiff", bg = 'white',
       dpi = 300, width = 12, height = 6)
