library(ggplot2); library(tidyverse); library(viridis); library(sf)
library(ggspatial); library(sp); library(ggmapinset); library(ggrepel)
library(scales); library(cowplot)

# Use spherical geometries.
sf_use_s2(T)

# Site info --------------------------------------------------------------------
 
sites <- readxl::read_excel("data/sample_sites.xlsx") %>% 
  group_by(species, region_revised) %>% 
  mutate(mLat = mean(latitude)) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  # arrange(mLat) %>% 
  arrange(latitude, .by_group = TRUE) %>% 
  arrange(mLat) %>% 
  mutate(sitenum = row_number(),
         alpha = case_when( region_revised %in% unique(region_revised)[seq(from = 1, to = length(unique(region_revised)), by = 2)] ~ 0.9,
                           !region_revised %in% unique(region_revised)[seq(from = 1, to = length(unique(region_revised)), by = 2)] ~ 1.0))



#  Arrange regions by average latitude.
lin <- sites %>% group_by(region_revised) %>% 
  summarize(mLat = mean(latitude)) %>% 
  arrange(mLat)


sites$region_revised <- factor(sites$region_revised, levels = lin$region_revised)
sites$species <- tools::toTitleCase(sites$species)

sitesf <- st_as_sf(x = sites, crs = 4326,
                    coords = c("longitude", "latitude"))

bound <- c(ymin = min(sites$latitude)-5, ymax = max(sites$latitude)+5, 
           xmin = min(sites$longitude)-5, xmax = max(sites$longitude)+5)
sf_use_s2(TRUE) # Necessary


# River shapefile --------------------------------------------------------------

# From the hydroRivers database: https://www.hydrosheds.org/products/hydrorivers
# Stitched together the NA and Arctic lines, set to a fixed with multipolygon
river_poly <- st_read("../map/rivers_lat25to70_ord4.shp")
keep_poly <- sf::st_is_within_distance(x = sitesf, y = river_poly, dist = 100, sparse = T)
keep_rivers <- river_poly[c(unique(unlist(keep_poly))),]
write_sf(keep_rivers, "data/keep_rivers_2.shp") 
keep_rivers <- st_read("data/keep_rivers_2.shp")
ggplot() + geom_sf(data = keep_rivers)

# Land features ----------------------------------------------------------------

# Outlines of the USA and Canada.
USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 1, path = "../map/"))
prov <- sf::st_as_sf(geodata::gadm(country = "CAN", level = 1, path = "../map/"))

# Multipolygon of land features.
land <- st_crop(rbind(USA, prov), y = bound) # For downstream plotting purposes.
land$NAME_1 <- factor(land$NAME_1, levels = c("British Columbia", "Alaska"))
rm(USA); rm(prov); gc() # Speeds this up.

# Inset map --------------------------------------------------------------------

(scbc <- sites %>% 
   filter(latitude > 49 & latitude < 56 & longitude < -120 & longitude > -135) %>% 
   st_as_sf(., crs = 4326, coords = c("longitude", "latitude")) %>% 
   st_union() %>% st_centroid() %>% st_cast(., "POINT") %>% st_make_valid())

# The inset will 2.5x larger than its base size and will be situated
# 1800km west and 600km south of the centroid of the point being inset.
scbc_insert <- configure_inset(
  shape_rectangle(centre = scbc, 
                  hwidth = 320, hheight = 400),
  scale = 2.5, units = "km",
  translation = c(-1400, -700)
)

coldf <- data.frame(
  region = unique(as.character(sites$region_revised)),
  alpha = c(rep(c(0.3, 1), length(unique(sites$region_revised))/2), 0.3),
  colour = hue_pal()(length(unique(sites$region_revised)))
)

(main <- ggplot(data = sitesf) +
    # ggnewscale::new_scale_fill() +
    geom_sf_inset(data = st_make_valid(land),
                  linewidth = 1/20,
                  aes(fill = COUNTRY),
                  show.legend = FALSE,
                  inherit.aes = FALSE) +
    scale_fill_manual(name = "Country", guide = "none",
                      values = c("gray95", "gray80")) +
    geom_sf_inset(data = keep_rivers,
                  fill = 'skyblue',
                  colour = 'skyblue',
                  linewidth = 1/10,
                  inherit.aes = FALSE,
                  show.legend = FALSE) +
    geom_sf_inset(data = st_make_valid(land), 
                  linewidth = 1/20,
                  fill = NA,
                  show.legend = FALSE) +
    geom_sf_inset(map_base = "clip", size = NA) +
    ggnewscale::new_scale_fill() +
    geom_inset_frame(lines.aes = list(linetype = 2)) +
    coord_sf_inset(xlim = c(min(sites$longitude), max(sites$longitude)), 
                   ylim = c(min(sites$latitude),  max(sites$latitude)),
                   inset = scbc_insert) +
    geom_label_repel(
      aes(
        x = after_stat(x_inset),
        y = after_stat(y_inset),
        fill = region_revised,
        label = sitenum,
        geometry = geometry
      ),
      colour = 'black',
      max.overlaps = Inf,
      label.size = 1/10,
      label.r = unit(2, "pt"),
      box.padding = unit(0.1, "pt"),
      label.padding = unit(1, "pt"),
      segment.size = 1/5,
      point.size = NA,
      min.segment.length = 0,
      seed = 240,
      show.legend = F,
      size = 2,
      stat = "sf_coordinates_inset"
    ) +
    scale_fill_manual(values = alpha(c(coldf$colour), coldf$alpha)) +
    # scale_fill_manual(values = alpha(c(sites$plot_col), sites$alpha)) +
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(fill = alpha("skyblue", 1/10)),
          panel.border = element_rect(color = "black", fill = NA),
          axis.ticks = element_line(linewidth = 1/5),
          legend.key = element_rect(fill = NA, color = NA),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text  = element_text(hjust = 0, size = 12),
          legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(override.aes = list(size = 2.5),
                               ncol = 1, reverse = TRUE,
                               keyheight = 1),
           fill_group2 = "none") +
    geom_point(data = sites %>% # Exclude points not in the inset for clarity.
                 filter(latitude > 48 & latitude < 54.6 & longitude < -121.5 & longitude > -132),
               aes(x = longitude, y = latitude,
                 fill = region_revised
               ), size = 1, 
               shape = 21,
               stroke = 1/10,
               colour = 'black',
               inherit.aes = FALSE,
               show.legend = TRUE) +
    ggnewscale::new_scale_fill() +
    scale_fill_manual(values = c(sites$plot_col)) +
    labs(x = NULL, y = NULL) +
    ggspatial::annotation_scale(location = "tl", 
                                width_hint = 1/10))

(fw <- main + facet_wrap(~ species, ncol = 1))

ggsave("plots/coho_chinook_lcwgs_map_hz.tiff", dpi = 300, width = 10, height = 8, bg = "white")

# Add another inset, but this one is of North America.
# Download low-res country outlines.
ca <- map_data("world", "Canada")
us <- map_data("world", "USA") 
me <- map_data("world", "Mexico")

# Make the inset plot by itself. 
(ins <- ggplot() +
    geom_polygon(data = us, aes(x = long, y = lat, group = group),
                 fill = "grey80", colour = "black", linewidth = 1/8) +
    geom_polygon(data = ca, aes(x = long, y = lat, group = group),
                 fill = "grey90", colour = "black", linewidth = 1/8) +
    geom_polygon(data = me, aes(x = long, y = lat, group = group),
                 fill = "grey70", colour = "black", linewidth = 1/8) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", 
                                      fill = NA, linewidth = 1/4),
          panel.background = element_rect(fill = "white"))  +
    annotate("rect", fill = NA, colour = "black",
             linewidth = 1/2,
             xmin = min(sites$longitude), xmax = max(sites$longitude),
             ymin = min(sites$latitude), ymax = max(sites$latitude)) +
    coord_map(ylim = c(72, 20),
              xlim = c(-57, -185)))


# -------------------------------------------------------------------------

(map_ins <- ggdraw(fw) +
   draw_plot(
     ins,
     x = 0.66,
     y = 0.65,
     width = 0.2,
     height = 0.15
   ))

ggsave("plots/coho_chinook_lcwgs_map_inset.tiff", dpi = 300, width = 14, height = 7, bg = "white")

