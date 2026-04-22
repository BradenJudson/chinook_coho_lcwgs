library(shiny); library(plotly); library(readr);  library(tidyverse)
library(ggnewscale); library(shinyjs); library(heatmaply); library(sf)

coastline <- st_read("data/simplified_coastline.shp")

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  useShinyjs(),
  
  titlePanel("Salmon lcWGS population structure PCA"),
  
  tabsetPanel(
    tabPanel("PCA", fluid = TRUE,
             
             # Establish sidebar tools that specify if the data are from 
             # genotype likelihoods or imputed genotypes, and whether the
             # data from each group are connected by lines.
             sidebarLayout(
               sidebarPanel(
                 radioButtons("species", "Species", choices = c("Chinook", "Coho")),
                 radioButtons("genotypes", "Genotypes",
                              choices = c("Imputed", "Likelihoods")),
                 radioButtons("stars",
                              "Connect points by population?",
                              choices = c("No", "Yes")),
                 # Input for selecting X variable
                 selectInput("xcol", "X-axis", choices = NULL),
                 # Input for selecting Y variable
                 selectInput("ycol", "Y-axis", choices = NULL),
                 selectInput("fill", "Fill", choices = c("population", "region"))
               ),
               
               # Panel for the main PCA figure, extended to be taller than default.
               mainPanel(
                 plotlyOutput("scatter_plot", height = '750px', width = 'auto')
               )
             )
    ),
    tabPanel("Map", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 selectInput("species", "Species", choices = c("Chinook", "Coho")),
                 radioButtons("metric", "Metric", choices = c("None",
                                                              "Genomic offset",
                                                              "Heterozygosity",
                                                              "Nucleotide diversity")),
                 radioButtons("data_type",
                              "Data type",
                              choices = c("Imputed genotypes", "Genotype likelihoods")),
                 radioButtons("off_ssp",
                              "Offset scenario",
                              choices = c("SSP2.6", "SSP4.5", "SSP8.5"))),
               mainPanel(
                 plotlyOutput("map", height = "700px", width = 'auto'))
             ))
    )
  )


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # Gray out certain options depending on previous selections.
  # I.e., cannot choose genomic offset scenario when visualizing site coordinates.
  observeEvent(input$metric, {
     
    if(input$metric == "Genomic offset") {
      enable("off_ssp")
    } else {
      disable("off_ssp")
    }
    
    if(input$metric != "None") {
      enable("data_type") }
    else {
      disable("data_type")
    }
    
  })
  
  pca_scores <- reactive({do.call("rbind", list(
    read.csv("data/chinook_imputed_pca_scores_n819.csv") %>% 
      mutate(species = "Chinook", genotypes = "Imputed"),
    read.csv("data/chinook_likelihoods_pca_scores_n819.csv") %>% 
      mutate(species = "Chinook", genotypes = "Likelihoods"),
    read.csv("data/coho_imputed_pca_scores_n650.csv") %>% 
      mutate(species = "Coho", genotypes = "Imputed"),
    read.csv("data/coho_likelihoods_pca_scores_n650.csv") %>% 
      mutate(species = "Coho", genotypes = "Likelihoods"))) %>% 
      dplyr::rename("region" = "region_revised") 
    })
  
  
  pca_dat <- reactive({
    
    df <- pca_scores()
    
    if(input$stars == "Yes") {
      
      p <- df %>% 
        filter(species == input$species & genotypes == input$genotypes)  %>% 
        group_by(population) %>% 
        summarise(across(starts_with("PC"), mean, .names = "mean_{.col}")) 
      
      j <- merge(df, p, by = "population") %>% 
        filter(species == input$species & genotypes == input$genotypes) 
      
      return(j)
      
    } else {
      
      j <- df %>% 
        filter(species == input$species)  %>% 
        filter(genotypes == input$genotypes)
      
      return(j)
      
    }
    
  })
  
  
  pvarexp <- reactive({
    
    if(input$species == "Chinook") {
      
      if(input$genotypes == "Imputed") {
        
        k <- read_tsv("data/chinook_imputed_n819.eigenval", col_names = "eigenval") %>% 
          mutate(axis = as.numeric(rownames(.)),
                 PC = paste0("PC", axis)) %>%
          arrange(axis) %>%
          mutate(var_exp = format(round(100*eigenval/sum(eigenval), 2), nsmall = 2))
        
          } else
            
        k <- data.frame(eigenval = eigen(read.table("data/chinook_angsd_n819.cov"))$values) %>% 
              rownames_to_column("PC") %>% mutate(PC = paste0("PC", PC)) %>% 
              mutate(var_exp = format(round(100*eigenval/sum(eigenval), 2), nsmall = 2))
          
           } else {
      
      if(input$genotypes == "Imputed") {
        
        k <- read_tsv("data/coho_imputed_pca_n650.eigenval", col_names = "eigenval") %>% 
          mutate(axis = as.numeric(rownames(.)),
                 PC = paste0("PC", axis)) %>%
          arrange(axis) %>%
          mutate(var_exp = format(round(100*eigenval/sum(eigenval), 2), nsmall = 2)) } else
            
        k <- data.frame(eigenval = eigen(read.table("data/coho_angsd_n650.cov"))$values) %>% 
              rownames_to_column("PC") %>% mutate(PC = paste0("PC", PC)) %>% 
              mutate(var_exp = format(round(100*eigenval/sum(eigenval), 2), nsmall = 2))
       
    }
    
    return(k)
    
  })
  
  
  # Dynamic variable selection for fill and both axes.
  observe({
    updateSelectInput(session, "xcol", selected = "PC1",
                      choices = colnames(pca_dat() %>% 
                                           dplyr::select(starts_with("PC"))))
    updateSelectInput(session, "ycol", selected = "PC2",
                      choices = colnames(pca_dat() %>% 
                                           dplyr::select(starts_with("PC"))))
    updateSelectInput(session, "fill", selected = "region",
                      choices = c("population", "region"))
  })
  
  output$scatter_plot <- renderPlotly({
    
    # Ensure variables are selected before plotting
    req(input$xcol, input$ycol) 
    
    # Construct base plot.
    # Optional geom_segment layer depending on arguments above.
    base_plot <- plotly::ggplotly(
      ggplot() +
        {if(input$stars == "Yes") geom_segment(data = pca_dat(),
                                               aes(x = .data[[paste0("mean_", input$xcol)]],
                                                   xend = .data[[input$xcol]],
                                                   y = .data[[paste0("mean_", input$ycol)]],
                                                   yend = .data[[input$ycol]],
                                                   colour = .data[[input$fill]]),
                                               size = 1, show.legend = FALSE,
                                               alpha = 0.8, inherit.aes = FALSE)} +
                                    geom_point(data = pca_dat(),
                                               shape = 21, stroke = 1/10, size = 2,
                                               aes(x = .data[[input$xcol]],
                                                   y = .data[[input$ycol]],
                                                   fill = .data[[input$fill]],
                                                   text = paste0("Region: ", region,"\n",
                                                                 "Population: ", population, "\n",
                                                                 "ID: ", id))) +
                                    theme_classic() +
                                    theme(legend.title = element_blank()) +
                                    scale_color_discrete(guide = "none") +
                                    guides(fill = guide_legend(override.aes = list(size = 2.5, shape = 21), ncol = 1), 
                                           colour = guide_legend(override.aes = list(size = 2.5, shape = 21), ncol = 1)) +
                                    labs(x = paste0(input$xcol, " (", pvarexp()$var_exp[which(pvarexp()$PC == input$xcol)] , "%)"),
                                         y = paste0(input$ycol, " (", pvarexp()$var_exp[which(pvarexp()$PC == input$ycol)] , "%)"),
                                         fill = NULL),
                                  tooltip = "text") 
    
    
    # Plotly merges the fill/colour aesthetics above, resulting in brackets,1 around
    # each legend label, so this loop awkwardly fixes that labelling issue.
    for (i in 1:length(base_plot$x$data)){
      if (!is.null(base_plot$x$data[[i]]$name)){
        base_plot$x$data[[i]]$name = gsub('^\\(|,\\d+\\)$', '', base_plot$x$data[[i]]$name)
      }
    }
    
    base_plot
    
    
  })
  

# ------------------------------------------------------------------------------

sites <- readxl::read_excel("data/sample_sites.xlsx")
  sites$population <- sites$site
  sites$region <- sites$region_revised
  sites$species <- tools::toTitleCase(sites$species)

  
offs <- do.call("rbind", list(
  read.csv("data/chinook_offset_19bio3SD_afsGT_2pca_n819.csv") %>% 
    mutate(species = "Chinook", genotypes = "Imputed genotypes"),
  read.csv("data/coho_offset_19bio3SD_afsGTs_2pca_n650.csv") %>% 
    mutate(species = "Coho", genotypes = "Imputed genotypes"),
  read.csv("data/chinook_offset_19bio3SD_afsGL_2pca_n819.csv") %>% 
    mutate(species = "Chinook", genotypes = "Genotype likelihoods"),
  read.csv("data/coho_offset_19bio3SD_afsGLs_2pca_n650.csv") %>% 
    mutate(species = "Coho", genotypes = "Genotype likelihoods")
)) %>% pivot_longer(cols = starts_with("go"), names_to = "scenario")

hets <- do.call("rbind", list(
  read.csv("data/chinook_imputed_het_n106.csv") %>% 
    mutate(species = "Chinook", genotypes = "Imputed genotypes"),
  read.csv("data/coho_imputed_het_n83.csv") %>% 
    mutate(species = "Coho", genotypes = "Imputed genotypes"),
  read.csv("data/chinook_lcwgs_hets_n106.csv", row.names = 1) %>% 
    mutate(species = "Chinook", genotypes = "Genotype likelihoods"),
  read.csv("data/coho_lcwgs_hets_n83.csv", row.names = 1) %>% 
    mutate(species = "Coho", genotypes = "Genotype likelihoods")
))

nuc <- do.call("rbind", list(
  read.csv("data/chinook_imputed_pi_n106.csv") %>% 
    mutate(species = "Chinook", genotypes = "Imputed genotypes"),
  read.csv("data/coho_imputed_pi_n83.csv") %>% 
    mutate(species = "Coho", genotypes = "Imputed genotypes"),
  read.csv("data/coho_likelihoods_pi_n83.csv") %>% 
    mutate(species = "Coho", genotypes = "Genotype likelihoods"),
  read.csv("data/chin_likelihoods_pi_n106.csv") %>% 
    mutate(species = "Chinook", genotypes = "Genotype likelihoods")
)) %>% dplyr::rename("site" = "pop")


gen_div <- merge(hets, nuc, 
                 by = c("species", "genotypes", "site"), 
                 all.x = T) %>% 
  pivot_longer(cols = c("meanPi", "het")) 


fill_dat <- reactive({
  
  if(input$metric == "None") {
    
    metdf <- sites %>% 
      mutate(species = tools::toTitleCase(species)) %>% 
      filter(species == input$species) %>% 
      dplyr::select(c("population", "latitude", "longitude", "region")) %>% 
      mutate(region = fct_reorder(region, latitude, .fun = mean)) %>% 
      mutate(metric = 1) 
    return(metdf)
    
  }
  
  if(input$metric == "Genomic offset") {
    
    SSP <- if(input$off_ssp == "SSP2.6") "go26" else if (input$off_ssp == "SSP4.5") "go45" else "go85" 
    
    metdf <- offs %>% 
      filter(species == input$species) %>% 
      filter(genotypes == input$data_type) %>% 
      filter(scenario == SSP) %>% 
      merge(., sites, by = c("site")) %>% 
      rename("metric" = value) %>% 
      dplyr::select(c("site", "latitude", "longitude", "metric", "region")) %>% 
      mutate(region = fct_reorder(region, latitude, .fun = mean))
    return(metdf)
    
  }
  
  if(input$metric %in% c("Heterozygosity", "Nucleotide diversity")) {
    
    metdf <- gen_div %>% 
      filter(species == input$species) %>% 
      filter(genotypes == input$data_type) %>% 
      filter(name == if(input$metric == "Heterozygosity") "het" else "meanPi") %>% 
      merge(., sites[sites$species == input$species, ], by = c("site")) %>% 
      rename("metric" = value) %>% 
      dplyr::select(c("site", "latitude", "longitude", "metric", "region")) %>% 
      mutate(region = fct_reorder(region, latitude, .fun = mean))
    return(metdf)
    
  }
  
})

  
  output$map <- renderPlotly({
    
     mp <- mean(fill_dat()$metric)
    
    (offmap <- plotly::ggplotly(ggplot() +
                                  geom_sf(data = coastline) +
                                  theme_classic() +
                                  {if (input$metric == "None") 
                                    geom_point(data = fill_dat(),
                                               aes(x = longitude, y = latitude, fill = region,
                                                   text = paste0(population, "\n", "Latitude: ", latitude,
                                                                 "\n", "Longitude: ", longitude,
                                                                 "\n", "Region: ", region)),
                                               size = 2, colour = 'black', shape = 21) else
                                                 geom_point(data = fill_dat(),
                                                            aes(x = longitude,
                                                                y = latitude,
                                                                size = metric,
                                                                fill = metric,
                                                                text = paste0(site, "\n",
                                                                              input$metric, ": ", round(metric, 3))),
                                                            shape = 21, stroke = 1/2)
                                    } +
                                  coord_sf(xlim = c(min(sites$longitude), max(sites$longitude)),
                                           ylim = c(min(sites$latitude),  max(sites$latitude))) +
                                  labs(x = NULL, y = NULL) +
                                  {if (input$metric != "None") scale_fill_gradient2(low = if(input$metric ==  "Genomic offset") 'skyblue' else "red2",
                                                                                    high = if(input$metric ==  "Genomic offset") 'red2' else "skyblue",,
                                                                                    midpoint = mp,
                                                                                    name = input$metric,
                                                                                    limits = if(input$metric == "Genomic offset") range(offs[offs$species == input$species, "value"]),
                                                                                    breaks = scales::breaks_pretty(n = 5)) } +
                                  scale_size_continuous(range = c(1,4),
                                                        name = 'Offset',
                                                        if(input$metric == "Genomic offset") range(offs[offs$species == input$species, "value"])) +
                                  theme(panel.grid = element_blank()),
                                tooltip = "text"))
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
