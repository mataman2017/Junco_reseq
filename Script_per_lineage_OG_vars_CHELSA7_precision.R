
setwd("/media/milab/Cronus_SSD/Q_niche_modeling")

# Save the entire workspace
#save.image(file ="JAVI_data/SDM_data_LIN_OGvars_CHELSA7_senglar.RData")
## Load the saved workspace
#load("JAVI_data/SDM_data_LIN_OGvars_CHELSA2.RData")
#load("JAVI_data/SDM_data_LIN_OGvars_CHELSA2_buffer250.RData")
#load("JAVI_data/SDM_data_LIN_OGvars_CHELSA4_buffer250.RData")
#load("JAVI_data/SDM_data_LIN_OGvars_CHELSA4_buffer150_senglar.RData")
#load("JAVI_data/SDM_data_LIN_OGvars_CHELSA6_senglar.RData")


## Directories ------------------------------------------------------------
### Define dicrectories in relation to project directory
Dir.Base <-  "./"
setwd(Dir.Base)
Dir.Data <- file.path(Dir.Base, "data")
Dir.Exports <- file.path(Dir.Base, "Exports")
### Create directories which aren't present yet
Dirs <- c(Dir.Data, Dir.Exports)
CreateDir <- sapply(Dirs, function(x) if(!dir.exists(x)) dir.create(x))

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "readr","dplyr","tidyr", "taxize", "taxadb", "magrittr", "readxl", "flexsdm", "terra", "tidysdm",
  "raster", "dismo", "rnaturalearthdata", "rnaturalearth", "ggplot2", "sf", "SDMtune", "dismo", "ecospat",
  "geobuffer", "rasterVis", "car", "beepr", "usdm", "ggdendro", "readr", "gridExtra", "grid",  "PresenceAbsence"
  # names of the packages required placed here as character objects
)

sapply(package_vec, install.load.package)

library(rJava)

## Functionality ---------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function
theme_set(theme_bw())
############### FUNCTIONS POLYGON ----

create_elevation_polygon <- function(elev_raster) {
  # Ensure the input is a SpatRaster
  if (!inherits(elev_raster, "SpatRaster")) {
    stop("Input must be a SpatRaster object.")
  }
  
  # Remove values lower than 0
  elev_raster[elev_raster < 0] <- NA
  
  # Create a binary mask
  elev_mask <- elev_raster
  elev_mask[!is.na(elev_mask)] <- 1
  # NA values remain as NA
  
  # Convert the mask to polygons
  poly <- as.polygons(elev_mask, dissolve = TRUE)
  
  # Aggregate polygons to dissolve them into a single polygon
  poly_dissolved <- terra::aggregate(poly)
  
  # Convert to sf object
  poly_sf <- st_as_sf(poly_dissolved)
}


############### FILTERING DEJU LINEAGES ########

# Leer los datos (ajusta el path a tu archivo)
rawdata <- read_tsv("data/junco_hyemalis_june.csv")

# Definir los `stateProvince` para cada grupo
hyemalis_states <- c("New Hampshire", "Massachusetts", "Yukon Territory", "New York", "Michigan", "Wisconsin", "Minnesota", "Quebec", "Alberta", 
                     "Ontario","Northwest Territories" , "Manitoba" , "Nunavut",  "Newfoundland and Labrador", "New Brunswick", "Prince Edward Island" ,
                     "Saint Pierre and Miquelon","Nova Scotia", "Vermont", "Saskatchewan", "Alaska", "Connecticut") 
oreganus_states <- c("California", "Oregon", "Washington", "Idaho", "British Columbia", "Montana")
caniceps_states <- c("Utah", "Colorado")
dorsalis_states <- c("Arizona", "New Mexico")

state_filters <- list(
  HYEMALIS = hyemalis_states,
  DORSALIS = dorsalis_states,
  CANICEPS = caniceps_states,
  OREGANUS = oreganus_states
)

# Iterar sobre cada elemento de `state_filters` para crear los dataframes
for (group_name in names(state_filters)) {
  # Filtrar el dataframe `rawdata` según los `stateProvince` de cada grupo
  group_data <- rawdata %>%
    filter(stateProvince %in% state_filters[[group_name]])
  
  # Comprobar si el filtro devolvió resultados
  if (nrow(group_data) > 0) {
    # Asignar el dataframe a una variable dinámica con el nombre del grupo
    assign(group_name, group_data)
  } else {
    warning(paste("No data found for group:", group_name))
  }
}

# Lista con los nombres de los dataframes ya creados: `DORSALIS`, `CANICEPS`, `HYEMALIS`, `OREGANUS`
group_names <- c("HYEMALIS", "CANICEPS", "DORSALIS", "OREGANUS")

# Iterate over each group
for (group in group_names) {
  if (exists(group)) {  # Check if the object exists
    # Convert to sf and process
    sf_object <- st_as_sf(get(group), coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant")
    
    sf_object <- st_transform(sf_object, crs = 3857)
    
    filtered_data <- sf_object %>%
      group_by(group = st_is_within_distance(geometry, dist = 1000)) %>%
      slice(1) %>%
      ungroup() %>%
      dplyr::select(-group)
    
    coords <- st_coordinates(filtered_data)
    filtered_data$decimalLongitude <- coords[, "X"]
    filtered_data$decimalLatitude <- coords[, "Y"]
    cleaned_data <- as.data.frame(filtered_data)
    
    assign(paste0("presencia_limpia_", group), cleaned_data)
  } else {
    warning(paste("Object", group, "does not exist. Skipping."))
  }
}

presencia_limpia_CANICEPS

# Mensaje de confirmación
cat("Se han creado los dataframes con nombres presencia_limpia_HYEMALIS, presencia_limpia_CANICEPS, presencia_limpia_DORSALIS, presencia_limpia_OREGANUS con las coordenadas preservadas.\n")

#write.csv(presencia_limpia_CANICEPS, file = 'presencia_limpia_CANICEPS.csv', row.names = FALSE)
#write.csv(presencia_limpia_DORSALIS, file = 'presencia_limpia_DORSALIS.csv', row.names = FALSE)
#write.csv(presencia_limpia_OREGANUS, file = 'presencia_limpia_OREGANUS.csv', row.names = FALSE)
#write.csv(presencia_limpia_HYEMALIS, file = 'presencia_limpia_HYEMALIS.csv', row.names = FALSE)
#
presencia_limpia_CANICEPS <- read.csv("data/CSV/presencia_limpia_CANICEPS.csv")
presencia_limpia_DORSALIS <- read.csv("data/CSV/presencia_limpia_DORSALIS.csv")
presencia_limpia_OREGANUS <- read.csv("data/CSV/presencia_limpia_OREGANUS.csv")
presencia_limpia_HYEMALIS <- read.csv("data/CSV/presencia_limpia_HYEMALIS.csv")

beepr::beep()

# Ha habido un problema y se me ha dejado las coords por el camino,  las he ido a buscar después de filtrar,
# en el script lo he arreglado, pero las guardadas no tienen cordenadas. 
# correr este script si trabajo con las guardadas:

# Iterar sobre cada grupo para hacer un `join` con `rawdata` usando `gbifID`
for (group in group_names) {
  # 1. Obtener el dataframe `presencia_limpia_XXX`
  cleaned_data <- get(paste0("presencia_limpia_", group))
  
  # 2. Seleccionar la columna `gbifID` de `cleaned_data`
  gbif_ids <- cleaned_data$gbifID
  
  # 3. Hacer un join con `rawdata` para obtener las entradas completas de `rawdata` correspondientes a los `gbifID` de `cleaned_data`
  full_data <- rawdata %>%
    filter(gbifID %in% gbif_ids)  # Filtrar `rawdata` con los `gbifID` de `cleaned_data`
  
  # 4. Asignar el nuevo dataframe con la información completa a una variable con el nombre `presencia_limpia_XXX_2`
  assign(paste0("presencia_limpia_", group, "_2"), full_data)
}


presencia_limpia_HYEMALIS_2 <- presencia_limpia_HYEMALIS_2 %>%  # I think it is to remove coastal oreganus
  filter(!(stateProvince == "Alaska" & decimalLongitude > -140))

presencia_limpia_OREGANUS_2 <- presencia_limpia_OREGANUS_2 %>%  # Remove occurrences in north to avoid hyemalis
  filter(!(decimalLatitude > 58.2)) 

  presencia_limpia_OREGANUS_2 <- presencia_limpia_OREGANUS_2 %>%   # Remove occurrences in Montana to avoid mearnsi
  filter(!(stateProvince == "Montana" & decimalLatitude < 45.807616))

presencia_limpia_OREGANUS_2 <- presencia_limpia_OREGANUS_2 %>%   # Remove occurrences in Montana to avoid mearnsi
  filter(!(stateProvince == "Montana" & decimalLongitude > -112.138139))

presencia_limpia_OREGANUS_2 <- presencia_limpia_OREGANUS_2 %>%   # Remove occurrences in Montana to avoid mearnsi
  filter(!(stateProvince == "Idaho" & decimalLatitude < 43.426945))

presencia_limpia_OREGANUS_2 <- presencia_limpia_OREGANUS_2 %>%   # Remove occurrences in Montana to avoid mearnsi
  filter(!(stateProvince == "Idaho" & decimalLongitude > -112.750940))

presencia_limpia_DORSALIS_2 <- presencia_limpia_DORSALIS_2 %>%   # Remove occurrences in Montana to avoid mearnsi
  filter(!( decimalLatitude < 32.568755))


#write.csv(presencia_limpia_CANICEPS_2, file = 'data/CSV/presencia_limpia_CANICEPS_2.csv', row.names = FALSE)
#write.csv(presencia_limpia_DORSALIS_2, file = 'data/CSV/presencia_limpia_DORSALIS_2.csv', row.names = FALSE)
#write.csv(presencia_limpia_OREGANUS_2, file = 'data/CSV/presencia_limpia_OREGANUS_2.csv', row.names = FALSE)
#write.csv(presencia_limpia_HYEMALIS_2, file = 'data/CSV/presencia_limpia_HYEMALIS_2.csv', row.names = FALSE)

presencia_limpia_CANICEPS_2 <- read.csv("data/CSV/presencia_limpia_CANICEPS_2.csv")
presencia_limpia_DORSALIS_2 <- read.csv("data/CSV/presencia_limpia_DORSALIS_2.csv")
presencia_limpia_OREGANUS_2 <- read.csv("data/CSV/presencia_limpia_OREGANUS_2.csv")
presencia_limpia_HYEMALIS_2 <- read.csv("data/CSV/presencia_limpia_HYEMALIS_2.csv")

# Get map data for the US
us_map <- map_data("state")

# Plot the map with occurrences
ggplot() +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "black") +
  geom_point(data = presencia_limpia_HYEMALIS_2, 
             aes(x = decimalLongitude, y = decimalLatitude),
             color = "red", size = 2) +
  coord_fixed(1.3) +  # Preserve aspect ratio
  labs(title = "Montana Filtered Occurrences",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

# Merge presencia_sf_DOR into presencia_sf_CAN
presencia_limpia_CANICEPS_2 <-  rbind(presencia_limpia_CANICEPS_2, presencia_limpia_DORSALIS_2)
nrow(presencia_limpia_CANICEPS_2)
# Confirm the merge
cat("Merged presencia_sf_DOR into presencia_sf_CAN.\n")

group_names <- c("HYEMALIS", "OREGANUS", "CANICEPS")

# Iterar sobre cada grupo para seleccionar las columnas y crear objetos `sf`
for (group in group_names) {
  # 1. Seleccionar solo las columnas necesarias del dataframe `presencia_limpia_XXX`
  # Obtener el dataframe correspondiente
  cleaned_data <- get(paste0("presencia_limpia_", group, "_2"))
  
  # Seleccionar solo las columnas `species`, `decimalLatitude`, `decimalLongitude`
  cleaned_data <- cleaned_data %>%
    dplyr::select(species, decimalLatitude, decimalLongitude)
  
  # 2. Asignar el dataframe limpio a la variable `presencia_limpia_XXX`
  assign(paste0("presencia_limpia_", group), cleaned_data)
  
  # 3. Crear el objeto `sf` a partir del dataframe limpio
  sf_object <- st_as_sf(cleaned_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant")
  
  # 4. Asignar el objeto `sf` a una variable dinámica con el nombre `presencia_sf_XXX`
  assign(paste0("presencia_sf_", group), sf_object)
}

presencia_sf_CANICEPS
#################### LOADING YEJU DATA -----

# Load presences from the tutorial "pr_junco" for the Junco phaeonotus
presencia_limpia_YEJU_2 <- read_delim("data/yellow_eye/Junco_phaeonotus_may_june_july_august.csv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

presencia_limpia_YEJU_2$decimalLongitude
presencia_limpia_YEJU_2 <- presencia_limpia_YEJU_2 %>%
  filter(!(decimalLongitude > -94.647756)) 

presencia_limpia_YEJU <- presencia_limpia_YEJU_2[, c("species", "decimalLatitude", "decimalLongitude")]

presencia_sf_YEJU = st_as_sf(presencia_limpia_YEJU, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant")
#################### LOADING YEJU DATA GOOD

define_and_process_yeju <- function(input_file, output_name) {
  # Load the presence data
  presence_data <- read_delim(input_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # Filter and retain necessary columns
  presence_data <- presence_data %>%
    filter(!(decimalLongitude > -94.647756)) %>%
    dplyr::select(species, decimalLatitude, decimalLongitude)
  
  # Convert to an sf object
  sf_object <- st_as_sf(presence_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant")
  
  # Transform to a projected coordinate system
  sf_object <- st_transform(sf_object, crs = 3857)
  
  # Filter using spatial distance (e.g., 1000 meters)
  filtered_data <- sf_object %>%
    group_by(group = st_is_within_distance(geometry, dist = 1000)) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(-group)
  
  # Reproject back to geographic CRS (WGS84)
  filtered_data <- st_transform(filtered_data, crs = 4326)
  
  # Extract correct longitude/latitude
  coords <- st_coordinates(filtered_data)
  filtered_data$decimalLongitude <- coords[, 1]  # X-coordinates (longitude)
  filtered_data$decimalLatitude <- coords[, 2]   # Y-coordinates (latitude)
  
  # Convert back to a data frame
  cleaned_data <- as.data.frame(filtered_data)
  
  # Save the cleaned data
  assign(output_name, cleaned_data, envir = .GlobalEnv)
  return(cleaned_data)
}

# Process YEJU dataset
yeju_cleaned <- define_and_process_yeju(
  input_file = "data/yellow_eye/Junco_phaeonotus_may_june_july_august.csv",
  output_name = "presencia_limpia_YEJU"
)

presencia_sf_YEJU = st_as_sf(yeju_cleaned, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant")


#################### Visualize the taxa occurrences ----

# visualize the species occurrences
world <- ne_countries(scale = "medium", returnclass = "sf")
unique(world$region_wb)
north_am <- world %>%
  filter(region_wb== "North America" |  region_wb== "Latin America & Caribbean")

ggplot(data = north_am) +
  geom_sf() +
  geom_sf(data = presencia_sf_YEJU, color = "gray5", alpha= 0.6)+
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-179.15, -50.12), ylim = c(10,80), expand = FALSE)


#ggplot(data = north_am) +
#  geom_sf() +
#  geom_sf(data = presencia_sf_DORSALIS, color = "gray5", alpha= 0.6)+
#  xlab("Longitude") + ylab("Latitude") +
#  coord_sf(xlim = c(-179.15, -50.12), ylim = c(10,80), expand = FALSE)

ggplot(data = north_am) +
  geom_sf() +
  geom_sf(data = presencia_sf_CANICEPS, color = "gray5", alpha= 0.6)+
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-179.15, -50.12), ylim = c(10,80), expand = FALSE)

ggplot(data = north_am) +
  geom_sf() +
  geom_sf(data = presencia_sf_OREGANUS, color = "gray5", alpha= 0.6)+
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-179.15, -50.12), ylim = c(10,80), expand = FALSE)

ggplot(data = north_am) +
  geom_sf() +
  geom_sf(data = presencia_sf_HYEMALIS, color = "gray5", alpha= 0.6)+
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-179.15, -50.12), ylim = c(10,80), expand = FALSE)

nrow(presencia_sf_CANICEPS)
nrow(presencia_sf_HYEMALIS)
nrow(presencia_sf_OREGANUS)
nrow(presencia_sf_YEJU)


##################### LOADING TODAY'S ENV DATA #####

# Past data from http://www.paleoclim.org/

lista.variables <- list.files(path="JAVI_data/X_INPUTS/A_today", pattern='*tif', full.names=TRUE)
variables <- rast(lista.variables)
variables <- variables[[c("bio_10" , "bio_15",  "bio_18",  "bio_2", "bio_4", "CHELSA_TraCE21k_glz_20_V1.edit","CHELSA_TraCE21k_scd_20_V1.edit", "ndvi_CRS", "srtm_CRS",   "tree_CRS")]]
names(variables) <- c("bio_10" , "bio_15",  "bio_18",  "bio_2", "bio_4", "glz","scd", "ndvi", "elev",   "tree")
num_variables <- nlyr(variables)
names(variables)

plot(variables)

# Only for running example purposes we are going to agreggate environmental data
variables_low <- aggregate(variables, fact= 3)
beepr::beep()
names(variables_low)

# Variables bio
names(variables)
variables_noGLZ <- variables_low[[c("bio_10", "bio_15", "bio_18", "bio_2", "bio_4","scd", "elev", "ndvi", "tree")]]
variables_bio <- variables[[c("bio_10", "bio_15", "bio_18", "bio_2", "bio_4","scd", "elev")]]
variables_bio_low <- variables_low[[c("bio_10", "bio_15", "bio_18", "bio_2", "bio_4","scd", "elev" )]] 
names(variables_noGLZ)
perimeter_present_simple <- create_elevation_polygon(variables_low$elev)
perimeter_present_simple_glz <- create_elevation_polygon(variables_low$glz)
elev_df_present_simple <- fortify(perimeter_present_simple)
glz_df_present_simple  <- fortify(perimeter_present_simple_glz)

# Extract values from the raster stack
# Extract values from all cells
values <- as.data.frame(values(variables_noGLZ, na.rm = TRUE))

# Calculate correlation matrix
cor_matrix <- cor(values, use = "complete.obs")
print(cor_matrix)

# Convert correlation matrix to a distance matrix
dist_matrix <- as.dist(1 - cor_matrix)

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "average")

# Plot the dendrogram
plot(hc, main = "Dendrogram of Variables Correlation", xlab = "", sub = "", cex = 0.9)

# Calculate VIF
vif_result <- vif(variables_bio_low)
print(variables_bio_low)

##################### LOADING PAST ENV DATA #####

#21kBP 
lista.variables_21kBP <- list.files(path="JAVI_data/Y_CHELSA2/F_21kBP",pattern='*tif', full.names=TRUE)
variables_21kBP <- rast(lista.variables_21kBP)
variables_21kBP <- resample(variables_21kBP, variables_low)
names(variables_21kBP)
names(variables_21kBP) <- c("bio_2","bio_4","bio_10", "bio_15", "bio_18", "elev","glz", "scd" )
GLZ_21kBP  <- variables_21kBP[[c("glz")]] 
variables_21kBP<- variables_21kBP[[c("bio_10", "bio_15", "bio_18", "bio_2", "bio_4", "elev", "scd")]]
perimeter_21kBP <- create_elevation_polygon(variables_21kBP$elev)
perimeter_21kBP_glz <- create_elevation_polygon(GLZ_21kBP)
elev_df_21kBP <- fortify(perimeter_21kBP)
glz_df_21kBP  <- fortify(perimeter_21kBP_glz)

#4kBP 
lista.variables_4kBP <- list.files(path="JAVI_data/Y_CHELSA2/B_4kBP",pattern='*tif', full.names=TRUE)
variables_4kBP <- rast(lista.variables_4kBP)
names(variables_4kBP) <- c("bio_2","bio_4","bio_10", "bio_15", "bio_18", "elev","glz", "scd" )
variables_4kBP <- resample(variables_4kBP, variables_bio_low )
GLZ_4kBP  <- variables_4kBP[[c("glz")]] 
variables_4kBP <- variables_4kBP[[c("bio_10", "bio_15", "bio_18", "bio_2", "bio_4", "elev", "scd")]]
perimeter_4kBP <- create_elevation_polygon(variables_4kBP$elev)
perimeter_4kBP_glz <- create_elevation_polygon(GLZ_4kBP)
elev_df_4kBP <- fortify(perimeter_4kBP)
glz_df_4kBP  <- fortify(perimeter_4kBP_glz)

#8kBP 
lista.variables_8kBP <- list.files(path="JAVI_data/Y_CHELSA2/C_8kBP",pattern='*tif', full.names=TRUE)
variables_8kBP <- rast(lista.variables_8kBP)
names(variables_8kBP) <- c("bio_2","bio_4","bio_10", "bio_15", "bio_18", "elev", "glz","scd" )
variables_8kBP <- resample(variables_8kBP, variables_bio_low )
GLZ_8kBP  <- variables_8kBP[[c("glz")]] 
variables_8kBP<- variables_8kBP[[c("bio_10", "bio_15", "bio_18", "bio_2", "bio_4", "elev","scd")]] 
perimeter_8kBP <- create_elevation_polygon(variables_8kBP$elev)
perimeter_8kBP_glz <- create_elevation_polygon(GLZ_8kBP)
elev_df_8kBP <- fortify(perimeter_8kBP)
glz_df_8kBP  <- fortify(perimeter_8kBP_glz)

#12kBP 
lista.variables_12kBP <- list.files(path="JAVI_data/Y_CHELSA2/D_12kBP",pattern='*tif', full.names=TRUE)
variables_12kBP <- rast(lista.variables_12kBP)
names(variables_12kBP) <- c("bio_2","bio_4","bio_10", "bio_15", "bio_18", "elev", "glz","scd" )
variables_12kBP <- resample(variables_12kBP, variables_bio_low )
GLZ_12kBP  <- variables_12kBP[[c("glz")]] 
variables_12kBP<- variables_12kBP[[c("bio_10", "bio_15", "bio_18", "bio_2", "bio_4", "elev","scd")]]
perimeter_12kBP <- create_elevation_polygon(variables_12kBP$elev)
perimeter_12kBP_glz <- create_elevation_polygon(GLZ_12kBP)
elev_df_12kBP <- fortify(perimeter_12kBP)
glz_df_12kBP  <- fortify(perimeter_12kBP_glz)

#16kBP 
lista.variables_16kBP <- list.files(path="JAVI_data/Y_CHELSA2/E_16kBP",pattern='*tif', full.names=TRUE)
variables_16kBP <- rast(lista.variables_16kBP)
names(variables_16kBP) <- c("bio_2","bio_4","bio_10", "bio_15", "bio_18", "elev", "glz", "scd" )
variables_16kBP <- resample(variables_16kBP, variables_bio_low )
GLZ_16kBP  <- variables_16kBP[[c("glz")]] 
variables_16kBP<- variables_16kBP[[c("bio_10", "bio_15", "bio_18", "bio_2", "bio_4", "elev", "scd")]]
perimeter_16kBP <- create_elevation_polygon(variables_16kBP$elev)
perimeter_16kBP_glz <- create_elevation_polygon(GLZ_16kBP)
elev_df_16kBP <- fortify(perimeter_16kBP)
glz_df_16kBP  <- fortify(perimeter_16kBP_glz)


##################### Occurrence filtering  ----

# Lista con los nombres de los dataframes `presencia_limpia_`
group_names_2 <- c("HYEMALIS", "CANICEPS", "OREGANUS", "YEJU")
group_names_2 <- c("HYEMALIS")
group_names_2 <- c("CANICEPS")
group_names_2 <- c("OREGANUS")
group_names_2 <- c("YEJU")

# Iterar sobre cada dataframe para aplicar el script
for (group in group_names_2) {
  # 1. Obtener el dataframe correspondiente `presencia_limpia_XXX`
  presence_data <- get(paste0("presencia_limpia_", group))
  
  # 2. Añadir una columna `id` con un identificador único a cada fila
  presence_data$id <- 1:nrow(presence_data)
  
  # 3. Aplicar `occfilt_env` para filtrar los datos
  filtered_data <- presence_data %>%
    occfilt_env(
      data = .,
      x = "decimalLongitude",
      y = "decimalLatitude",
      id = "id",
      nbins = 9,  # Asegúrate de tener `num_variables` definido en tu entorno
      env_layer  = variables_noGLZ # Asegúrate de tener `variables` definido en tu entorno
    ) %>%
    # 4. Hacer un `left_join` con el dataframe original para recuperar las columnas adicionales
    left_join(presence_data, by = c("id", "decimalLongitude", "decimalLatitude"))
  
  # 5. Añadir la columna `pr_ab` con valor 1
  filtered_data$pr_ab <- 1
  
  # 6. Asignar el nuevo dataframe a una variable dinámica con el nombre `junco_pf_[NOMBRE]`
  assign(paste0("junco_pf_", group), filtered_data)
  beepr::beep()
}

junco_pf_CANICEPS
junco_pf_OREGANUS
junco_pf_HYEMALIS
junco_pf_YEJU

# Mensaje de confirmación
cat("Se han creado los dataframes `junco_pf_HYEMALIS`, `junco_pf_CANICEPS`, `junco_pf_DORSALIS`, `junco_pf_OREGANUS` y `junco_pf_YEJU` con el filtrado ambiental aplicado.\n")

beepr::beep()

##################### Delimit calibration area ----
# Asumiendo que `variables` y `north_am` ya están definidos en tu entorno

group_names_300km <- c("HYEMALIS")
group_names_250km <- c("OREGANUS", "CANICEPS")
group_names_150km <- c("YEJU")

# Iterar sobre cada grupo para aplicar el área de calibración usando `presencia_limpia_`
for (group in group_names_300km) {
  # 1. Obtener el dataframe `presencia_limpia_XXX` correspondiente
  presence_data <- get(paste0("presencia_limpia_", group))
  
  # 2. Calcular el área de calibración usando `calib_area`
  ca <- calib_area(
    data = presence_data,
    x = "decimalLongitude",
    y = "decimalLatitude",
    method = c('buffer', width = 300000),  # Buffer de 250 km
    crs = crs(variables)  # Definir el sistema de coordenadas según `variables`
  )
  
  # 3. Convertir el área de calibración a `sf`
  ca_sf <- sf::st_as_sf(ca)
  
  # 4. Convertir `north_am` a `SpatVector` para recortar el área
  north_am_corr <- terra::vect(north_am)
  
  # 5. Recortar el área de calibración usando `crop`
  ca_crop <- crop(north_am_corr, ca)
  
  # 6. Convertir el área recortada de nuevo a `sf`
  ca_sf_crop <- sf::st_as_sf(ca_crop)
  
  # 7. Asignar los objetos resultantes a nombres dinámicos
  assign(paste0("ca_sf_", group), ca_sf)
  assign(paste0("ca_sf_crop_", group), ca_sf_crop)
}

# Iterar sobre cada grupo para aplicar el área de calibración usando `presencia_limpia_`
for (group in group_names_250km) {
  # 1. Obtener el dataframe `presencia_limpia_XXX` correspondiente
  presence_data <- get(paste0("presencia_limpia_", group))
  
  # 2. Calcular el área de calibración usando `calib_area`
  ca <- calib_area(
    data = presence_data,
    x = "decimalLongitude",
    y = "decimalLatitude",
    method = c('buffer', width = 250000),  # Buffer de 250 km
    crs = crs(variables)  # Definir el sistema de coordenadas según `variables`
  )
  
  # 3. Convertir el área de calibración a `sf`
  ca_sf <- sf::st_as_sf(ca)
  
  # 4. Convertir `north_am` a `SpatVector` para recortar el área
  north_am_corr <- terra::vect(north_am)
  
  # 5. Recortar el área de calibración usando `crop`
  ca_crop <- crop(north_am_corr, ca)
  
  # 6. Convertir el área recortada de nuevo a `sf`
  ca_sf_crop <- sf::st_as_sf(ca_crop)
  
  # 7. Asignar los objetos resultantes a nombres dinámicos
  assign(paste0("ca_sf_", group), ca_sf)
  assign(paste0("ca_sf_crop_", group), ca_sf_crop)
}

# Iterar sobre cada grupo para aplicar el área de calibración usando `presencia_limpia_`
for (group in group_names_150km) {
  # 1. Obtener el dataframe `presencia_limpia_XXX` correspondiente
  presence_data <- get(paste0("presencia_limpia_", group))
  
  # 2. Calcular el área de calibración usando `calib_area`
  ca <- calib_area(
    data = presence_data,
    x = "decimalLongitude",
    y = "decimalLatitude",
    method = c('buffer', width = 150000),  # Buffer de 250 km
    crs = crs(variables)  # Definir el sistema de coordenadas según `variables`
  )
  
  # 3. Convertir el área de calibración a `sf`
  ca_sf <- sf::st_as_sf(ca)
  
  # 4. Convertir `north_am` a `SpatVector` para recortar el área
  north_am_corr <- terra::vect(north_am)
  
  # 5. Recortar el área de calibración usando `crop`
  ca_crop <- crop(north_am_corr, ca)
  
  # 6. Convertir el área recortada de nuevo a `sf`
  ca_sf_crop <- sf::st_as_sf(ca_crop)
  
  # 7. Asignar los objetos resultantes a nombres dinámicos
  assign(paste0("ca_sf_", group), ca_sf)
  assign(paste0("ca_sf_crop_", group), ca_sf_crop)
}

group_names_2 <- c("CANICEPS", "OREGANUS", "HYEMALIS", "YEJU")

# Crear un nuevo bucle para convertir los objetos `ca_sf_crop_XXX` en `SpatVector` y aplicar `crop`
for (group in group_names_2) {
  # 1. Obtener el objeto `sf` de área de calibración recortada (`ca_sf_crop_XXX`)
  ca_sf_crop <- get(paste0("ca_sf_crop_", group))
  
  # 2. Convertir el objeto `sf` a `SpatVector` para usar `crop`
  ca_crop_spat <- vect(ca_sf_crop)
  
  # 3. Usar `crop` con el `SpatVector` y la capa `variables$tree_CRS`
  ca_crop <- crop(ca_crop_spat, variables$tree)
  
  # 4. Asignar el objeto recortado a una variable dinámica `ca_crop_[NOMBRE]`
  assign(paste0("ca_crop_", group), ca_crop)
}

beepr::beep()

ggplot(data = north_am) +
  geom_sf() +
  geom_sf(data = presencia_sf_OREGANUS, color = "gray5", alpha= 0.6)+
  geom_sf(data= st_geometry(ca_sf_crop_OREGANUS), color= "red4") +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-179.15, -50.12), ylim = c(20,80), expand = FALSE)

ggplot(data = north_am) +
  geom_sf() +
  geom_sf(data = presencia_sf_YEJU, color = "gray5", alpha= 0.6)+
  geom_sf(data= st_geometry(ca_sf_crop_YEJU), color= "red4") +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-179.15, -50.12), ylim = c(10,70), expand = FALSE)

ggplot(data = north_am) +
  geom_sf() +
  geom_sf(data = presencia_sf_HYEMALIS, color = "gray5", alpha= 0.6)+
  geom_sf(data= st_geometry(ca_sf_crop_HYEMALIS), color= "red4") +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-179.15, -50.12), ylim = c(20,80), expand = FALSE)

##################### Creating background points DEJU ----

# Iterar sobre cada grupo para aplicar el muestreo de background
for (group in group_names_2) {
  # 1. Obtener los datos correspondientes `junco_pf_XXX` y `ca_crop_XXX`
  junco_pf <- get(paste0("junco_pf_", group))
  ca_crop <- get(paste0("ca_sf_crop_", group))  # Área de calibración recortada
  
  # 2. Crear el set de pseudo-ausencias usando `sample_background`
  set.seed(10)
  bg <- flexsdm::sample_background(
    data = junco_pf,
    x = "decimalLongitude",
    y = "decimalLatitude",
    n = nrow(junco_pf) * 10,  # Número de background: 10 veces el número de presencias
    method = "random",
    rlayer = variables$tree,  # Usar la capa raster correspondiente
    calibarea = ca_crop
  )
  
  # 3. Seleccionar las columnas necesarias de `junco_pf` y modificar la columna `pr_ab`
  junco_pf_1 <- junco_pf %>% 
    dplyr::select(-species, id) %>% 
    mutate(pr_ab = 1)
  
  # 4. Combinar presencias y background en un nuevo dataframe `junco_pa`
  junco_pa <- bind_rows(junco_pf_1, bg)
  
  # 5. Asignar el nuevo dataframe a una variable dinámica con el nombre `junco_pa_[NOMBRE]`
  assign(paste0("junco_pa_", group), junco_pa)

  beepr::beep()

}
graphics.off(); dev.list()

dev.off()
beepr::beep()
beepr::beep()
beepr::beep()
png(file="CANDOR.png",    width=600, height=350)
plot(variables_low$tree, col="gray80", legend=FALSE, axes=FALSE)
plot(crop(ca_crop_CANICEPS, variables_low$tree), add=TRUE)
points(presencia_limpia_CANICEPS_2[,c("decimalLongitude", "decimalLatitude")], col = "#00000480")
points(junco_pf_CANICEPS[,c("decimalLongitude", "decimalLatitude")], col = "#5DC86180")

plot(variables_low$tree, col="gray80", legend=FALSE, axes=FALSE)
plot(crop(ca_crop_OREGANUS, variables_low$tree), add=TRUE)
points(presencia_limpia_OREGANUS_2[,c("decimalLongitude", "decimalLatitude")], col = "#00000480")
points(junco_pf_OREGANUS[,c("decimalLongitude", "decimalLatitude")], col = "#5DC86180")

plot(variables_low$tree, col="gray80", legend=FALSE, axes=FALSE)
plot(crop(ca_crop_HYEMALIS, variables_low$tree), add=TRUE)
points(presencia_limpia_HYEMALIS_2[,c("decimalLongitude", "decimalLatitude")], col = "#00000480")
points(junco_pf_HYEMALIS[,c("decimalLongitude", "decimalLatitude")], col = "#5DC86180")

plot(variables_low$tree, col="gray80", legend=FALSE, axes=FALSE)
plot(crop(ca_crop_YEJU, variables_low$tree), add=TRUE)
points(presencia_limpia_YEJU_2[,c("decimalLongitude", "decimalLatitude")], col = "#00000480")
points(junco_pf_YEJU[,c("decimalLongitude", "decimalLatitude")], col = "#5DC86180")

##################### Partition data for evaluating models #####

# Iterar sobre cada grupo para aplicar `part_random` y crear `junco_pa2_[NOMBRE]`
for (group in group_names_2) {
  # 1. Obtener el dataframe `junco_pa_[NOMBRE]`
  junco_pa <- get(paste0("junco_pa_", group))
  
  # 2. Aplicar `part_random` para crear particiones k-fold con replicates
  junco_pa2 <- part_random(
    data = junco_pa,
    pr_ab = "pr_ab",
    method = c(method = "rep_kfold", folds = 5, replicates = 10)  # Usar 5 folds y 10 replicates
  )
  
  # 3. Asignar el nuevo dataframe a una variable dinámica `junco_pa2_[NOMBRE]`
  assign(paste0("junco_pa2_", group), junco_pa2)
}

##################### Extracting environmental values ----

# Iterar sobre cada grupo para aplicar el script de extracción de variables ambientales
for (group in group_names_2) {
  # 1. Obtener el dataframe `junco_pa2_[NOMBRE]` correspondiente
  junco_pa2 <- get(paste0("junco_pa2_", group))
  
  # 2. Extraer variables ambientales para cada punto de presencia/pseudo-ausencia
  presencias_sobre_variables <- extract(
    x = variables_noGLZ,  # Raster de variables ambientales
    y = junco_pa2[ , c("decimalLongitude", "decimalLatitude")]  # Coordenadas de los puntos
  )
  
  # 3. Convertir a dataframe para facilitar la manipulación
  presencias_sobre_variables <- data.frame(presencias_sobre_variables)
  
  # 4. Unir las variables ambientales extraídas con el dataframe `junco_pa2`
  presencia <- cbind(junco_pa2, presencias_sobre_variables)
  
  # 5. Quitar registros con valores nulos en `tree_CRS` (u otra variable de referencia)
  presencia <- presencia[!(is.na(presencia$tree)), ]
  
  # 6. Convertir a `tibble` y seleccionar las columnas deseadas
  pr_junco <- as_tibble(presencia) %>% 
    dplyr::select(-c(ID))  # Eliminar la columna `ID` si existe
  
  # 7. Asignar el nuevo dataframe a una variable dinámica `pr_junco_[NOMBRE]`
  assign(paste0("pr_junco_", group), pr_junco)
}

beepr::beep()
beepr::beep()
beepr::beep()

pr_junco_CANICEPS
pr_junco_OREGANUS
pr_junco_HYEMALIS
pr_junco_YEJU


##################### Multivariate Environmental Similarity Surface  ----

# Convert SpatRasters to RasterStacks for compatibility with the mess() function
variable_sets_MESS <- list(
  "time_present_simple" = as(variables_bio_low, "Raster"),
  "time_21kBP" = as(variables_21kBP, "Raster"),
  "time_16kBP" = as(variables_16kBP, "Raster"),
  "time_12kBP" = as(variables_12kBP, "Raster"),
  "time_8kBP" = as(variables_8kBP, "Raster"),
  "time_4kBP" = as(variables_4kBP, "Raster")
)

names(variable_sets_MESS$time_present)
names(variable_sets_MESS$time_4kBP)
names(variable_sets_MESS$time_8kBP)
names(variable_sets_MESS$time_12kBP)
names(variable_sets_MESS$time_16kBP)
names(variable_sets_MESS$time_21kBP)

group_names_2 <- c("HYEMALIS", "CANICEPS", "OREGANUS", "YEJU")
group_names_2 <- c("YEJU")
group_names_2 <- c("CANICEPS")

for (group in group_names_2) {
  # Get the filtered presences for the calibration set
  calibracion_data <- get(paste0("pr_junco_", group))
  
  # Extract coordinates for the presence points
  coords <- calibracion_data %>% dplyr::select(decimalLongitude, decimalLatitude)
  
  # Use a valid environmental layer for calibration (e.g., time_4kBP)
  calibration_period <- "time_present_simple"  # Update as needed to a valid key in variable_sets_MESS
  if (!calibration_period %in% names(variable_sets_MESS)) {
    stop(paste("Error: Calibration period", calibration_period, "is missing in variable_sets_MESS."))
  }
  
  # Extract environmental values at presence points for calibration data
  env_values <- terra::extract(as(variable_sets_MESS[[calibration_period]], "SpatRaster"), coords)
  
  if (is.null(env_values)) {
    stop(paste("Error: Environmental values could not be extracted for group:", group))
  }
  
  env_values <- env_values[, -1]  # Remove the ID column from `extract`
  
  # Iterate over each environmental variable set to apply the MESS analysis
  for (period_name in names(variable_sets_MESS)) {
    env_layers <- variable_sets_MESS[[period_name]]
    
    # Perform the MESS analysis
    mess_map <- mess(env_layers, env_values)
    
    # Save the result to a dynamic variable
    assign(paste0("mess_", period_name, "_", group), mess_map, envir = .GlobalEnv)
    
    # Confirmation message
    cat(paste("Análisis MESS completado para el grupo:", group, "en el periodo:", period_name, "\n"))
  }
}

## Guardar todos los MESS en un PDF (una página por mapa)
period_keys   <- c("time_present_simple","time_4kBP","time_8kBP","time_12kBP","time_16kBP","time_21kBP")
period_labels <- c("Present","4kya","8kya","12kya","16kya","21kya")
taxa_keys     <- c("HYEMALIS","OREGANUS","CANICEPS","YEJU")
taxa_labels   <- c("hyemalis","oreganus","caniceps","phaeonotus")

if (!dir.exists(Dir.Exports)) dir.create(Dir.Exports, recursive = TRUE)
pdf(file.path(Dir.Exports, "MESS_maps.pdf"), width = 8, height = 6)
par(mar = c(2,2,2,1))

for (i in seq_along(period_keys)) {
  for (j in seq_along(taxa_keys)) {
    nm <- paste0("mess_", period_keys[i], "_", taxa_keys[j])
    if (!exists(nm)) next
    r <- get(nm)
    r <- if (inherits(r, "Raster")) terra::rast(r) else r
    terra::plot(r, legend = TRUE, axes = FALSE, box = FALSE,
                main = paste0(period_labels[i], " — ", taxa_labels[j]))
    # contorno opcional (si tienes 'north_am'):
    try(terra::plot(terra::vect(north_am), add = TRUE, border = "grey25", lwd = 0.6), silent = TRUE)
  }
}
dev.off()

##################### Maxent with vegetation-----
names(variables_noGLZ)
group_names_2 <- c("YEJU")
group_names_2 <- c("HYEMALIS", "CANICEPS", "OREGANUS")
group_names_2 <- c("HYEMALIS", "CANICEPS", "OREGANUS", "YEJU")
# Iterar sobre cada grupo para aplicar el modelo `Maxent`
for (group in group_names_2) {
  # 1. Obtener el dataframe `pr_junco_[NOMBRE]`
  pr_junco <- get(paste0("pr_junco_", group))
  
  # 2. Separar presencias y pseudoausencias
  junco_maxent <- pr_junco %>% 
    filter(pr_ab == 1) %>%  # Filtrar solo presencias
    na.omit()
  
  backg_maxent <- pr_junco %>% 
    filter(pr_ab == 0)  # Filtrar solo pseudoausencias
  
  # 3. Preparar el objeto SWD para Maxent
  data_SWD <- prepareSWD(
    species = paste("Junco", group),
    p = junco_maxent %>% dplyr::select(decimalLongitude, decimalLatitude),
    a = backg_maxent %>% dplyr::select(decimalLongitude, decimalLatitude),
    env = variables_noGLZ  # Raster de variables ambientales (asegúrate de definirlo en tu entorno)
  )
  
  # 4. Crear folds para la partición de los datos
  folds <- SDMtune::randomFolds(
    data_SWD, 
    k = 4, 
    only_presence = TRUE, 
    seed = 25
  )
  
  # 5. Entrenar el modelo Maxent usando `train()`
  model <- train(
    method = "Maxent",
    data = data_SWD,
    folds = folds
  )
  
  # 6. Asignar los objetos generados a nombres dinámicos
  assign(paste0("data_SWD_", group), data_SWD)  # Objeto SWD para cada grupo
  assign(paste0("folds_", group), folds)  # Folds generados
  assign(paste0("maxent_model_", group), model)  # Modelo Maxent ajustado
  beepr::beep()
}

# Iterar sobre cada grupo para calcular la importancia de variables y la validación del modelo
for (group in group_names_2) {
  # 1. Obtener el modelo Maxent `maxent_model_[NOMBRE]`
  model <- get(paste0("maxent_model_", group))
  
  # 2. Calcular la importancia de las variables
  vi <- maxentVarImp(model)
  
  # 3. Mostrar la importancia de las variables
  print(paste("Importancia de las variables para el grupo:", group))
  print(vi)
  # 4. Validación del modelo (Training y Testing AUC)
  cat("\nValidación del modelo para el grupo:", group, "\n")
  cat("Training AUC: ", SDMtune:::auc(model), "\n")
  cat("Testing AUC: ", SDMtune:::auc(model, test = TRUE), "\n\n")
}

# Conjunto de variables ambientales a usar en la predicción (definido en tu entorno)
variables_to_use <- c("bio_10", "bio_15", "bio_18", "bio_2", "bio_4",  "scd", "ndvi", "elev", "tree")

# Iterar sobre cada grupo para aplicar la predicción con el modelo Maxent
for (group in group_names_2) {
  # 1. Obtener el modelo Maxent `maxent_model_[NOMBRE]`
  model <- get(paste0("maxent_model_", group))
  
  # 2. Realizar la predicción en el espacio usando las condiciones ambientales actuales
  map_present_total <- predict(
    model,
    data = variables_low[[variables_to_use]],  # Subconjunto de variables ambientales de `variables_low`
    type = "cloglog"  # Tipo de predicción, puedes cambiar a "logistic" o "raw" si lo prefieres
  )
  
  # 3. Asignar la predicción espacial a una variable dinámica `map_present_total_[NOMBRE]`
  assign(paste0("map_present_total_", group), map_present_total)
  
  # Mensaje de confirmación para cada grupo
  cat(paste("Predicción espacial completada para el grupo:", group, "\n"))
  
  map_df <- as.data.frame(map_present_total, xy = TRUE, na.rm = TRUE)
  # Rename the raster value column to 'suitability' for clarity
  colnames(map_df)[3] <- "suitability"
  assign(paste0("map_df_present_", group), map_df)
}

cl <- c("#FDE725", "#B3DC2B", "#6DCC57", "#36B677", "#1F9D87", "#25818E", "#30678D", "#3D4988", "#462777", "#440154")


plotPred(map_present_total_HYEMALIS,
         lt = "Habitat\nsuitability - HYEMALIS",
         colorramp = cl)

plotPred(map_present_total_CANICEPS,
         lt = "Habitat\nsuitability - CANICEPS",
         colorramp = cl)

plotPred(map_present_total_OREGANUS,
         lt = "Habitat\nsuitability - OREGANUS",
         colorramp = cl)

plotPred(map_present_total_YEJU,
         lt = "Habitat\nsuitability - YEJU",
         colorramp = cl)

# Rename the raster value column to 'suitability' for clarity
colnames(map_df)[3] <- "suitability"
head(map_df)

ggplot() +
  geom_raster(
    data = map_df,
    aes(x = x, y = y, fill = map_df[[3]])
  ) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    na.value = "transparent",
    guide = guide_colorbar(
      direction = "horizontal",    # Barra horizontal
      barwidth = unit(10, "cm"),   # Ajusta el ancho
      barheight = unit(0.5, "cm"), # Ajusta la altura
      title.position = "top",      # Título arriba de la barra
      title.hjust = 0.5,           # Centrar título
      frame.colour = "black",      # Color del recuadro
      frame.linewidth = 0.1        # Grosor de la línea del recuadro
    )
  ) +
  coord_sf() +
  theme_minimal() +
  geom_sf(data = elev_df_present_simple, fill = NA, color = "black", size = 0.5) +
  geom_sf(data = glz_df_present_simple, fill = "#A8D6FA", color = NA, alpha = 1) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10)
  )

  
###### PLOT AND SAVE WITH MESS MASK AND GLACIER
beep()
beep()
beep()
beep()
beep()

output_dir <- "Figures10.0/"  # Specify the directory to save the maps
dir.create(output_dir, showWarnings = FALSE)  # Create directory if it doesn't exist
# Loop over the group names
group_names <- c("HYEMALIS", "CANICEPS", "YEJU", "OREGANUS")
group_names <- c("YEJU")
group_names <- c("CANICEPS")

for (group in group_names) {
  # Construct the name of the raster
  raster_name <- paste0("map_present_total_", group)
  
  # Check if the raster exists
  if (exists(raster_name, envir = .GlobalEnv)) {
    # Retrieve the raster
    pred_raster <- get(raster_name, envir = .GlobalEnv)
    
    # Convert the raster to a data frame
    map_df <- as.data.frame(pred_raster, xy = TRUE, na.rm = TRUE)
    
    # Rename the third column to 'suitability'
    if (ncol(map_df) >= 3) {
      colnames(map_df)[3] <- "suitability"
    } else {
      stop(paste("Error: Raster data frame has fewer than 3 columns for", raster_name))
    }
    
    # Create the plot using your provided ggplot code
    plot <- ggplot() +
      geom_raster(data = map_df, aes(x = x, y = y, fill = suitability)) +
      scale_fill_gradient(
        low = "white", high = "darkred",
        na.value = "transparent"
      ) +
      coord_sf() +
      theme_minimal() +
      geom_sf(data = elev_df_present_simple, fill = NA, color = "black", size = 0.5) +
      geom_sf(data = glz_df_present_simple, fill = "#A8D6FA", color = NA, alpha = 1) +
      theme(
        panel.grid = element_blank(),        # Remove the grid
        axis.text = element_blank(),         # Remove x and y axis text
        axis.title = element_blank(),        # Remove x and y axis labels
        axis.ticks = element_blank(),        # Remove axis ticks
        legend.position = "none"             # Remove the legend
      )
    
    # Save the plot to a file
    output_file <- file.path(output_dir, paste0("map_present_allVARS_", group, ".png"))
    ggsave(filename = output_file, plot = plot, width = 8, height = 6, dpi = 300)
    
    # Confirmation message
    cat(paste("Plot saved for group:", group, "->", output_file, "\n"))
    
  } else {
    cat(paste("The raster", raster_name, "does not exist in the global environment.\n"))
  }
}

##################### Test TSS n' boyce -----
# Inicializamos listas para almacenar los umbrales óptimos, TSS y Boyce Index
optimal_thresholds <- list()
TSS_results <- list()
boyce_indices <- list()

# Supongamos que 'group_names_2' contiene los grupos que quieres evaluar, por ejemplo:
# group_names_2 <- c("HYEMALIS", "CANICEPS", "OREGANUS", "YEJU")

for (group in group_names_2) {
  # Obtener datos observados y predicción del modelo
  pr_junco <- get(paste0("pr_junco_", group))
  map_present_total <- get(paste0("map_present_total_", group))
  
  # Extraer los valores de predicción en las ubicaciones de presencia
  predicted_values <- terra::extract(map_present_total, pr_junco[, c("decimalLongitude", "decimalLatitude")])
  pr_junco$predicted <- predicted_values$mean
  pr_junco <- pr_junco %>% filter(!is.na(predicted))
  
  # Preparar datos para calcular TSS
  data_PA <- data.frame(
    id = 1:nrow(pr_junco),
    observed = pr_junco$pr_ab,
    predicted = pr_junco$predicted
  )
  
  thresholds_seq <- seq(0, 1, by = 0.001)
  eval_results <- PresenceAbsence::optimal.thresholds(
    DATA = data_PA,
    threshold = thresholds_seq,
    opt.methods = c("MaxSens+Spec")
  )
  optimal_threshold <- eval_results$predicted[1]
  optimal_thresholds[[group]] <- optimal_threshold
  
  cm <- PresenceAbsence::cmx(
    DATA = data_PA,
    threshold = optimal_threshold
  )
  
  sensitivity <- PresenceAbsence::sensitivity(cm, st.dev = FALSE)
  specificity <- PresenceAbsence::specificity(cm, st.dev = FALSE)
  TSS <- sensitivity + specificity - 1
  TSS_results[[group]] <- TSS
  
  cat("Grupo:", group, "\n")
  cat("Umbral Óptimo:", round(optimal_threshold, 3), "\n")
  cat("TSS:", round(TSS, 3), "\n")
  
  # ---- Cálculo del índice de Boyce ----
  # Extraer todos los valores de predicción del raster (excluyendo NA)
  all_pred <- terra::values(map_present_total)
  all_pred <- all_pred[!is.na(all_pred)]
  
  # Valores de predicción en las ubicaciones de presencia
  pres_pred <- pr_junco$predicted
  
  # Calcular el índice de Boyce usando ecospat.boyce
  # Ajusta 'window.w' según convenga (por ejemplo, 0.2 puede ser un buen punto de partida)
  boyce_result <- ecospat::ecospat.boyce(
    fit = all_pred,
    obs = pres_pred,
    nclass = 0,         # Let the function determine the number of classes
    window.w = 0.2,     # Adjust this parameter as needed
    PEplot = FALSE      # Set to TRUE if you want to see the PE plot
  )
  
  boyce_index <- boyce_result$cor
  boyce_indices[[group]] <- boyce_index
  
  cat("Índice de Boyce:", round(boyce_index, 3), "\n\n")
}

##################### Predicción de modelo simple MAXENT -----

# Preparar y entrenar el modelo Maxent para cada grupo usando el subconjunto de variables
for (group in group_names_2) {
  # Obtener el dataframe `pr_junco_[NOMBRE]` correspondiente al grupo
  pr_junco <- get(paste0("pr_junco_", group))  # Asegúrate de que `pr_junco_[NOMBRE]` esté definido en el entorno
  
  # Separar presencias y pseudoausencias
  junco_maxent <- pr_junco %>%
    filter(pr_ab == 1) %>%
    na.omit()  # Filtrar solo presencias y eliminar `NA`s
  
  backg_maxent <- pr_junco %>%
    filter(pr_ab == 0)  # Filtrar solo pseudoausencias
  
  # Preparar el objeto SWD con el subconjunto de variables `variables_bio` como antes
  data_SWD_simple <- prepareSWD(
    species = paste("Junco", group),  # Usar el nombre del grupo para la especie
    p = junco_maxent %>% dplyr::select(decimalLongitude, decimalLatitude),
    a = backg_maxent %>% dplyr::select(decimalLongitude, decimalLatitude),
    env = variables_bio_low  # Subconjunto de variables ambientales
  )
  
  # Crear folds para la partición de los datos usando `randomFolds`
  folds_simple <- randomFolds(
    data_SWD_simple, 
    k = 4, 
    only_presence = TRUE, 
    seed = 25
  )
  
  # Entrenar el modelo Maxent usando `train()` con el subconjunto de variables
  model_simple <- train(
    method = "Maxent",
    data = data_SWD_simple,
    folds = folds_simple
  )
  
  # Asignar el objeto SWD y el modelo a nombres dinámicos
  assign(paste0("data_SWD_simple_", group), data_SWD_simple)  # Objeto SWD con variables bio
  assign(paste0("maxent_model_simple_", group), model_simple)  # Modelo ajustado con variables bio
  
  # Mensaje de confirmación para cada grupo
  cat(paste("Modelo Maxent con variables bio completado para el grupo:", group, "\n"))
}

#group_names <- c("HYEMALIS", "CANICEPS", "YEJU", "OREGANUS")
#group_names <- c("YEJU")
#group_names <- c("HYEMALIS")

# Iterar sobre cada grupo para calcular la importancia de variables y la validación del modelo
for (group in group_names_2) {
  # 1. Obtener el modelo Maxent `maxent_model_[NOMBRE]`
  model <- get(paste0("maxent_model_simple_", group))
  
  # 2. Calcular la importancia de las variables
  vi <- maxentVarImp(model)
  
  # 3. Mostrar la importancia de las variables
  print(paste("Importancia de las variables para el grupo:", group))
  print(vi)
  # 4. Validación del modelo (Training y Testing AUC)
  cat("\nValidación del modelo para el grupo:", group, "\n")
  cat("Training AUC: ", SDMtune:::auc(model), "\n")
  cat("Testing AUC: ", SDMtune:::auc(model, test = TRUE), "\n\n")
}


# Realizar la predicción en el espacio usando `variables_bio_low` para cada grupo y asignar a nombres dinámicos
for (group in group_names_2) {
  # Obtener el modelo entrenado para cada grupo
  model_simple <- get(paste0("maxent_model_simple_", group))
  
  # Realizar la predicción en el espacio usando `variables_bio_low`
  map_present_simple <- predict(
    model_simple,
    data = variables_bio_low,  # Subconjunto de variables ambientales `variables_bio_low`
    type = "cloglog"  # Tipo de predicción (puedes cambiar a "logistic" o "raw" si lo prefieres)
  )
  
  # Asignar la predicción a un nombre dinámico
  assign(paste0("map_present_simple_", group), map_present_simple)  # Predicción espacial con variables bio
  
  # Mensaje de confirmación para cada grupo
  cat(paste("Predicción espacial con variables bio completada para el grupo:", group, "\n"))
  
  map_df <- as.data.frame(map_present_simple, xy = TRUE, na.rm = TRUE)
  # Rename the raster value column to 'suitability' for clarity
  colnames(map_df)[3] <- "suitability"
  assign(paste0("map_df_present_simple_", group), map_df)
}

################# TSS SIMPLE 
# He movido este arriba para podeer maskear suitability con MTSS
group_names_2 <- c("HYEMALIS", "CANICEPS", "YEJU", "OREGANUS")

# Create a list to store the optimal thresholds for each group (simple model)
optimal_thresholds_simple <- list()

for (group in group_names_2) {
  # Retrieve the observed data for the group
  pr_junco <- get(paste0("pr_junco_", group))
  
  # Retrieve the predicted spatial raster for the group (simple model)
  map_present_simple <- get(paste0("map_present_simple_", group))
  
  # Extract predicted suitability values at observation points
  predicted_values <- terra::extract(map_present_simple, pr_junco[, c("decimalLongitude", "decimalLatitude")])
  
  # Use the 'mean' column for predicted suitability
  pr_junco$predicted <- predicted_values$mean
  
  # Remove rows with NA predicted values
  pr_junco <- pr_junco %>% filter(!is.na(predicted))
  
  # Prepare data for TSS calculation
  data_PA <- data.frame(
    id = 1:nrow(pr_junco),
    observed = pr_junco$pr_ab,
    predicted = pr_junco$predicted
  )
  
  # Generate a sequence of thresholds from 0 to 1
  thresholds_seq <- seq(0, 1, by = 0.001)
  
  # Use PresenceAbsence package to find the optimal threshold for TSS (simple model)
  eval_results_simple <- PresenceAbsence::optimal.thresholds(
    DATA = data_PA,
    threshold = thresholds_seq,
    opt.methods = c("MaxSens+Spec")
  )
  
  # Extract the optimal threshold for the simple model
  optimal_threshold_simple <- eval_results_simple$predicted[1]
  
  # Save the optimal threshold in the list
  optimal_thresholds_simple[[group]] <- optimal_threshold_simple
  
  # Create a confusion matrix using the optimal threshold
  cm_simple <- PresenceAbsence::cmx(
    DATA = data_PA,
    threshold = optimal_threshold_simple
  )
  
  # Calculate Sensitivity and Specificity
  sensitivity_simple <- PresenceAbsence::sensitivity(cm_simple, st.dev = FALSE)
  specificity_simple <- PresenceAbsence::specificity(cm_simple, st.dev = FALSE)
  
  # Calculate TSS for the simple model
  TSS_simple <- sensitivity_simple + specificity_simple - 1
  
  # Print results for each group
  cat("Group:", group, "\n")
  cat("Optimal Threshold (Simple Model):", round(optimal_threshold_simple, 3), "\n")
  cat("TSS (Simple Model):", round(TSS_simple, 3), "\n")
  
  # ---- Cálculo del índice de Boyce para el modelo simple ----
  # Extraer todos los valores de predicción del raster (excluyendo NA)
  all_pred_simple <- terra::values(map_present_simple)
  all_pred_simple <- all_pred_simple[!is.na(all_pred_simple)]
  
  # Valores de predicción en las ubicaciones de presencia
  pres_pred_simple <- pr_junco$predicted
  
  # Calcular el índice de Boyce usando ecospat.boyce
  boyce_result_simple <- ecospat::ecospat.boyce(
    fit = all_pred_simple,
    obs = pres_pred_simple,
    nclass = 0,         # Deja que la función determine el número de clases
    window.w = 0.2,     # Ajusta este parámetro según la distribución de tus predicciones
    PEplot = FALSE      # Cambia a TRUE si deseas visualizar el gráfico de PE
  )
  
  boyce_index_simple <- boyce_result_simple$cor
  cat("Boyce Index (Simple Model):", round(boyce_index_simple, 3), "\n\n")
  
  # Opcional: guardar los valores para análisis posterior
  assign(paste0("TSS_simple_", group), TSS_simple)
  assign(paste0("optimal_threshold_simple_", group), optimal_threshold_simple)
  assign(paste0("boyce_index_simple_", group), boyce_index_simple)
}

cat("
Group: HYEMALIS 
Optimal Threshold (Simple Model): 0.466 
TSS (Simple Model): 0.496 

Group: CANICEPS 
Optimal Threshold (Simple Model): 0.264 
TSS (Simple Model): 0.74 

Group: OREGANUS 
Optimal Threshold (Simple Model): 0.288 
TSS (Simple Model): 0.593 

Group: YEJU con buffer 150km
Optimal Threshold (Simple Model): 0.16 
TSS (Simple Model): 0.758 
    ")


###### PLOT AND SAVE WITH MESS MASK AND GLACIER SIMPLE

output_dir <- "/media/milab/Cronus_SSD/Q_niche_modeling/Figures10.0"  # Specify the directory to save the maps
dir.create(output_dir, showWarnings = TRUE)  # Create directory if it doesn't exist
# Loop over the group names
group_names <- c("HYEMALIS", "CANICEPS", "YEJU", "OREGANUS")

for (group in group_names) {
  # Construct the name of the raster
  raster_name <- paste0("map_present_simple_", group)
  
  # Check if the raster exists
  if (exists(raster_name, envir = .GlobalEnv)) {
    # Retrieve the raster
    pred_raster <- get(raster_name, envir = .GlobalEnv)
    
    # Convert the raster to a data frame
    map_df <- as.data.frame(pred_raster, xy = TRUE, na.rm = TRUE)
    
    # Rename the third column to 'suitability'
    if (ncol(map_df) >= 3) {
      colnames(map_df)[3] <- "suitability"
    } else {
      stop(paste("Error: Raster data frame has fewer than 3 columns for", raster_name))
    }
    
    # Create the plot using your provided ggplot code
    plot <- ggplot() +
      geom_raster(data = map_df, aes(x = x, y = y, fill = suitability)) +
      scale_fill_gradient(
        low = "white", high = "darkred",
        na.value = "transparent"
      ) +
      coord_sf() +
      theme_minimal() +
      geom_sf(data = elev_df_present_simple, fill = NA, color = "black", size = 0.5) +
      geom_sf(data = glz_df_present_simple, fill = "#A8D6FA", color = NA, alpha = 1) +
      theme(
        panel.grid = element_blank(),        # Remove the grid
        axis.text = element_blank(),         # Remove x and y axis text
        axis.title = element_blank(),        # Remove x and y axis labels
        axis.ticks = element_blank(),        # Remove axis ticks
        legend.position = "none"             # Remove the legend
      )
    
    # Save the plot to a file
    output_file <- file.path(output_dir, paste0("map_present_simple_", group, ".pdf"))
    ggsave(filename = output_file, plot = plot, width = 8, height = 6, device = "pdf")

    # Confirmation message
    cat(paste("Plot saved for group:", group, "->", output_file, "\n"))
    
  } else {
    cat(paste("The raster", raster_name, "does not exist in the global environment.\n"))
  }
}

###### PLOT AND SAVE WITH MESS MASK AND GLACIER SIMPLE and water(FILTERED BY MTSS)

# -------------------------------
# PREPARE THE WATER BODIES VECTOR
# -------------------------------

# Load the lakes shapefile using terra's vect()
water_bodies <- vect("/media/milab/Cronus_SSD/Q_niche_modeling/data/ne_10m_lakes/ne_10m_lakes.shp")

# Validate geometries & repair invalid ones
water_bodies_valid <- makeValid(water_bodies)

# Convert to a Spatial object for compatibility with raster package functions
water_sp <- as(water_bodies_valid, "Spatial")

# (Optional) Reproject water_sp to match your prediction raster CRS.
# Here we extract the CRS from a known prediction raster:
pred_temp <- get("map_present_simple_HYEMALIS", envir = .GlobalEnv)
# Convert to RasterLayer if needed so that CRS extraction works:
pred_temp <- raster(pred_temp)
pred_crs <- proj4string(pred_temp)
water_sp <- spTransform(water_sp, CRS(pred_crs))
water_vect <- water_sp
# -------------------------------
# SET OUTPUT DIRECTORY
# -------------------------------
output_dir <- "Figures10.0/"  # Specify the directory to save the maps
dir.create(output_dir, showWarnings = FALSE)  # Create directory if it doesn't exist

# -------------------------------
# LOOP OVER GROUPS AND APPLY WATER MASK
# -------------------------------
group_names <- c("HYEMALIS", "CANICEPS", "YEJU", "OREGANUS")

for (group in group_names) {
  # Construct the name of the raster
  raster_name <- paste0("map_present_simple_", group)
  
  # Check if the raster exists in the global environment
  if (exists(raster_name, envir = .GlobalEnv)) {
    # Retrieve the prediction raster and convert it to a RasterLayer for compatibility.
    pred_raster <- raster(get(raster_name, envir = .GlobalEnv))
    
    # -------------------------------
    # APPLY WATER MASK
    # -------------------------------
    # Crop water bodies to the extent of the prediction raster
    water_clip <- crop(water_sp, extent(pred_raster))
    
    # Rasterize the cropped water bodies. We use 'field = 1' (adjust if needed).
    water_mask <- raster::rasterize(water_clip, pred_raster, field = 1)
    
    # Set cells to NA in pred_raster wherever water_mask is not NA
    pred_raster[!is.na(water_mask)] <- NA
    
    # -------------------------------
    # CONVERT THE RASTER TO A DATA FRAME
    # -------------------------------
    map_df <- as.data.frame(pred_raster, xy = TRUE, na.rm = TRUE)
    
    # Rename the third column to 'suitability'
    if (ncol(map_df) >= 3) {
      colnames(map_df)[3] <- "suitability"
    } else {
      stop(paste("Error: Raster data frame has fewer than 3 columns for", raster_name))
    }
    
    # Retrieve the MTSS threshold for the current group
    threshold_simple <- optimal_thresholds_simple[[group]]
    
    # Apply the MTSS threshold filter: Set suitability values below the threshold to NA (mask them)
    map_df$suitability[map_df$suitability < threshold_simple] <- NA
    
    # -------------------------------
    # CREATE AND SAVE THE PLOT
    # -------------------------------
    plot <- ggplot() +
      geom_raster(data = map_df, aes(x = x, y = y, fill = suitability)) +
      scale_fill_gradient(
        low = "white", high = "darkred",
        na.value = "transparent"  # Masked areas will be transparent
      ) +
      coord_sf() +
      theme_minimal() +
      geom_sf(data = elev_df_present_simple, fill = NA, color = "black", size = 0.5) +
      geom_sf(data = glz_df_present_simple, fill = "#A8D6FA", color = NA, alpha = 1) +
      theme(
        panel.grid = element_blank(),        # Remove the grid
        axis.text = element_blank(),           # Remove x and y axis text
        axis.title = element_blank(),          # Remove x and y axis labels
        axis.ticks = element_blank(),          # Remove axis ticks
        legend.position = "none"               # Remove the legend
      )
    
    # Save the plot to a file
    output_file <- file.path(output_dir, paste0("map_present_simple_filtered_", group, ".pdf"))
    ggsave(filename = output_file, plot = plot, width = 8, height = 6, device = "pdf")
    
    # Confirmation message
    cat(paste("Filtered plot saved for group:", group, "->", output_file, "\n"))
    
  } else {
    cat(paste("The raster", raster_name, "does not exist in the global environment.\n"))
  }
}

##################### Predicción de modelo pasado MAXENT -----

# Definir los conjuntos de variables ambientales
variable_sets <- list(
  "21kBP" = variables_21kBP,
  "4kBP" = variables_4kBP,
  "8kBP" = variables_8kBP,
  "12kBP" = variables_12kBP,
  "16kBP" = variables_16kBP
)

variable_sets <- list(
  "21kBP" = variables_4kBP,
  "8kBP" = variables_8kBP,
  "12kBP" = variables_12kBP
)

variable_sets <- list(
  "4kBP" = variables_4kBP,
  "16kBP" = variables_16kBP
)

group_names_2 <- c("HYEMALIS", "CANICEPS", "OREGANUS", "YEJU")
#group_names_2 <- c("YEJU")
group_names_2 <- c("HYEMALIS")

# Función `pred` para realizar predicciones con diferentes conjuntos de variables
pred <- function(group, variable_set_name, variable_set) {
  # Obtener el modelo simple del grupo
  model <- get(paste0("maxent_model_simple_", group))
  
  # Realizar la predicción usando el conjunto de variables proporcionado
  pred_map <- predict(
    model,
    data = variable_set,  # Conjunto de variables (variables_21kBP, variables_pleis12ka, o variables_midH)
    type = "cloglog"  # Tipo de predicción (puedes cambiar a "logistic" o "raw" si lo prefieres)
  )
  
  # Asignar la predicción a un nombre dinámico
  assign(paste0("map_", variable_set_name, "_", group), pred_map, envir = .GlobalEnv)  # Guardar en el entorno global
  
  # Mensaje de confirmación
  cat(paste("Predicción completada para el grupo", group, "usando el conjunto de variables", variable_set_name, "\n"))
}

# Aplicar la función `pred()` para cada grupo y cada conjuinto de variables
for (group in group_names_2) {
  for (variable_set_name in names(variable_sets)) {
    pred(group, variable_set_name, variable_sets[[variable_set_name]])
  }
}



##################### creation of map_DF to plot in ggplot -----

#optimal_thresholds_simple$HYEMALIS <- 0.296

optimal_thresholds_simple$HYEMALIS
optimal_thresholds_simple$CANICEPS
optimal_thresholds_simple$DORSALIS
optimal_thresholds_simple$YEJU

# Load the lakes shapefile using terra's vect()
water_bodies <- vect("/media/milab/Cronus_SSD/Q_niche_modeling/data/ne_10m_lakes/ne_10m_lakes.shp")

# Validate geometries & repair invalid ones
water_bodies_valid <- makeValid(water_bodies)

# Convert to a Spatial object for compatibility with raster package functions
water_sp <- as(water_bodies_valid, "Spatial")
water_vect <- vect(water_sp)

time_periods <- c("present_simple", "4kBP", "8kBP", "12kBP", "16kBP", "21kBP")
group_names <- c("HYEMALIS", "CANICEPS", "OREGANUS", "YEJU")

time_periods <- c("present_simple")
group_names <- c("HYEMALIS")

group_names <- c("YEJU")

for (group in group_names) {
  for (time in time_periods) {
    # Define the names of the prediction and MESS raster objects dynamically
    if (time == "present_simple") {
      pred_name <- paste0("map_present_simple_", group)
      mess_name <- paste0("mess_time_present_simple_", group)
    } else {
      pred_name <- paste0("map_", time, "_", group)
      mess_name <- paste0("mess_time_", time, "_", group)
    }
    
    # Check if the prediction raster object exists
    if (exists(pred_name, envir = .GlobalEnv)) {
      # Retrieve the prediction raster object (this version works fine)
      pred_raster <- get(pred_name, envir = .GlobalEnv)
      
      # -------------------------------
      # (OPTIONAL) APPLY WATER MASK
      # -------------------------------
      # If you want to integrate a water mask later, you can use a block like this.
      # Make sure that your water bodies are a terra SpatVector (e.g., water_vect)
      # and that they have been reprojected to match pred_raster.
      #
      # Uncomment and adjust the following block step by step if desired:
      #
       {
         water_clip <- terra::crop(water_vect, terra::ext(pred_raster))
         if (!is.null(water_clip) && length(water_clip) > 0) {
           water_mask <- terra::rasterize(water_clip, pred_raster, field = 1)
           # Replace any NaN values with NA
           water_mask[ is.nan(water_mask) ] <- NA
           if (!all(is.na(terra::values(water_mask)))) {
             pred_raster[ !is.na(water_mask) ] <- NA
           } else {
             cat("Warning: water_mask for", pred_name, "has no valid values. Skipping water masking.\n")
           }
         } else {
           cat("Warning: No water features intersect", pred_name, "\n")
         }
       }
      
      # -------------------------------
      # APPLY MESS MASK
      # -------------------------------
      if (exists(mess_name, envir = .GlobalEnv)) {
        mess_raster <- get(mess_name, envir = .GlobalEnv)
        # If mess_raster is a RasterLayer, convert it to a terra SpatRaster
        if (inherits(mess_raster, "RasterLayer")) {
          mess_raster <- terra::rast(mess_raster)
        }
        # Ensure matching extent/resolution; resample if necessary
        if (!terra::compareGeom(pred_raster, mess_raster, stopOnError = FALSE)) {
          mess_raster <- terra::resample(mess_raster, pred_raster)
        }
        # Calculate global statistics and determine threshold
        mess_stats <- terra::global(mess_raster, fun = c('min', 'max'), na.rm = TRUE)
        threshold <- (mess_stats$min) / 3
        # Apply the mask: set cells in pred_raster to NA where mess_raster < threshold
        pred_raster[ mess_raster < threshold ] <- NA
      } else {
        cat(paste("The MESS raster object", mess_name, "does not exist in the global environment.\n"))
        next
      }
      

      # -------------------------------
      # APPLY MTSS THRESHOLD 
      # -------------------------------
      
      ## Retrieve the MTSS threshold for the current group
      #if (exists("optimal_thresholds_simple", envir = .GlobalEnv)) {
      #  threshold_mtss <- optimal_thresholds_simple[[group]]
      #  pred_raster[ pred_raster < threshold_mtss ] <- NA
      #}
      #
      # -------------------------------
      # CONVERT MASKED RASTER TO DATA FRAME
      # -------------------------------
      map_df <- as.data.frame(pred_raster, xy = TRUE, na.rm = TRUE)
      if (ncol(map_df) >= 3) {
        colnames(map_df)[3] <- "suitability"
      } else {
        stop(paste("Error: Raster data frame has fewer than 3 columns for", pred_name))
      }
      
      # Save the data frame with a dynamic name in the global environment
      assign(paste0("map_df_", time, "_", group), map_df, envir = .GlobalEnv)
      
    } else {
      cat(paste("The raster object", pred_name, "does not exist in the global environment.\n"))
    }
  }
}


################## plotting all ----
# -------------------------------
# SET OUTPUT DIRECTORY
# -------------------------------
output_dir <- "Figures10.0/"  # Specify the directory to save the maps
dir.create(output_dir, showWarnings = FALSE)  # Create directory if it doesn't exist

# Definir los nombres de los grupos y periodos de tiempo con la capitalización correcta
group_names <- c("HYEMALIS", "CANICEPS", "OREGANUS", "YEJU")

time_periods <- c("present_simple", "4kBP", "8kBP","12kBP", "16kBP", "21kBP")  # Nombres con minúsculas para MidH y Pleis12ka


# Create a container to store all the plots
plots <- list()

for (group in group_names) {
  for (time in time_periods) {
    # Define the prediction name
    if (time == "present_simple") {
      pred_name <- paste0("map_present_simple_", group)
    } else {
      pred_name <- paste0("map_", time, "_", group)
    }
    
    # Check if the prediction exists
    if (exists(pred_name, envir = .GlobalEnv)) {
      pred_raster <- get(pred_name, envir = .GlobalEnv)
      
      # Define the associated data frame name
      map_df_name <- paste0("map_df_", time, "_", group)
      
      # Check if the data frame exists
      if (exists(map_df_name, envir = .GlobalEnv)) {
        map_df <- get(map_df_name, envir = .GlobalEnv)
        
        # Debug: Print map_df structure
        cat("Debugging map_df structure for:", map_df_name, "\n")
        str(map_df)
        
        if (ncol(map_df) < 3) {
          stop(paste("Error: Less than 3 columns in", map_df_name))
        }
        
        # Dynamically select the 3rd column for plotting
        map_df$selected_value <- map_df[[3]]
        
        # Define elevation and glacier datasets dynamically
        elev_data_name <- paste0("elev_df_", time)
        glz_data_name <- paste0("glz_df_", time)
        
        # Check if the elevation and glacier datasets exist
        if (!exists(elev_data_name, envir = .GlobalEnv)) {
          stop(paste("Error: Elevation dataset", elev_data_name, "does not exist"))
        }
        if (!exists(glz_data_name, envir = .GlobalEnv)) {
          stop(paste("Error: Glacier dataset", glz_data_name, "does not exist"))
        }
        
        # Retrieve the elevation and glacier datasets
        elev_data <- get(elev_data_name, envir = .GlobalEnv)
        glz_data <- get(glz_data_name, envir = .GlobalEnv)
        
        # Create the plot
        plot <- ggplot() +
          geom_raster(data = map_df, aes(x = x, y = y, fill = selected_value)) +
          scale_fill_gradient(
            name = "Selected Value",
            low = "white", high = "darkred",
            na.value = "transparent"
          ) +
          theme_minimal() +
          geom_sf(data = elev_data, fill = NA, color = "black", size = 0.5) +
          geom_sf(data = glz_data, fill = "#A8D6FA", color = NA, alpha = 1) +
          theme(
            panel.grid = element_blank(),        # Remove the grid
            axis.text = element_blank(),         # Remove x and y axis text
            axis.title = element_blank(),        # Remove x and y axis labels
            axis.ticks = element_blank(),        # Remove axis ticks
            legend.position = "none"             # Remove the legend
          )
        
        # Store the plot
        plots[[paste(group, time, sep = "_")]] <- plot
      } else {
        cat(paste("The object", map_df_name, "does not exist in the global environment.\n"))
      }
    } else {
      cat(paste("The object", pred_name, "does not exist in the global environment.\n"))
    }
  }
}

# Get the total number of plots
length(plots)

plots$OREGANUS_present_simple
plots$OREGANUS_8kBP
plots$OREGANUS_16kBP
plots$OREGANUS_21kBP

plots$CANICEPS_present_simple
plots$CANICEPS_4kBP
plots$CANICEPS_8kBP
plots$CANICEPS_16kBP
plots$CANICEPS_21kBP

plots$HYEMALIS_present_simple
plots$HYEMALIS_8kBP
plots$HYEMALIS_16kBP
plots$HYEMALIS_21kBP

plots$YEJU_present_simple
plots$YEJU_8kBP
plots$YEJU_16kBP
plots$YEJU_21kBP

######### Transpose and plot in PDF (not recomendable too big to edit)----
# it is better to plot the grid in TIFF and plot each map individualy in a vectorized format so t

# Create the new order explicitly by transposing the hierarchy
new_order <- c()
for (time in time_periods) {
  for (group in group_names) {
    new_order <- c(new_order, paste(group, time, sep = "_"))
  }
}

# Reorder the plots list based on the new order
plots <- plots[new_order]

# Verify the new order
names(plots)

# Define the output PDF filename and properties
pdf(file = "Figures7.0/habitat_suitability_predictions_no_labels_messTHIRD_mask_water.pdf",
    width = 10,       # Width in inches
    height = 14,      # Height in inches
    useDingbats = FALSE)  # Avoids font issues in Adobe Illustrator

# Create the grid of plots, ensuring a tight layout
plot_grid <- arrangeGrob(
  grobs = plots,      # List of reordered plots
  nrow = 6,          # Set rows
  ncol = 4,          # Set columns
  top = NULL,        # Remove top labels
  bottom = NULL,     # Remove bottom labels
  left = NULL,       # Remove left labels
  right = NULL,      # Remove right labels
  padding = unit(0, "line"),  # Minimize padding
  heights = unit(rep(1, 6), "null")  # Ensure rows take up equal space with minimal gaps
)

# Draw the grid of plots
grid.newpage()
grid.draw(plot_grid)

# Close the PDF device to save the file
dev.off()


######### Transpose and plot in TIFF ----

# Definir la carpeta de salida
output_folder_tiff <- "Figures10.0/"

# Asegurar que la carpeta existe
if (!dir.exists(output_folder_tiff)) {
  dir.create(output_folder_tiff, recursive = TRUE)
}

# Definir el nombre del archivo TIFF
tiff_filename <- paste0(output_folder_tiff, "habitat_suitability_predictions_no_labels_messTHIRD_mask_water_mask.tiff")

# Crear el nuevo orden transpuesto
new_order <- c()
for (time in time_periods) {
  for (group in group_names) {
    new_order <- c(new_order, paste(group, time, sep = "_"))
  }
}

# Reordenar la lista de gráficos en base al nuevo orden
plots <- plots[new_order]

# Abrir el dispositivo gráfico para guardar en TIFF
tiff(filename = tiff_filename, width = 10, height = 14, units = "in", res = 300, compression = "lzw")

# Crear la cuadrícula de gráficos asegurando una disposición compacta
plot_grid <- arrangeGrob(
  grobs = plots, 
  nrow = 6, 
  ncol = 4,
  top = NULL,
  bottom = NULL,
  left = NULL,
  right = NULL,
  padding = unit(0, "line"),  # Minimizar espacios
  heights = unit(rep(1, 6), "null")  # Asegurar filas iguales
)

# Dibujar la cuadrícula en el TIFF
grid.newpage()
grid.draw(plot_grid)

# Cerrar el dispositivo gráfico para guardar el archivo
dev.off()

cat("✅ Se ha guardado la cuadrícula completa en formato TIFF en:", tiff_filename, "\n")



##### PLOT ALL in PDF individually to allow the file to be manegable

# Definir la carpeta de salida
output_folder <- "Figures9.0/individual_plots/"

# Crear la carpeta si no existe
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Guardar cada gráfico en un archivo PDF individual
for (plot_name in names(plots)) {
  # Definir el nombre del archivo con prefijo "PLOTS_"
  file_name <- paste0(output_folder, "PLOTS_", plot_name, ".pdf")
  
  # Guardar el plot en PDF con Cairo para optimización
  cairo_pdf(file = file_name, width = 8, height = 6, family = "Arial", pointsize = 10)
  print(plots[[plot_name]])
  dev.off()
}

cat("✅ Todos los gráficos han sido guardados en formato PDF en la carpeta:", output_folder, "\n")


beep()
beep()
beep()
beep()
