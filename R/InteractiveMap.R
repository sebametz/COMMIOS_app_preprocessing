# return a leaflet object based on current context database
interactiveMap <- function(df, Lat, Long, label) {
  
  if(nrow(df) == 0){
    map <- leaflet() %>%
      addProviderTiles(providers$Esri.WorldTopoMap) %>%
      setMaxBounds(lng1 = -15, lat1 = 35, lng2 = 15, lat2 = 65)
    return(map)
  } 
  
  # Convert Lat, Long, and label to strings if they are not already
  df_n <- df |>
    mutate(Long = .data[[Long]], Lat = .data[[Lat]], label = .data[[label]])
  
  # Create the map with markers
  map <- leaflet(df_n) %>%
    addProviderTiles(providers$Esri.WorldTopoMap) %>%
    setMaxBounds(lng1 = -15, lat1 = 35, lng2 = 15, lat2 = 65) %>%
    addMarkers(lng = ~ Long, lat = ~ Lat,
               clusterOptions = markerClusterOptions(disableClusteringAtZoom = 13), 
               label = ~ lapply(label, htmltools::HTML))
  map
}

