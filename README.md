# monarch roost analysis

This project aims to assess how monarch Fall migration through the Midwest flyway has changed over time. The conservation status of monarch butterflies is complicated by the dynamic, long-distance migration to and from overwintering colonies, as well as regional variation in climate, alterations in landscape composition and configuration, milkweed reductions due (in part) to pesticide use, etc. What is clear is that the number of overwintering monarchs in Mexico has been declining steadily over the past few decades (assuming these estimates approximate abundance well). It is possible that overwintering colonies are declining due to localized declines in the summer populations, or increasing mortality associated with migration. The current assumption is that the majority of monarchs survive the migration, and thus declines in the overwintering colonies reflect global declines in summer populations. 

This study hopes to make a meaningful contribution to the conversation about monarch migration, while also elucidating some basic biology of monarch roosting behavior. We leverage a vast database of citizen science observations of monarch roosts collected by Journey North participants during the Fall migration. Observations include roost size estimates, coordinates, and time stamps. The way Journey North operates is that participating scientists report early  migration observations and are further notified of the coming wave of monarchs. So we do not have continuously monitored roosting sites, but rather observations that are distributed over space and time to capture the peak of the migrating wave of monarchs across the flyway. We opted for a spatially varying coefficient model that allows the effect of time to vary spatially. This way, we can smooth roost size across the flyway, and examine how this spatial plane evolves over the duration of the study. Information on roost size and thus trends borrow from adjacent sites, and, to my understanding, uncensused areas / locations are essentially interpolated, weighted by neighboring observations. 

Getting back to the biology, there are few patterns we might expect that would lend support for different hypotheses. If summer monarch populations are declining, then fewer monarchs participate in migration and thus fewer will arrive at their overwintering sites in Mexico, lending support for the "local population declines hypothesis". If summer populations are stable, and if monarchs are experiencing increasing mortality during migration each year (due to increasing OE disease, changes in weather/climate, etc.) will generate latitudinal clines in roost size trends and should explain monarch declines at the overwintering sites, lending support to the "migrational mortality hypothesis". We will first examine the patterns and then seek to identify the processes underlying such patterns.

The Objectives of the analysis proceed as follows. The first objective is purely description--we will use a space-time model to get estimates of spatial variation in roost size trends. We will compare this with a fixed effect model (year as fixed effect vs SVC) and an intercept model to see whether there are any trends in roost size and, if so, which model best captures changes in roost size. For Objective 2, we seek to explain spatial variation in trends using the hypothesis weather like wind speed, wind direction (tail vs head winds),temperature, precipitation, etc., and landscape features like NDVI (a proxy for nectar availability) both influence roost size and might, due to changes in these variables over time, explain spatio-temporal roost size trends. We extract information on weather and landscape conditions at a particular site and compare via model selection.


1. [curate and inspect.R](curate and inspect.R) takes the original [Journey North data file](monarch roost_original.xlsx) and performs the following:
  + takes the maximum observation per location per year--including multiple observations from a single observer for a given year may equate to taking the mean of the observations for that year, which would reduce the effective estimate from that locale. I don't think we want this because we are trying to explain variation in the size of the wave of monarchs migrating through. 
  + tidies up column names / formatting / etc.
  + creates some GIS utilities/base naps and saves them as [GIS utilities.RData](GIS utilities.RData)
  + inspects the amount of data available and the spatial distribution of observations, creates supplementary figures
  + writes data file to [monarch roost_curated.csv](monarch roost_curated.csv)

2. Weather variables: [visual crossing weather extraction.R](visual crossing weather extraction.R) scrapes and tidies a dataframe of weather observations up to 7 days prior to each observation, saving observation IDs in [weather.rda](weather.rda), raw downloaded data in [vcweather.rda](vcweather.rda), and tidied csv in [vcweather.csv](vcweather.csv). Weather data was further curated to produce 7day averages in [weather_curated.R](weather_curated.R), and this was written to [weather_curated.csv](weather_curated.csv).

3. NDVI was extracted in [ndvi extract_same day.R](ndvi extract_same day.R) to produce [ndvi_same day.csv](ndvi_same day.csv)

4. [fixed model list.R](fixed model list.R) creates a list of all models to be compared via model selection

5. [model selection_final.R](model selection_final.R) runs all models, compares, conducts model averaging, etc.  The following is created/produced:
  + modeling mesh
  + SPDE object
  + model output
  + WAIC table
  + plot results of full model in [plots](plots)
  + etc.
  