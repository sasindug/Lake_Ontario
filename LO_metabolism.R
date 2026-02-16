##### Glider Metabolism Script #####
##### Sasindu Gunawardana 2026 #####

# Required packages
library(rLakeAnalyzer)
library(LakeMetabolizer)
library(dplyr)
library(limSolve)
library(stringr)
library(readr)
library(geosphere)
library(lutz)
library(chron)
library(imputeTS)
library(streamMetabolizer)
library(lubridate)
library(bayestestR)

# Set working directory to a folder on your computer that contains the raw data 
# from the glider DAC and wind data from the NOAA national buoy data center. 
# Save optional data here too. 

##### EDIT CODE BELOW #####
setwd("")
##### EDIT CODE ABOVE #####

# Read in the raw data from DAC. Need to use the readr package to increase
# speed in working with big data. Most computers should have enough RAM to 
# handle the glider data, but close other programs to conserve working memory. 

# Make sure your raw data is downloaded as .csvp from the DAC. This will ensure
# that all column information is in the first row. 
##### EDIT CODE BELOW #####
raw.data<-read_csv("nokomis-20141030T1510-delayed_7615_516f_51a7.csv")
##### EDIT CODE ABOVE #####

# Each glider dataset is formatted differently. We first need to edit the column
# names to the same format as this script. View column names and change the ones
# that are needed. 
colnames(raw.data)

##### EDIT CODE BELOW #####

raw.data<-raw.data%>%
  rename(DateTime = `precise_time (UTC)`, depth = `depth (m)`, 
         temperature = `temperature (Celsius)`, 
         chlor = `CPHL (ug/l)`,
         precise_lat = `precise_lat (degree_north)`, 
         precise_lon = `precise_lon (degree_east)`,
         sci_oxy4_oxygen = `sci_oxy4_oxygen (uM)`,
         conduc = `conductivity (S m-1)`,
         cdom = `CDOM (ppb)`,
         salinity = `salinity (1)`,
         density = `density (kg m-3)`)

##### EDIT CODE ABOVE #####

# Sometimes there is depth data that is negative
raw.data<-raw.data%>%mutate(depth = case_when(depth < 0 ~ 0, depth > 0 ~ depth))

# This script requires a date time format of yyyy-mm-dd hh:mm:ss. This is the 
# default format for the glider data and the default format in R. We just need
# to specify that these data are formatted this way. 
raw.data$DateTime<-as.POSIXct(raw.data$DateTime, tz = "GMT")

# The metabolism script has the 'day' starting at 9AM. Although a full 24 hours
# is not required, it is important that there is a light and dark period in the 
# time series. Cut the data to start at 9AM. You can also use this to filter the 
# data to only the time period you are interested in

# Use this to view the time extent of the glider mission. 
min(raw.data$DateTime)
max(raw.data$DateTime)

write.csv(data.frame(Start = min(raw.data$DateTime), End = max(raw.data$DateTime)), "Time_Range.txt")

##### EDIT CODE BELOW #####
raw.data<-raw.data%>%
  filter(DateTime >= as.POSIXct("2014-10-31 09:00:00", tz = "GMT"))

### Ascent and descent speed calculations
speeddata <- raw.data %>% select(profile_id, depth, DateTime)
speed_data <- speeddata %>%
  group_by(profile_id) %>%
  arrange(DateTime) %>%
  mutate(
    depth_diff = lead(depth) - depth,
    time_diff = as.numeric(difftime(lead(DateTime), DateTime, units = "secs")),
    speed = depth_diff / time_diff,
    direction = ifelse(depth_diff > 0, "ascending", "descending")
  ) %>%
  ungroup()

mean_speeds_mission <- speed_data %>%
  filter(!is.na(speed)) %>%
  group_by(direction) %>%
  summarize(mean_speed = mean(abs(speed), na.rm = TRUE))

mean_speeds_mission
##### EDIT CODE ABOVE #####


# Format the raw wind data from NOAA for the script. Note: wind data is not 
# formatted the same way for each buoy. Must check formats on a dataset basis. 
# Some glider missions span large distances and cover areas represented by 
# different buoys. This script takes that into account automatically using 
# latitudes and longitudes. You must read in at least one wind dataset and 
# assign it to 'W' and 'Wind'. All extra datasets (up to 7) can be formatted 
# as 'W1' and 'Wind1' and so on without additional changes to the script.
# Wind speed needs to be in m/s.
# Wind data set has numeric and character data, so we have to comvert them into
# numeric columns.


##### EDIT CODE BELOW #####
# Read in the data for the first buoy
Wind<-read_csv("Wind.csv")

# Specify the spatial location of the first buoy
W<-data.frame(precise_lat = 46.814, precise_lon = -91.829)
##### EDIT CODE ABOVE #####

# Format wind data date time column
Wind<-Wind[-1,]%>%rename(YY = '#YY', Speed = WSPD)
Wind <- as.data.frame(sapply(Wind, as.numeric ))
Wind$DateTime<-as.POSIXct(paste(Wind$YY,"-",Wind$MM,"-",Wind$DD," ",Wind$hh,
                                ":",Wind$mm,":00", sep = ""), tz = "GMT")

# Different buoys have different time steps, likewise there can be missing data
# generate a time series by the time step of the buoy (e.g., 10 min)
##### EDIT CODE BELOW #####
wdatetime<-data.frame(DateTime = seq.POSIXt(from = min(Wind$DateTime, na.rm = T), 
                                            to = max(Wind$DateTime, na.rm = T),
                                            by='10 min'))
##### EDIT CODE ABOVE #####

# Join the wind data with the time series and filter it to the start of your
# glider dataset
Wind<-left_join(wdatetime, Wind, by = "DateTime")
##### EDIT CODE BELOW #####
Wind<-Wind%>%filter(DateTime >= as.POSIXct("2014-10-31 09:00:00", tz = "GMT"))
##### EDIT CODE ABOVE #####

# We need to measure time in the units of day from the start of our time series.
# We will do this for the glider data below. Add a column with time as 'days'.
# Make sure the time step ('by = ') is the same as the dataset. For example, 
# 600 seconds in every 10 min time step divided by 86400 seconds in a day.
##### EDIT CODE BELOW #####
Wind$days<-seq(from = 0.375, by = 600/86400, length.out = nrow(Wind))
##### EDIT CODE ABOVE #####


# Add multiple buoys as necessary. Remove the code below if not in use. 

###### EDIT CODE ABOVE #####


# This next bit of code formats the glider data for the analysis. Glider 
# datasets contain missing data whether it is just the odd second or a large
# chunk. We need to account for this by creating a time series from the min to
# max of the glider dataset by seconds. 

datetime<-data.frame(DateTime = seq.POSIXt(from = min(raw.data$DateTime, na.rm = T),
                                           to = max(raw.data$DateTime, na.rm = T), by='sec'))

# The glider data also contains repeated measurements within one second. This 
# next code averages these repeated measurments.
data<-raw.data%>%group_by(DateTime)%>%
  summarize(profile_id = min(profile_id, na.rm = T), 
            precise_lat = mean(precise_lat, na.rm = T),
            precise_lon = mean(precise_lon, na.rm = T),
            sci_oxy4_oxygen = mean(sci_oxy4_oxygen, na.rm = T),
            depth = mean(depth, na.rm = T), 
            temperature = mean(temperature, na.rm = T),
            chl = mean(chlor, na.rm = T),
            conduc = mean(conduc, na.rm = T),
            cdom = mean(cdom, na.rm =T),
            salinity = mean(salinity, na.rm = T),
            density = mean(density, na.rm = T)) %>%as.data.frame()

# Profile ID is important and sometimes there is INF data in here. We need to 
# ensure that this column is complete. Change any INF data to NA. 
data$profile_id<-ifelse(is.infinite(data$profile_id), NA, data$profile_id)

# Final join the data by the complete time series and assign any missing data to 
# NA. 
seq.data<-left_join(datetime, data, by = "DateTime")

##### EDIT/CHECK CODE BELOW #####
# The units of DO differ with the gliders. Need to select the code that formats
# the appropriate units. 

# If DO is given as a %, you need to convert that percent to a concentration
# to do this you need to know the atmospheric pressure in your region in atm. We 
# have continuous pressure in hPa from the wind dataset which we can average and
# multiply by a conversion factor

Pressure<-1024*0.0009869233

# Using pressure and temperature the DO concentration at 100% saturation is 
# determined (DOsat) and used to calculate the concentration of DO from the 
# measured %DO data (sci_oxy4_oxygen)

# seq.data$DOsat<-(exp(7.7117 - 1.31403 * log(seq.data$temperature + 45.93)))*Pressure*(((1-(exp((11.8571-3840.70/(seq.data$temperature+273.15)-216961/(seq.data$temperature+273.15)^2))/0.994621))*(1-(0.000975-0.00001426*seq.data$temperature+0.00000006436*seq.data$temperature^2)*0.994621))/((1-exp((11.8571-3840.70/(seq.data$temperature+273.15)-216961/(seq.data$temperature+273.15)^2)))*(1-0.000975-0.00001426*seq.data$temperature+0.00000006436*seq.data$temperature^2)))   
# seq.data$sci_oxy4_oxygen<-seq.data$DOsat*(seq.data$sci_oxy4_saturation/100)

# Sometimes the concentration is in umol. This is a much simpler conversion. If
# it is measured in umol, the dataset often has %DO reported as well. If not, 
# it will need to be calculated using some edits to the code above. 
seq.data$sci_oxy4_oxygen<-seq.data$sci_oxy4_oxygen*32/1000
# temp doesnt line up with concentration
seq.data$temperature<-na_interpolation(seq.data$temperature)
seq.data$DOsat<-(exp(7.7117 - 1.31403 * log(seq.data$temperature + 45.93)))*Pressure*(((1-(exp((11.8571-3840.70/(seq.data$temperature+273.15)-216961/(seq.data$temperature+273.15)^2))/0.994621))*(1-(0.000975-0.00001426*seq.data$temperature+0.00000006436*seq.data$temperature^2)*0.994621))/((1-exp((11.8571-3840.70/(seq.data$temperature+273.15)-216961/(seq.data$temperature+273.15)^2)))*(1-0.000975-0.00001426*seq.data$temperature+0.00000006436*seq.data$temperature^2)))   
seq.data$sci_oxy4_saturation<-(seq.data$sci_oxy4_oxygen/seq.data$DOsat)*100


##### EDIT/CHECK CODE ABOVE #####

# This section of code calculates the average lat lon per day
# used to map metabolism
latlong <- seq.data%>%mutate(Year = str_sub(DateTime,1,4),
                             Month = str_sub(DateTime,6,7),
                             Day = str_sub(DateTime,9,10), 
                             Hour = str_sub(DateTime,12,13))

latlong<-latlong%>%group_by(Year, Month, Day, Hour)%>%
  summarize(lat = mean(precise_lat, na.rm = T),long = mean(precise_lon, na.rm = T))
latlong$Date<-paste(latlong$Year,latlong$Month,latlong$Day, sep = "-")
latlong$DateTime<-paste(latlong$Date," ", latlong$Hour, ":00", sep = "")

write.csv(latlong, "mission_latlong.csv")

# Like the wind data, we need to measure time in the units of 'days'. Add column
# to the glider dataset that measures time in the units of days by second
seq.data$Days<-seq(from = 0.375, by = 1/86400, length.out = nrow(datetime))

# Latitude, longitude, depth and profile id need to be complete, interpolate any 
# missing values. 
seq.data$profile_id<-na_locf(seq.data$profile_id)
seq.data$precise_lat<-na_locf(seq.data$precise_lat)
seq.data$precise_lon<-na_locf(seq.data$precise_lon)
seq.data$depth.in<-na_interpolation(seq.data$depth)

# Remove unused dataframes from the R memory. 
rm(data, datetime, raw.data)

# Next, we need to add date as a unique identifier with the code below. Note:
# day is defined as 9:00 AM to 9:00 AM, so it gets a little complicated. 

seq.data$Date<-str_sub(seq.data$DateTime,1,10)
seq.data$row<-row.names(seq.data)
day.int<-seq(from = 0.375, to = length(unique(seq.data$Date))+0.375, by = 1)
day.int.mat<-as.data.frame(cbind(day.int,day.int))
colnames(day.int.mat)<-c("day.int","day.int.lag")
day.int.mat$day.int.lag<-lead(day.int.mat$day.int.lag)
day.int.mat<-na.omit(day.int.mat)
day.table<-data.frame(row = NA, TID=NA)

for(i in 1:nrow(day.int.mat)){
  
  result<-seq.data%>%dplyr::select(Days, row)%>%filter(Days >= day.int.mat[i,"day.int"]&Days <= day.int.mat[i,"day.int.lag"])%>%
    mutate(TID = i)%>%dplyr::select(row, TID)
  
  day.table[result$row,]<-result
  day.table<-na.omit(day.table)
  
}

seq.data<-left_join(seq.data, day.table, by = "row")

# Remove unused dataframes from the R memory. 
rm(day.int.mat, day.table, result)

# Due to the movement of the glider, there is a time lag in DO and depth, we 
# need to determine that time lag and correct for it by comparing the ascending 
# and descending data. The code below lags the series, compares the two datasets
# and determines the lag that has the lowest rmse between profiles. It is
# usually around 10 seconds

# Separate profiles based on ascending vs descending
seq.data<-seq.data%>%mutate(Even = ifelse(profile_id %% 2 == 0, "Y", "N"))
dec.data<-seq.data%>%filter(Even == "N")%>%dplyr::select(-Even)
asc.data<-seq.data%>%filter(Even == "Y")%>%dplyr::select(-Even)

final.table<-data.frame(Time = 0:40)
for(x in unique(dec.data$profile_id)){
  
  dec.filter<-seq.data%>%filter(profile_id == x)
  asc.filter<-seq.data%>%filter(profile_id == (x+1))
  
  data<-rbind(dec.filter, asc.filter)
  head(data)
  
  ID<-paste("ID",x, sep = "")
  output.table<-data.frame(ID = 0:40, rme = NA)
  colnames(output.table)<-c("Time", ID)
  
  for(b in seq(0, 40, by = 1)){
    
    data$DOL<-lead(data$sci_oxy4_oxygen, b)
    
    v<-data[data$profile_id == x+1, "depth.in"]
    q<-data[data$profile_id == x+1, "DOL"]
    
    w<-data[data$profile_id == x, "depth.in"]
    e<-data[data$profile_id == x, "DOL"]
    
    dq<-approxfun(v, q)
    aq<-approxfun(w, e)
    
    t<-seq(0,floor(max(seq.data$depth.in)), by = 1)
    
    rme<-dq(t)-aq(t)
    rme<-rme^2
    rme<-sqrt((sum(rme, na.rm = T))/(floor(max(seq.data$depth.in))+1-sum(is.na(rme))))
    rme
    
    result<-c(b, rme)
    row <- which(output.table$Time==b) 
    output.table[row,]<-result
    
    
    #plot(t, dq(t), xlab = "Depth (m)", ylab = "DO (mg/L)")
    #lines(t, aq(t))
  }
  
  output.table<-data.frame(output.table[,-1])
  colnames(output.table)<-c(ID)
  
  final.table<-cbind(final.table, output.table)
  
}

offset<-rowMeans(final.table, na.rm = T)
mean(offset)
offset<-as.data.frame(cbind(Time = 0:40, offset))
offset.adj<-offset[offset$offset == min(offset$offset, na.rm = T), "Time"]
#plot(offset$Time, offset$offset, type = "o", pch = 19, col = "red")

# The above code will automatically determine the best offset time. The code
# below will apply that time to the oxygen data. 
input.data<-seq.data%>%select(-Even, -depth.in)
input.data$sci_oxy4_oxygen<-lead(input.data$sci_oxy4_oxygen, offset.adj)
input.data$sci_oxy4_saturation<-lead(input.data$sci_oxy4_saturation, offset.adj)

write.csv(input.data, "C:/Users/sasin/Sync/Sasindu/Trent Metabolism/Lake Superior/Outputs/Input tables/input.data.csv")

# Calculate the mean depth for each day
daily_mean_depth <- input.data %>%
  mutate(date = as.Date(DateTime)) %>%  # Convert DateTime to date
  group_by(date) %>%                    # Group by date
  summarize(mean_depth = mean(depth, na.rm = TRUE))

write.csv(daily_mean_depth, "C:/Users/sasin/Sync/Sasindu/Trent Metabolism/Lake Superior/daily_max_depth/mean_depth/daily_mean_depth_2014_10_30.csv")


# Remove unused dataframes/vectors from R. 
rm(asc.data, asc.filter, data, dec.data, dec.filter, final.table, offset, output.table,
   seq.data, b, day.int, e, i, ID, offset.adj, q, result, rme, row, t, v, w, x, aq, dq,
   wdatetime)

### Here we make a ggplot to see temperature profile of each day with depth
filtered_data <- input.data %>% filter(depth <= 60) # to make graph clear

# Plot the graph
ggplot(filtered_data, aes(x = temperature, y = depth)) +
  geom_point(color = "blue", size = 1) +
  scale_y_reverse() +
  labs(x = "Temperature",y = "Depth (m)") +
  facet_wrap(~ Date) +
  theme_minimal(base_family = "Times New Roman") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12)) 



# Define metadepth function. It used the meta.depths() function from 
# rlakeAnalyzer but it is edited to return only the epilimnion depth. Also, if 
# it returns NA and is fully mixed, it uses a depth of 15m. If it returns a very
# shallow value, it uses a depth of 3m. 
m.d<-function(x,y){d<-meta.depths(x,y)
if(is.na(d[1])){15}else{if(d[1]< 3){3}else{d[1]}}}

# Set up output table for the metabolism part of the script. 
output.table<-data.frame(TID = unique(input.data$TID), O2.mean = NA, GPP.mean = NA, ER.mean = NA,
                         O2.sd = NA, GPP.sd = NA, ER.sd = NA, R2 = NA, MetaDepth = NA, Chl = NA, 
                         conduc = NA, O2.Sat = NA, temperature = NA, salinity = NA, K = NA)

# This loop will automatically calculate daily metabolism from the glider data
for(i in unique(input.data$TID)){
  
  # Filter data for each day
  data<-input.data%>%filter(TID == i)
  
  # This if statement checks if there is missing glider data. If more than 80% 
  # is missing per day, metabolism is not calculated
  if(sum(is.na(data$depth))/nrow(data) > 0.8)
  {
    result<-c(i, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA) } else{
      # Interpolate depth
      data$depth<-na_interpolation(data$depth)
      # Round depth to integer
      data<-data%>%mutate(Depth.I = floor(depth))
      # Create a new unique ID column
      data$FID <- factor(data$profile_id):factor(data$Depth.I)
      # Select the columns to be analyzed
      vars<-c("Days","Depth.I","profile_id", "precise_lat","precise_lon",
              "sci_oxy4_oxygen", "sci_oxy4_saturation","temperature", "chl", "conduc", "salinity")
      # Summarize data at a 1 meter interval for each profile
      data<-data%>%group_by(FID)%>%
        summarize_at(vars, mean, na.rm = T)%>%arrange(Days)%>%dplyr::select(-FID)%>%as.data.frame()%>%na.omit()
      
      # Calculate the epilimnion depth for each profile and determine the 25% 
      # quartile epilimnion depth from all of the profiles 
      metadepth.list<-data%>%group_by(profile_id)%>%summarize(depths = m.d(temperature, Depth.I))%>%as.data.frame()
      metadepth<-as.numeric(quantile(metadepth.list$depths, prob = 0.25, na.rm = T))
      # Filter the dataset to remove all of the data below the epilimnion depth
      # and average each variable across the epilimnion
      metab.data<-data%>%filter(Depth.I <= metadepth)%>%group_by(profile_id)%>%
        summarize(Day = mean(Days), O2mean = mean(sci_oxy4_saturation, na.rm = T), chl = mean(chl, na.rm = T), 
                  conduc = mean(conduc, na.rm = T), salinity = mean(salinity, na.rm =T),
                  O2 = auc(Depth.I, sci_oxy4_oxygen)/(max(Depth.I)-min(Depth.I)),
                  O2satper = auc(Depth.I, sci_oxy4_saturation)/(max(Depth.I)-min(Depth.I)), temperature = mean(temperature,na.rm =T))%>%
        arrange(Day)%>%mutate(O2sat = O2/(O2satper/100))%>%dplyr::select(-O2satper)%>%as.data.frame()
      
      
      # To determine the change in oxygen due to mixing we need to know the 
      # oxygen right below (+2 meters) the epilimnion for each profile
      mix<-data%>%filter(Depth.I >= metadepth&Depth.I <= metadepth+2)%>%group_by(profile_id)%>%
        summarize(D.O2 = auc(Depth.I, sci_oxy4_oxygen)/(max(Depth.I)-min(Depth.I)))%>%as.data.frame()
      metab.data<-left_join(metab.data, mix, by = "profile_id")
      # Next we need to model PAR based on the spatial location of the glider
      tz<-tz_lookup_coords(lat = mean(data$precise_lat, na.rm = T), lon = mean(data$precise_lon, na.rm = T), method = "accurate")
      hm <- merge(0:23, seq(0, 45, by = 15))
      datetime <- merge(unique(input.data$Date), chron(time = paste(hm$x, ':', hm$y, ':', 0)))
      colnames(datetime) <- c('date', 'time')
      datetime$dt <- as.POSIXct(paste(datetime$date, datetime$time))
      datetime <- datetime[order(datetime$dt), ]
      datetime$dt <- as.POSIXct(datetime$dt, tz=tz)
      row.names(datetime) <- NULL
      
      datetime$Days <- seq(from = 0, to = (nrow(datetime)-1)*(900/86400), by = (900/86400))
      datetime$solar.time<-calc_solar_time(datetime$dt, longitude = mean(data$precise_lon, na.rm = T))
      
      datetime$PAR<-calc_light(datetime$solar.time, latitude = mean(data$precise_lat, na.rm = T), longitude = mean(data$precise_lon, na.rm = T))
      # The PAR time series is then converted into a function in order to 
      # estimate the light at specific points in time (units of 'days')
      fpar<-approxfun(datetime$Days, datetime$PAR, method = "linear", rule = 2)
      
      # Here, we determine the closest buoy (up to 7) and create a function to 
      # estimate wind at specific time points
      Wind.Sel<-data.frame(W = distm(c(mean(data$precise_lon), mean(data$precise_lat)), c(W$precise_lon, W$precise_lat), fun = distHaversine),
                           W1 = tryCatch(distm(c(mean(data$precise_lon), mean(data$precise_lat)), c(W1$precise_lon, W1$precise_lat), fun = distHaversine), error=function(err) NA),
                           W2 = tryCatch(distm(c(mean(data$precise_lon), mean(data$precise_lat)), c(W2$precise_lon, W2$precise_lat), fun = distHaversine), error=function(err) NA),
                           W3 = tryCatch(distm(c(mean(data$precise_lon), mean(data$precise_lat)), c(W3$precise_lon, W3$precise_lat), fun = distHaversine), error=function(err) NA),
                           W4 = tryCatch(distm(c(mean(data$precise_lon), mean(data$precise_lat)), c(W4$precise_lon, W4$precise_lat), fun = distHaversine), error=function(err) NA),
                           W5 = tryCatch(distm(c(mean(data$precise_lon), mean(data$precise_lat)), c(W5$precise_lon, W5$precise_lat), fun = distHaversine), error=function(err) NA),
                           W6 = tryCatch(distm(c(mean(data$precise_lon), mean(data$precise_lat)), c(W6$precise_lon, W6$precise_lat), fun = distHaversine), error=function(err) NA),
                           W7 = tryCatch(distm(c(mean(data$precise_lon), mean(data$precise_lat)), c(W7$precise_lon, W7$precise_lat), fun = distHaversine), error=function(err) NA),
                           W8 = tryCatch(distm(c(mean(data$precise_lon), mean(data$precise_lat)), c(W8$precise_lon, W8$precise_lat), fun = distHaversine), error=function(err) NA),
                           W9 = tryCatch(distm(c(mean(data$precise_lon), mean(data$precise_lat)), c(W9$precise_lon, W9$precise_lat), fun = distHaversine), error=function(err) NA))
      
      Wind.Col<-colnames(Wind.Sel)
      WS<-Wind.Col[c(which(Wind.Sel == min(Wind.Sel, na.rm = T)))]
      
      if(WS == "W") {
        fwind <- approxfun(Wind$days, Wind$Speed, method = "linear", rule = 2)} else {
          if(WS == "W1") {fwind <- approxfun(Wind1$days, Wind1$Speed, method = "linear", rule = 2)} else {
            if(WS == "W2") {fwind <- approxfun(Wind2$days, Wind2$Speed, method = "linear", rule = 2)} else {
              if(WS == "W3") {fwind <- approxfun(Wind3$days, Wind3$Speed, method = "linear", rule = 2)} else {
                if(WS == "W4") {fwind <- approxfun(Wind4$days, Wind4$Speed, method = "linear", rule = 2)} else {
                  if(WS == "W5") {fwind <- approxfun(Wind5$days, Wind5$Speed, method = "linear", rule = 2)} else {
                    if(WS == "W6") {fwind <- approxfun(Wind6$days, Wind6$Speed, method = "linear", rule = 2)} else {
                      if(WS == "W7") {fwind <- approxfun(Wind7$days, Wind7$Speed, method = "linear", rule = 2)} else {
                        if(WS == "W8") {fwind <- approxfun(Wind8$days, Wind8$Speed, method = "linear", rule = 2)} else {
                          fwind <- approxfun(Wind7$days, Wind7$Speed, method = "linear", rule = 2)}}}}}}}}} 
      
      # Use the functions to estimate PAR and wind by day time
      metab.data$PAR<-fpar(metab.data$Day)
      metab.data$Wind<-fwind(metab.data$Day)
      
      # Calculate PAR parameters for the metabolism equation
      day.x<-seq(from = floor(min(metab.data$Day))+0.375, to = floor(min(metab.data$Day))+1.375, by = 0.005)
      E.par<-fpar(day.x)
      E<-auc(day.x,E.par)
      
      
      # Estimate gas exchange using the wind data. This is just a simple gas
      # exchange equation. Also estimate the amount of mixing.  
      metab.data<-metab.data%>%mutate(K = k.vachon.base(metab.data$Wind, lake.area = 82103000000, params=c(2.51,1.48,0.39)))%>%
        mutate(Fatm = K*(O2sat - O2), Mix.O2 = (D.O2-O2)*0.000055)
      # Extract the average chlorophyll, conductivity, cdom,  %DO, temperature and salinity per day
      chl<-mean(metab.data$chl, na.rm = T)
      conduc<-mean(metab.data$conduc, na.rm = T)
      cdom<-mean(metab.data$cdom, na.rm = T)
      O2mean<-mean(metab.data$O2mean, na.rm = T)
      temperature<-mean(metab.data$temperature, na.rm = T)
      salinity <-mean(metab.data$salinity, na.rm = T)
      K <- mean(metab.data$K, na.rm = T)
      # This if statement checks if there is missing DO data. If there is too 
      # much missing data, the rest of the metabolism script is not run. 
      if(sum(is.na(metab.data$O2))/nrow(metab.data) > 0.5 || 
         sum(is.na(metab.data$Mix.O2))/nrow(metab.data) > 0.5)
      {
        result<-c(i, NA, NA, NA, NA, NA, NA, NA, metadepth, chl, conduc, O2mean, temperature, salinity, K) } else{
          
          fFatm<-approxfun(metab.data$Day, metab.data$Fatm)
          fMix<-approxfun(metab.data$Day, metab.data$Mix.O2)
          fPAR<-approxfun(metab.data$Day, metab.data$PAR)
          
          metab.data<-na.omit(metab.data)
          t0<-metab.data[1,c("Day","PAR","Fatm","Mix.O2")]
          metab.data<-metab.data%>%rowwise()%>%mutate(
            aucFatm = auc(seq(t0[,"Day"], Day, length.out = 1000), 
                          fFatm(seq(t0[,"Day"], Day, length.out = 1000)))/floor(metadepth),
            aucMix = auc(seq(t0[,"Day"], Day, length.out = 1000), 
                         fMix(seq(t0[,"Day"], Day, length.out = 1000)))/floor(metadepth),
            aucPAR = auc(seq(t0[,"Day"], Day, length.out = 1000), 
                         fPAR(seq(t0[,"Day"], Day, length.out = 1000)))/E)%>%
            mutate(Obio = O2 - aucFatm - aucMix, dt = Day - t0[,"Day"])%>%as.data.frame()%>%na.omit()
          
          
          # Matrix algebra to find the soloution to the barone equation
          A<- matrix(c(rep(1,nrow(metab.data)),metab.data$aucPAR,metab.data$dt),
                     nrow = nrow(metab.data), ncol = 3, byrow = F) 
          
          b<-matrix(c(metab.data$Obio),
                    nrow = nrow(metab.data), ncol = 1, byrow = F)
          # This part is really important as it constrains GPP to be positive
          # and ER to be negative
          G<-diag(3)
          G[3,3]<-c(-1)
          
          H<-c(0,0,0)
          # Solve the matrix algebra equation
          xs<-as.data.frame(xsample(A = A, B = b, G = G, H = H)$X)
          colnames(xs)<-c("O2","GPP","ER")
          xs<-xs%>%filter(GPP > 0&ER < 0)%>%summarise(O2.mean = mean(O2), GPP.mean = mean(GPP), ER.mean = mean(ER),
                                                      O2.sd = sd(O2), GPP.sd = sd(GPP), ER.sd = sd(ER))
          # Plug the coefficients back into the equation to see how well they
          # predict the observed data
          metab.data$pred<-xs$O2.mean + xs$GPP.mean*metab.data$aucPAR + xs$ER.mean*metab.data$dt
          R2<-cor(metab.data$Obio, metab.data$pred)^2
          # Export the results
          result<-c(i, xs$O2.mean, xs$GPP.mean, xs$ER.mean,
                    xs$O2.sd, xs$GPP.sd, xs$ER.sd, R2, metadepth, chl, conduc, O2mean,temperature, salinity, K)
        }
    }
  row <- which(output.table$TID==i) 
  output.table[row,]<-result
  
}

output.table$Date<-unique(input.data$Date)

write.csv(output.table, "C:/Users/sasin/Sync/Sasindu/Trent Metabolism/Meta_metabolism/Comtarison_epi_meta/epi_2014_10_30.csv")


write.csv(output.table, "C:/Users/sasin/Sync/Sasindu/Trent Metabolism/Lake Superior/Vachan_P_Output/output_2014_10_30.csv")

write.csv(input.data, "C:/Users/sasin/Sync/Sasindu/Trent Metabolism/Lake Superior/Outputs/input_2014_10_30.csv")


### DO Sat data
# Surface DO_Sat data #####
surface.data <- input.data %>%
  filter(depth <= 3)

# Convert datetime to date and aggregate by day
daily_do_sat <- surface.data %>%
  mutate(date = as.Date(DateTime)) %>%
  group_by(date) %>%
  summarize(daily_DO_Sat = mean(sci_oxy4_saturation, na.rm = TRUE))

write.csv(daily_do_sat, "C:/Users/sasin/Sync/Sasindu/Trent Metabolism/Lake Superior/DO_Sat/surface_do_sat_2014_10_30.csv")

