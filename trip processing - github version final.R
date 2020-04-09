######################
######################
#
# Greenspace visit GPS data processing
# Author: Meghann Mears meghann.mears@googlemail.com
# Github repository: 
#
######################
######################

library(tidyverse) 
library(sf)
library(data.table) # for rbindlist function
library(lubridate)

######################
# Data loading
 
# Import and pre-process GPS data file
# This script assumes that the GPS data is in csv format and
# that the following fields are present:
# device ID, x & y coordinates (projected system), timestamp

gps_dat <- read_csv(file.choose()) # GPS data file

# Extract/rename relevant fields - change field names as required
gps_dat <- tibble(device_id = gps_dat$DEVICE_ID_FIELD,
			easting = gps_dat$EASTING_FIELD,
			northing = gps_dat$NORTHING_FIELD,
			timestamp = as.POSIXct(gps_dat$TIMESTAMP_FIELD, format="%Y-%m-%d %H:%M:%S", 
				origin = "1970-01-01", tz="UTC")) %>%
	distinct() %>% # Remove duplicates
	arrange(device_id, timestamp) # Sort


# Load other GIS data necessary for analysis.
# These are assumed to be in a format already containing geometry
# data e.g. shapefile, GeoJSON.
# Change CRS as required (27700 = OSGB 36).

# Greenspace boundaries
gs_boundaries <- st_read(file.choose()) %>%
	st_set_crs(27700) %>% # setting CRS is needed if using shapefiles
	select("SITE_ID_FIELD") %>% # remove all columns except site ID
	rename("site_id" = "SITE_ID_FIELD") # rename site ID column

# Greenspace entrance points
gs_entrance <- st_read(file.choose()) %>%
	st_set_crs(27700)

# Buildings inside greenspaces
gs_buildings <- st_read(file.choose()) %>%
	st_set_crs(27700)

# User home coordinates - only required if distance from home to be calculated

users <- read_csv("E:/Users/megha/IWUN/WP3/shmapped/users.csv")

users <- tibble(device_id = users$DEVICE_ID_FIELD,
			easting = users$EASTING_FIELD,
			northing = users$NORTHING_FIELD) %>%
	distinct() %>%
	st_as_sf(coords=c("easting", "northing"), crs=27700, remove = FALSE)


#######################
# Divide GPS points into trips

# Data frame to hold trip separation data
new_trip <- tibble(

	# Identify changes of device id
	phones = c(FALSE, with(gps_dat, (device_id[-1] != device_id[-nrow(gps_dat)]))),

	# Identify gaps of >10 mins
	gaps = c(FALSE, with(gps_dat, difftime(timestamp[-1], timestamp[-nrow(gps_dat)], units="mins") > 10))

)

# Cumulative counter of new trips (i.e. where one or more of the criteria are met)
new_trip <- new_trip %>%
	mutate(trip_id = cumsum(gaps == TRUE | phones == TRUE))

# Assign back to gps_dat as trip ID
gps_dat$trip_id <- new_trip$trip_id


##################################
# Data cleaning

# Remove unreasonable altitudes - outside of Sheffield's max/min altitudes,
# with 10m allowance for error
gps_dat <- gps_dat[gps_dat$altitude <= (592+30) & gps_dat$altitude >= (19-18),]

# Remove very short trips - <70 seconds (or with only 1-2 points)
gps_dat_summary <- gps_dat %>%
	group_by(trip_id) %>%
	summarise(duration = difftime(max(timestamp), min(timestamp), units="secs"),
		point_count = n()) 

gps_dat <- gps_dat %>%
	left_join(gps_dat_summary, by="trip_id")

gps_dat <- gps_dat[gps_dat$duration >= 70 & gps_dat$point_count > 2,] 


# Break into ‘quality segments’ where the speed between consecutive data points
# is not faster than would be expected for non-vehicular travel. Compare each pair
# of quality segments and remove the shorter. Repeat until the trip is a single
# quality segment.

# Function to compare each data point (within trips) with the next 
# to identify where the distance travelled is greater than a permissible
# maximum speed (defaults to 7mps) + error distance (defaults to 30m).
create_segments <- function(trip, maxspeed = 7, error = 30) {

	# Initialise vars
	trip$distance <- 0
	trip$time_diff <- 0	

	# 3D Euclidean distance from previous point
	trip$distance[2:nrow(trip)] <- sqrt( (trip$northing[1:nrow(trip)-1] - trip$northing[2:nrow(trip)])^2 +
		(trip$easting[1:nrow(trip)-1] - trip$easting[2:nrow(trip)])^2 +
		(trip$altitude[1:nrow(trip)-1] - trip$altitude[2:nrow(trip)])^2 )

	# Time difference from previous point
	trip$time_diff[2:nrow(trip)] <- trip$timestamp[2:nrow(trip)] - trip$timestamp[1:nrow(trip)-1] 

	# Maximum allowable distance from previous point
	trip$max_distance <- (trip$time_diff * 7) + 30

	# If true, distance is more than allowable
	trip$split <- trip$distance > trip$max_distance 
	
	trip$segments <- cumsum(trip$split)	

	return(trip)
	}

# Function to compare length of time of two segments (identified by create_segments())
compare_segments <- function(segment1, segment2) {

	# Calculate length of time of each segment
	time_seg1 <- max(segment1$timestamp) - min(segment1$timestamp)
	time_seg2 <- max(segment2$timestamp) - min(segment2$timestamp)

	# Compare time of each segment (and hope it's not equal because I can't 
	# figure out a better way to handle that). Return the longer.
	if (time_seg1 > time_seg2) {
		return(segment1)
	} else {
		return(segment2)
	}
}

# Recursive function to fully clean trips via segmentation.
# Takes the parameters for distance_trip (but defaults set).
clean_trip <- function(trip, maxspeed = 7, error = 30) {
	
	# Get trip segments based on maximum permissible distance/speed test
	trip <- create_segments(trip, maxspeed=maxspeed, error=error)
	
	# Get rid of segments with only one GPS data point
	segment_summary <- trip %>%
		group_by(segments) %>%
		summarise(count = n())
	
	if (length(segment_summary$segments[segment_summary$count == 1]) > 0) {
		trip <- trip[-which(trip$segments %in% segment_summary$segments[segment_summary$count == 1]),]
	}

	# If the trip has more than one segment (i.e. places where distances are
	# travelled unreasonably fast...
	if (length(unique(trip$segments)) > 1) {

		# Compare the first two segments and keep the longer one.
		trip <- rbind(compare_segments(trip[trip$segments==unique(trip$segments)[1],], trip[trip$segments==unique(trip$segments)[2],]), trip[trip$segments>unique(trip$segments)[2],])

		# If there are other segments after the first two (i.e. still >1)...
		if (length(unique(trip$segments)) > 1) {

			# Run this function again. Until there is only one segment.
			trip <- clean_trip(trip, maxspeed=maxspeed, error=error)
		}
	}
	return(trip)
}

# Prepare dat to hold cleaned/segmented data
gps_dat_clean <- vector("list", length = max(gps_dat$trip_id)+1)

# Iterate through trips and segment. Takes a while - progress reported.
print(paste0("Cleaning ", length(unique(gps_dat$trip_id)), " trips"))
j <- 0 # used for reporting progress
for (i in unique(gps_dat$trip_id)) {     

	trip_dat <- gps_dat[gps_dat$trip_id == i,]
	gps_dat_clean[[i+1]] <- clean_trip(trip_dat, error=15)

	j <- j + 1
	if (j %% 1000 == 0) print(j)
}

# Un-list
gps_dat_clean <- as_tibble(rbindlist(gps_dat_clean, use.names=TRUE)) 

# Remove trips that, after cleaning, have duration < 70 or < 2 points (the latter shouldn't happen)
gps_dat_clean_summary <- gps_dat_clean %>%
	group_by(trip_id) %>%
	summarise(duration = difftime(max(timestamp), min(timestamp), units="secs"),
		point_count = n()) 

gps_dat_clean <- gps_dat_clean %>%
	left_join(gps_dat_clean_summary, by="trip_id")

gps_dat_clean <- gps_dat_clean[gps_dat_clean$duration.y > 70 & gps_dat_clean$point_count.y > 2,] 


############################
# Trip interpolation

# Function: data for each individual trip is passed in for interpolation.
# Returns a data frame of interpolated coordinates and timestamps.
# Assumes trips are from gps_dat with column names as above.
# time_gap parameter can be passed in but defaults to 10 seconds.
trip_interpolation <- function(trip_dat, time_gap = 10) {
	N <- trip_dat$northing
	E <- trip_dat$easting
	T <- trip_dat$timestamp # as.POSIXct(trip_dat$RECORDED_A, format="%Y-%m-%d %H:%M:%S")
	trip_id <- trip_dat$trip_id[1]
	device_id <- trip_dat$device_id[1]

	# If the trip lasts >time_gap secs, interpolate points at every time_gap secs(ish)
	if (difftime(max(T), min(T), unit="secs") > time_gap) {

		# Compute number of time points to be interpolated to...
		z <- ceiling(as.numeric(difftime(max(T), min(T), unit="secs"))/time_gap) + 1

		# ...and compute time points to interpolate to
		i_p<-data.frame(time = seq(from = min(T), to = max(T), length.out = z),
			trip_id = trip_id, device_id=device_id)

		# If easting needs interpolation (i.e. >1 unique easting), then interpolate
		if(length(unique(E))>1) {
			i_p$easting <- approx(T, E, xout=i_p$time)$y
		} else {
			i_p$easting<- E[1]
		}

		# If northing needs interpolation (i.e. >1 unique northing), then interpolate
		if(length(unique(N))>1) {
			i_p$northing <- approx(T, N, xout=i_p$time)$y
		} else {
			i_p$northing<- N[1]
		}

	# If trip lasts <time_gap secs... ### REALLY YOU COULD JUST DELETE THEM AT THIS POINT ###
	} else {
		# If there is more than one data point for the trip, take the first and last
		if (length(N) > 1) {
			i_p <- data.frame(northing = N[c(1,length(N))], easting = E[c(1,length(E))], 
				time = T[c(1,length(T))], trip_id=trip_id, device_id=device_id)
		# If there is only one data point for the trip, use that
		} else {
			i_p <- data.frame(northing = N[1], easting = E[1], time = T[1], trip_id=trip_id, device_id=device_id)
		}
	}

	return(i_p)
}

# Prepare list to hold output from trip_interpolation() function
interpolated_points <- vector("list", length = max(gps_dat_clean$trip_id)+1)

# Iterate through trips, using trip_interpolation() function on each
# NB this takes some time - progress printed every 1000 trips
# NB also it produces warnings that can be ignored
print(paste0("Interpolating ", length(unique(gps_dat_clean$trip_id)), " trips"))
j <- 0 # used for reporting progress
for (i in unique(gps_dat_clean$trip_id)) {     
	trip_dat <- gps_dat_clean[gps_dat_clean$trip_id == i,] 
	interpolated_points[[i+1]] <- trip_interpolation(trip_dat)
	j <- j + 1
	if (j %% 1000 == 0) print(j)
}

# Un-list
interpolated_points <- as_tibble(rbindlist(interpolated_points, use.names=TRUE)) 


###################################
# Smoothing

# NB this stage was not performed due to curving of lines introducing
# erroneous 'entries' into greenspaces. It is blocked out using:
# if (FALSE) { }
# However the following sections still use a variable named 'smoothed_points'
# so the following line puts the interpolated_points variable into it
# with geometry added - comment out this line if smoothing is performed.
smoothed_points <- interpolated_points %>%
	st_as_sf(coords=c("easting", "northing"), crs=27700, remove = FALSE)


if (FALSE) { # START BLOCKING OUT

# Function for smoothing the points (after densification, see below). 
# This is a modified and reduced version of smoothr's smooth_ksmooth function.
# (The error handling has also been removed.)
# Takes a matrix of two columns containing x and y coordinates.
# Smoothness probably needs to be quite high after densification,
# but use n = 1L.
correct_ksmooth <- function(x, smoothness = 1, n = 1L, max_distance, trip_id) {

	# Make copy of input in case of needing to pass onward
	x_keep <- x

	# Calculate bandwidth
	d_orig <- smoothr:::point_distance(x)
	bandwidth <- smoothness * (max(d_orig)*0.1)

	# Calculate pad for start and end (I think this helps the first and last
	# smoothed points be close to the raw data.)
   	pad <- list(start = rbind(x[1, ], 2 * x[1, ] - x[2, ]), 
		end = rbind(x[nrow(x), ], 2 * x[nrow(x), ] - x[nrow(x) - 1, ]))
	pad$start <- pad$start[2:1, ]

	# Calculate smoothed values
	n_pts <- nrow(x)
	x <- rbind(pad$start, x, pad$end)
	d <- c(0, cumsum(smoothr:::point_distance(x)))
	x_smooth <- stats::ksmooth(d, x[, 1], n.points = length(d)-2, 
		kernel = "normal", bandwidth = bandwidth)
	y_smooth <- stats::ksmooth(d, x[, 2], n.points = length(d)-2, 
		kernel = "normal", bandwidth = bandwidth)
	x <- cbind(x_smooth$y, y_smooth$y)[2:(length(x_smooth$y)-1), ]

	# Replace first and last values with originals
	x[1, ] <- pad$start[nrow(pad$start), ]
	x[nrow(x), ] <- pad$end[1, ]

	# If there are NAs in the smoothed data - this is likely due to the bandwidth
	# not being adequate (a problem if there are lots of closely-spaced points and
	# then one very large gap). The 0.75*max(distance) seems to stop this.
	if(sum(is.na(x)) > 0) {
		print(paste0("NAs produced, trip ", trip_id, ", returning original data"))				
		return(x_keep)
	}
	return(x)
}

# initialise list to hold results
smoothed_points <- vector("list", length = length(unique(interpolated_points$trip_id)))

# Takes a long time, print progress
print(paste0("Smoothing ", length(unique(interpolated_points$trip_id)), " trips"))
j <- 0
# for each trip....
for (i in 1:length(unique(interpolated_points$trip_id))) {
	# Extract trip data
	trip_dat <- interpolated_points[interpolated_points$trip_id==unique(interpolated_points$trip_id)[i],]

	# If there is any movement in the trip:
	if (length(unique(trip_dat$easting)) > 1 | length(unique(trip_dat$northing)) > 1 ) {

		# Smooth spatial columns
		trip_smooth <- as.data.frame(correct_ksmooth(as.matrix(cbind(trip_dat$easting,trip_dat$northing)), smoothness=2, n=1L, trip_id = trip_dat$trip_id[1]))

		# Add back other data columns
		trip_smooth$time <- trip_dat$time
		trip_smooth$trip_id <-unique(trip_dat$trip_id)
		trip_smooth$device_id <-unique(trip_dat$device_id)

	# Otherwise, don’t smooth, just use the interpolated data
	} else {
		trip_smooth <- data.frame(V1 = trip_dat$easting, V2 = trip_dat$northing, time = trip_dat$time, trip_id = trip_dat$trip_id, device_id = trip_dat$device_id)
	}

	# Put trip result in list
	smoothed_points[[i]] <- trip_smooth
	j <- j + 1
	if (j %% 1000 == 0) print(j)
}

# Un-list and add geometry
smoothed_points <- as_tibble(rbindlist(smoothed_points, use.names=TRUE)) 
names(smoothed_points)[1:2] <- c('easting','northing')
smoothed_points <- st_as_sf(smoothed_points, coords=c("easting", "northing"), crs=27700, remove = FALSE)


} # END BLOCKING OUT


##################################
# Identify periods of time spent outside of sites

# This takes each interpolated point in trips and calculated the
# proportion of the five minutes centred on that point that was spent
# within site boundaries.
# If <50%, the point is dropped. If this happens in the middle of a trip,
# the trip is split into two.
# If the trip is <5 minutes, the whole trip average is used - i.e. the
# whole trip is excluded if <50% within a greenspace.

# The same is also done for speed - where average speed is >7mps over 
# 5 minutes it is split into two.

# This is achieved using a function for calculating the moving average.
# For points at the start and end of trips where there isn't enough
# 'tail' for the true moving average to be calculated, the first and last
# values where there is enough tail are used. 
# Number of points to be averaged over can be passed in (n); default is 
# equivalent to five mins sampled every 15 secs.
ma <- function(x, n = 21) {
	res <- stats::filter(x, rep(1 / n, n), sides = 2)
	first_non_NA <- min(which(!is.na(res)))
	res[1:(first_non_NA-1)] <- res[first_non_NA]
	last_non_NA <- max(which(!is.na(res)))
	res[(last_non_NA+1):length(res)] <- res[last_non_NA]
	return(res)
}

# Intersect with greenspace boundaries (and remove duplicates due to points in multiple sites)
smoothed_points <- st_join(smoothed_points, gs_boundaries, left = TRUE) %>% 
	distinct(time, trip_id, device_id, easting, northing, .keep_all=TRUE)

# Add flag for if in site
smoothed_points$in_site <- !is.na(smoothed_points$site_id)

# Calculate speed since last point
smoothed_points$speed <- 0 
smoothed_points$speed[2:nrow(smoothed_points)] = with(smoothed_points, 
		# Euclidean distance (m)....
		(sqrt((northing[-1] - northing[-nrow(smoothed_points)])^2 + 
		(easting[-1] - easting[-nrow(smoothed_points)])^2) / 
		# over time (sec) = speed (mps)
		as.numeric(difftime(time[-1], time[-nrow(smoothed_points)], units="secs"))))

# Calculate moving average per trip.
# If trip is too short to have a moving average,
# use whole trip average instead.

# First, calculate number of points to average over for whether in greenspace:
time_gap_in_gs <- 15 # set to the time_gap parameter passed into trip_interpolation 
ma_length_in_gs <- 5*60 # set to time to average over, in seconds
n_in_gs <- ceiling(ma_length_in_gs / time_gap_in_gs) + 1

# And number to average over for speed:
time_gap_speed <- 15 # set to the time_gap parameter passed into trip_interpolation 
ma_length_speed <- 5*60 # set to time to average over, in seconds
n_speed <- ceiling(ma_length_speed / time_gap_speed) + 1

# Initialise column in tibble
smoothed_points$ma_in_gs <- NA
smoothed_points$ma_speed <- NA

# Iterate over trips
# This takes a while - progress printed
print(paste0("Calculating moving averages for ", length(unique(smoothed_points$trip_id)), " trips"))
j <- 0
for (i in unique(smoothed_points$trip_id)) {
	# Extract trip data
	trip_points <- smoothed_points[smoothed_points$trip_id == i,]

	# If the trip lasted longer than the moving average length for whether in greenspace,
	# calculate the moving average proportion of time spent in greenspace(s)
	if (dim(trip_points)[1] >= n_in_gs) {
		trip_points_ma_in_gs <- ma(trip_points$in_site, n_in_gs)
	# Else calculate proportion for entire trip 
	} else {
		trip_points_ma_in_gs <- mean(trip_points$in_site)
	}

	# If the trip lasted longer than the moving average length for speed,
	# calculate the moving average speed
	if (dim(trip_points)[1] >= n_speed) {
		trip_points_ma_speed <- ma(trip_points$speed, n_speed)
	# Else calculate proportion for entire trip 
	} else {
		trip_points_ma_speed <- mean(trip_points$speed)
	}

	# Put moving average back into tibble
	smoothed_points$ma_in_gs[smoothed_points$trip_id == i] <- trip_points_ma_in_gs
	smoothed_points$ma_speed[smoothed_points$trip_id == i] <- trip_points_ma_speed

	j <- j + 1
	if (j %% 1000 == 0) {
		print(j)
	}
} 

# Calculate new trip ID due to trips being split by moving average <50%
new_trip_2 <- tibble(
	# Check for original trip_id changes
	orig_new_trip = c(FALSE, with(smoothed_points, trip_id[-1] != trip_id[-nrow(smoothed_points)])),
	# Check for places where the moving average goes from <50% to >50%
	reenter_gs = c(FALSE, with(smoothed_points, ma_in_gs[-nrow(smoothed_points)] < 0.5 & ma_in_gs[-1] >= 0.5)),
	# Check for places where the moving average speed goes from >7mps to <7mps (~15mph)
	speed_reduces = c(FALSE, with(smoothed_points, ma_speed[-nrow(smoothed_points)] <= 7 & ma_speed[-1] > 7))
	) %>%
	mutate(new_trip_id = cumsum(orig_new_trip == TRUE | reenter_gs == TRUE | speed_reduces == TRUE))
smoothed_points$new_trip_id <- new_trip_2$new_trip_id

# remove points with moving average < 0.5
smoothed_points <- smoothed_points[smoothed_points$ma_in_gs >= 0.5 & smoothed_points$ma_speed <= 7,]


#######################################################
# Crop trip starts and ends outside of site boundaries

# Add uid for sorting
smoothed_points$uid <- 1:nrow(smoothed_points)

# Create column to store points to trim at start of trips
smoothed_points$start_trim <- FALSE
new_trip_flag <- FALSE

# Iterate through rows
for (i in 1:nrow(smoothed_points)) {
	# If first row, or a new trip id, reset flag to true
	if (i == 1) {
		new_trip_flag <- TRUE
	} else if (i > 1 & smoothed_points$new_trip_id[i] != smoothed_points$new_trip_id[i-1]) {
		new_trip_flag <- TRUE
	}

	# If trip has not yet been into a site, and this row and the next are still not 
	# in a site, set trim column to true 
	if (new_trip_flag == TRUE & is.na(smoothed_points$site_id[i]) & is.na(smoothed_points$site_id[i+1])) {
		smoothed_points$start_trim[i] <- TRUE
	# Once trip has gone into a site, set flag to false
	} else {
		new_trip_flag<- FALSE
	}
}

# Reverse order of tibble and repeat the process to trim end of trips
smoothed_points <- smoothed_points %>%
	arrange(desc(uid))

smoothed_points$end_trim <- FALSE
new_trip_flag <- FALSE

for (i in 1:nrow(smoothed_points)) {
	if (i == 1) {
		new_trip_flag <- TRUE
	} else if (i > 1 & smoothed_points$new_trip_id[i] != smoothed_points$new_trip_id[i-1]) {
		new_trip_flag <- TRUE
	}

	if (new_trip_flag == TRUE & is.na(smoothed_points$site_id[i])& is.na(smoothed_points$site_id[i+1])) {
		smoothed_points$end_trim[i] <- TRUE
	} else {
		new_trip_flag<- FALSE
	}
}

# ...and order properly again
smoothed_points <- smoothed_points %>%
	arrange(uid)

# Trim starts and ends
smoothed_points <- smoothed_points[smoothed_points$start_trim == FALSE & 
	smoothed_points$end_trim == FALSE,]


######################################
# Turn interpolated trips into line geometries

smoothed_trips <- smoothed_points %>%
	group_by(new_trip_id) %>% 
	summarise(device_id = first(device_id), 
		start_time = min(time),
		end_time = max(time),
		orig_trip_id = first(trip_id),
		do_union=FALSE) %>%  
	st_cast("LINESTRING") 

# Calculate basic trip data
smoothed_trips <- smoothed_trips %>%
	add_column(trip_len = st_length(smoothed_trips),
		duration = difftime(smoothed_trips$end_time, smoothed_trips$start_time, units="mins"),
		path_pts = mapply(smoothed_trips$geometry, FUN=length)/2) %>%
	mutate(speed_mps = as.numeric(trip_len) / (as.numeric(duration)*60))

# Once again, remove trips with duration <70 or path points <=2. #Also with length = 0.
smoothed_trips <- smoothed_trips[smoothed_trips$duration >= (70/60) & smoothed_trips$path_pts > 2,] 
#smoothed_trips <- smoothed_trips[-which(as.double(st_length(smoothed_trips$geometry)) == 0),]


######################################################################
# Calculate distance inside greenspace(s) and number of greenspace(s) visited

# NB! Percent in greenspace may be very slightly over 1. This is due to 
# the length of intersects ending up slightly different to the length of
# the original line. The difference is only in the order of magnitude of a
# couple of cm maximum.

# Dissolve boundaries for pct in gs in case of slight overlaps.
# (0 width buffer prevents topography errors)
gs_boundaries_dissolve <- st_buffer(st_combine(gs_boundaries), dist=0)

# Intersect interpolated trips and calculate intersect lengths. 
trips_gs_isect <- smoothed_trips %>%
	st_set_precision(100) %>% # precision of 1cm - this may not be high enough for all datasets
	st_intersection(gs_boundaries) 

trips_gs_dissolve_isect <- smoothed_trips %>%
	st_set_precision(100) %>% # precision of 1cm - this may not be high enough for all datasets
	st_intersection(gs_boundaries_dissolve) 

# Then summarise by new_trip_id and compute trip-level variables.
# NB make sure this has the same number of rows as smoothed_trips - 
# if any are missing, increase the precision in the previous step 
# by e.g. a factor of ten.
trips_gs_isect <- trips_gs_isect %>%
	st_drop_geometry() %>%
	group_by(new_trip_id) %>% # 
	summarise(gs_count = n_distinct(site_id))

trips_gs_dissolve_isect <- trips_gs_dissolve_isect %>%
	mutate(length = st_length(trips_gs_dissolve_isect)) %>%
	st_drop_geometry() %>%
	group_by(new_trip_id) %>% # 
	summarise(in_gs_len = sum(length))

# Join back to the trips tibble
smoothed_trips <- smoothed_trips %>%
	left_join(trips_gs_isect, by = 'new_trip_id') %>%
	left_join(trips_gs_dissolve_isect, by = 'new_trip_id') %>%
	mutate(pct_in_gs = as.numeric(in_gs_len) / as.numeric(trip_len)) 

# For sites that don’t include any movement - add in data manually.
# (They are definitely in greenspaces otherwise would have been removed
# previously - but the intersect doesn’t product a result for the trip.)
smoothed_trips[which(as.numeric(smoothed_trips$trip_len) == 0),c('in_gs_len')] <- 0
smoothed_trips[which(as.numeric(smoothed_trips$trip_len) == 0),c('gs_count','pct_in_gs')] <- 1

##########################################################
# Check for missing data at start/end of trips
# Uses distance of trip start and end to greenspace entrance point

# Get start and end points
n_row <- nrow(smoothed_points)
start_pts <- smoothed_points[c(TRUE, 
	smoothed_points$new_trip_id[2:n_row] != smoothed_points$new_trip_id[1:(n_row-1)]),] 
end_pts <- smoothed_points[c(smoothed_points$new_trip_id[1:(n_row-1)] != 
	smoothed_points$new_trip_id[2:n_row], TRUE),]

# Remove points from trips that have been removed following casting to lines
start_pts <- start_pts[start_pts$new_trip_id %in% smoothed_trips$new_trip_id,]
end_pts <- end_pts[end_pts$new_trip_id %in% smoothed_trips$new_trip_id,]

# Initialise columns to hold results
start_pts$boundary_distance <- NA
end_pts$boundary_distance <- NA

# Iterate through trips to calculate distance - progress reported
print(paste0("Calculating distance to entrances for ", length(unique(smoothed_trips$new_trip_id)), " trips"))

for (i in 1:nrow(smoothed_trips)) {

	# Find nearest entrance point
	start_nearest <- st_nearest_feature(start_pts$geometry[i], gs_boundaries$geometry)[1]
	end_nearest <- st_nearest_feature(end_pts$geometry[i], gs_boundaries$geometry)[1]

	# Calculate distance to nearest entrance point
	start_pts$boundary_distance[i] <- st_distance(start_pts$geometry[i], 
		gs_boundaries$geometry[start_nearest])[1,1]
	end_pts$boundary_distance[i] <- st_distance(end_pts$geometry[i], 
		gs_boundaries$geometry[end_nearest])[1,1]

	if (i %% 1000 == 0) {
		print(i)
	}
}

# Put distances back into tibble
smoothed_trips$start_entrance_distance <- start_pts$boundary_distance
smoothed_trips$end_entrance_distance <- end_pts$boundary_distance


##############################################
# Intersect with buildings inside greenspaces

# Intersect trips with buildings, calculate intersect length, then
# summarise by new_trip_id and calculate trip-level variable.
trips_buildings_isect <- smoothed_trips %>%
	st_set_precision(100) %>%
	st_intersection(gs_buildings) 
trips_buildings_isect <- trips_buildings_isect %>%
	mutate(isect_len = st_length(trips_buildings_isect)) %>%
	st_drop_geometry() %>%
	group_by(new_trip_id) %>%
	summarise(building_isect_len = sum(isect_len))

# Repeat, but with the points, to get approximate time spent in buildings
points_buildings_isect <- st_intersection(smoothed_points, gs_buildings)
points_buildings_isect <- points_buildings_isect %>%
	st_drop_geometry() %>%
	group_by(new_trip_id) %>%
	summarise(pts_in_buildings = n())

# Add back to tibble
smoothed_trips <- smoothed_trips %>%
	left_join(trips_buildings_isect, by="new_trip_id") %>%
	mutate(building_isect_pct = building_isect_len / trip_len) %>%
	left_join(points_buildings_isect, new="new_trip_id") %>%
	mutate(building_pts_pct = pts_in_buildings / path_pts)

# Add 0 for trips where no intersect occurs
smoothed_trips[is.na(smoothed_trips$building_isect_len), c("building_isect_len", "building_isect_pct")] <- 0
smoothed_trips[is.na(smoothed_trips$pts_in_buildings), c("pts_in_buildings", "building_pts_pct")] <- 0


#############################
# Distance from user's home to start of trip

# Initialise column
start_pts$home_distance <- NA

# Iterate through trips - progress reported
print(paste0("Calculating distance to home for ", length(unique(smoothed_trips$new_trip_id)), " trips"))

for (i in 1:nrow(start_pts)) {
	# Calculate distance to user's home postcode location
	if (start_pts$device_id[i] %in% users$device_id) {
		start_pts$home_distance[i] <- st_distance(start_pts$geometry[i], 
			users$geometry[which(users$device_id == start_pts$device_id[i])])
	}

	if (i %% 1000 == 0) {
		print(i)
	}
}

smoothed_trips$home_distance <- start_pts$home_distance


###############################
# Trip validity/certainty flags

smoothed_trips <- smoothed_trips %>%
	mutate(flag_short_distance = as.numeric(trip_len) < 25,
		flag_incomplete_data = (start_entrance_distance > 25 | end_entrance_distance > 25),
		flag_high_speed = speed_mps > 5,
		flag_buildings = (as.numeric(building_isect_pct) > 0.5 | building_pts_pct > 0.5),
		flag_outside_gs = pct_in_gs < 0.5) %>%
	mutate(flag_count = flag_short_distance + flag_incomplete_data + flag_high_speed +
		flag_buildings + flag_outside_gs)


####################################
# Convert to format suitable for Arc

trips_shp <- smoothed_trips

# Columns not of standard numeric/char type need setting to numeric/character
trips_shp$duration <- as.numeric(trips_shp$duration)
trips_shp$start_time <- as.character(trips_shp$start_time)
trips_shp$end_time <- as.character(trips_shp$end_time)
trips_shp$trip_len <- as.numeric(trips_shp$trip_len)
trips_shp$in_gs_len <- as.numeric(trips_shp$in_gs_len)
trips_shp$building_isect_len <- as.numeric(trips_shp$building_isect_len)
trips_shp$building_isect_pct <- as.numeric(trips_shp$building_isect_pct)


names(trips_shp) <- c("new_trip_id", "device_id", "start_time", "end_time", "orig_id", "geometry",            
	"trip_len", "duration", "path_pts", "speed_mps", "gs_count",            
	"in_gs_len", "pct_in_gs", "start_dist", "end_dist", "in_bld_len", "pct_in_bld",  
	"in_bld_pts", "pct_pts_bd", "home_dist", "short_dist", "incomplete", "high_speed",     
	"in_bldgs", "outside_gs", "flag_count")


# write to shapefile....
st_write(trips_shp, dsn="EXPORT_FILE_PATH.shp")
st_write(trips_shp, dsn="C:/export/trips_final.shp")


