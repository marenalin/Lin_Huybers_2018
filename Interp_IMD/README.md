%%Contents

Cos_interp.m is the main function used to interpolate the TRMM data according to the IMD station counts. 

Post_process_figs.m takes all the data generated through interpolation to produce each of the figures and statistics in the manuscript.

Both functions are commented, and please feel free to contact my first name dot last name at gmail.com if you have further questions. 

%%Data 
%%%%%%%%%%
%%%%%%%%%%

IMD_rf provides the IMD 0.25 degree by 0.25 degree data. If this is to be used outside of replicating our results here, it must be purchased from the India Meteorological Department. Both rf25 (daily rainfall amount in mm, interpolated from station data collected at 8:30 am) and stn25 (count of stations within a gridbox) are arranged according to location and time. Locations are given by longitude and latitude coordinates, available as lonr and latr, respectively. Refer to Pai et al. publications for descriptions of these data. Time are given as ttm (fractions of the year), y25, m25, d25, which are the year, month, and day respectively. 

IMD_rf_1x1 provides the IMD 1 degree by 1 degree data. Rainfall and station count are given as ind_rf_cut and ind_stn_cut, respectively. These data are described by the Rajeevan et al. publications. Their dimensions are time and location. Time rows correspond to yy, mm, dd, which are the year, month, and day of each row of data. The locations are latm_cut and lonm_cut, the latitude and longitude. 

TRMM_mic are all available time-points of TRMM microwave-only data over the South Asia region, including India. The structure are arranged as vals (rainfall amount values), yy, mm, dd, hh, which are the year, month, day, and hour of the microwave data value in TRMM_mic(1).vals. The locations are given by lat_TRMM and lon_TRMM. 

See cos_interp.m, the Matlab function used to interpolate the TRMM microwave data using the 0.25 degree station counts, in order to understand the specifics of these data. 


%%TRMM_interp_3am
%%%%%%%%%%
%%%%%%%%%%

%% Cos_interp.m in the main folder is the function used to interpolate each daily field of rainfall data in TRMM_mic according to corresponding arrangement of stations given in IMD_rf. The output of this function is given in the folder TRMM_interp_3am. Each output .mat file is named as TRMM_interp_3am_XXXX, in which the year corresponds to the TRMM year interpolated according to all IMD years. 

%%TRMM_interp_3am_XXXX contains one variable, rf25_interp, which is year XXXX TRMM rainfall interpolated according to the daily summer 1951-2016 station counts from the IMD 0.25 degree data. The location dimension of length 4964 corresponds to the latr and lonr variables of IMD_rf. The dimension 8052 corresponds to the indices l1 of time variables y25, m25, d25, found using lines 96-98 of cos_interp.m:

    l1=find(m25>5 & m25<10);
    l2=find(y25>1950);
    l1=intersect(l1,l2);
 
%%TRMM_interp_1x1
%%%%%%%%%%
%%%%%%%%%%

%%TRMM_interp_cos_1x1_XXXX contains one variable, rf1_interp, which is year XXXX TRMM rainfall interpolated according to the daily summer 1951-2015 station counts from the IMD 1 degree data. The dimension are time by location, where locations are given by lonm_cut and latm_cut in IMD_rf_1x1.mat in the ./Data folder. See Cos_interp_1x1.m and TRMM_1x1_main.m for their provenance and computation. Keep in mind that we have TRMM data at 0.25 degree resolution, but the IMD 1 degree  by 1 degree data and the station counts are at 16x lower resolution. Therefore, treating the 0.25 degree TRMM gridboxes as stations required specifying exact locations when all we had were station counts at the 1 degree resolution. See the Rajeevan et al. texts for the exact interpolation method we followed, in addition to the threshold distances we employed.  


%%Interp_IMD
%%%%%%%%%%
%%%%%%%%%%

%%Each of IMD_interp_rep_n_XXXX takes XXXX year of 0.25 degree IMD rainfall values as truth and reinterpolates it according to the 1951-2016 station counts. These values of course already have the artifacts of interpolation, but this shows how much further interpolation can distort values after taking these fields as truth. 
