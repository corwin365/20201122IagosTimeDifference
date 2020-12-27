clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IAGOS travel time study, master processing script
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/11/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS 
%all paths and arbitrary choices made in the analysis are set here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%which parts to run?
%this is a list of programme components in order - 1 to run, 0 to skip
%we can save time if earliers parts don't need to be re-run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MasterSettings.Run = [0, ...  %generate airport geolocation dataset
                      0, ...  %find and store all flights between regions
                      1, ...  %rearrange data into routes
                      0, ...  %plot airport metadata
                      0];     %plot flight paths used


%general settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%path to the IAGOS data tree
%this continues all the individual netCDF flight files from the Data Portal
MasterSettings.DataDir = [LocalDataDir,'/IAGOS/Timeseries/'];

%airports to include
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%North America
MasterSettings.Airports.NA = {'ATL','BOS','BWI','CDW','CLE','CLT','CVG','DRM','DTW','EWR', ...
                              'FOK','IAD','JFK','LUK','MKE','MRB','ORD','PHL','PNE','YMX', ...
                              'YQB','YUL','YYR','YYZ','YZD'};
%Europe
MasterSettings.Airports.Eur = {'AGA','AGP','AHO','AMM','AMS','ATH','AYT','BCN','BEY','BOD', ...
                               'BRE','BRU','BTS','BUD','CAI','CDG','CGN','CIA','CRL','DBV', ...
                               'DLM','DME','DRS','DUS','ESB','FCO','FKB','FRA','GHF','GRO', ...
                               'HAJ','HAM','HEL','HER','HSK','IST','LCA','LEI','LEJ','LGW', ...
                               'LHR','LIS','LNZ','LYS','MAD','MAN','MLA','MRS','MUC','MXP', ...
                               'NCE','NUE','ORY','OST','OTP','PMI','PRG','PSA','PUY','RHO', ...
                               'RIX','RLG','SDV','SKG','SNN','SPM','STN','SXB','SZG','SZW', ...
                               'TLS','TLV','TOJ','TXL','UTC','VIE','ZNV','ZQL','ZRH'};

%data cleansing and prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                             

%exclusion zone around start and end airports
%used to avoid unusual behaviour around takeoff and landing due to e.g. 
%air-traffic control
MasterSettings.DepArrExclusion = 200; %km

%time spacing of downsampled full routes, used to generate maps
MasterSettings.ResampleTime = 2; %minutes

%minimum number of flights on a route to use it in our analysis
MasterSettings.MinFlights = 10;

%allowable range of flight times relative to median for route
%this is to exclude unusual flights due to e.g. rerouting
MasterSettings.RelativeTime = [0.90,1.10]; %10% allowance, roughly

%mesh for geographic maps
MasterSettings.Maps.Lon = -100:1:25;
MasterSettings.Maps.Lat =   20:1:80;

%time handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time period we will study
MasterSettings.TimeScale = datenum(1994,1,1):1:datenum(2020,12,31);

%seasons to use
%these don't have to be actual seasons - they could be any arbitrary set of days-of-year
%the programme will use the names of the sub-structures as the "season" names
MasterSettings.Seasons.DJF = date2doy(datenum(2000,12,1):datenum(2001, 3,1)-1);
% % MasterSettings.Seasons.MAM = date2doy(datenum(2000, 3,1):datenum(2000, 6,1)-1);
% % MasterSettings.Seasons.JJA = date2doy(datenum(2000, 6,1):datenum(2000, 9,1)-1);
% % MasterSettings.Seasons.SON = date2doy(datenum(2000, 9,1):datenum(2000,12,1)-1);
MasterSettings.Seasons.All = 1:1:366;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROGRAMME BEGINS HERE. 
%YOU SHOULDN'T NEED TO MODIFY ANYTHING BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. generate geolocation arrays of airports 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(1) == 1;

  %notification
  disp('----------> Geolocating airports')
  
  %set needed variables
  Settings = struct();
  Settings.DataDir = MasterSettings.DataDir;
  
  %call routine
  aa_get_codes(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
else
  
  %notification
  disp('x-x-x-x-x-> Airport geolocation SKIPPED')
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. find and store all flights between the airports 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(2) == 1;

  %notification
  disp('----------> Computing flight metadata')
  
  %set needed variables
  Settings = struct();
  Settings.TimeScale = MasterSettings.TimeScale;
  Settings.DataDir = MasterSettings.DataDir;
  Settings.MinDist = MasterSettings.DepArrExclusion;
  Settings.ResampleTime = MasterSettings.ResampleTime;
  Settings.NA = MasterSettings.Airports.NA;
  Settings.Eur = MasterSettings.Airports.Eur;
  
  %call routine
  bb_prep_data(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Flight metadata computation SKIPPED')
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. identify routes in the data which have enough flights for 
%fair normalisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(3) == 1;

  %notification
  disp('----------> Splitting data into routes ')
  
  %set needed variables
  Settings = struct();
  Settings.MinFlights = MasterSettings.MinFlights;
  Settings.RelativeTime = MasterSettings.RelativeTime;
  Settings.Seasons = MasterSettings.Seasons;
  
  %call routine
  cc_routesplit(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Route splitting SKIPPED')
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. plot airport info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(4) == 1;

  %notification
  disp('----------> Plotting airport info ')
  
  %set needed variables
  Settings = struct();
  Settings.NA = MasterSettings.Airports.NA;
  Settings.Eur = MasterSettings.Airports.Eur;  
  Settings.Seasons = MasterSettings.Seasons;
  
  %call routine
  dd_airportplot(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Airport info plotting SKIPPED')
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. plot route info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(5) == 1;

  %notification
  disp('----------> Plotting route info ')
  
  %set needed variables
  Settings = struct();
  Settings.NA = MasterSettings.Airports.NA;
  Settings.Eur = MasterSettings.Airports.Eur; 
  Settings.Maps = MasterSettings.Maps;
  Settings.Seasons = MasterSettings.Seasons;
  
  %call routine
  ee_routeplot(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Route info plotting SKIPPED')
  
end

