clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IAGOS travel time study, master processing script
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/11/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS 
%all paths and arbitrary choices made in the analysis are set here
%these are the bits you want to change!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%which parts to run?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this is a list of programme components in order - 1 to run, 0 to skip
%this is just to save runtime - if you want a full reset, set all to 1.
%IF IN DOUBT, RUN THEM ALL, as there may be some unexpected dependencies
%two-letter prefixes refer to the actual functions called - the naming
%convention is related to the order I wrote them and has no further meaning

MasterSettings.Run = [0, ...  %aa: generate airport geolocation dataset
                      0, ...  %bb: find and store all flights between regions
                      0, ...  %cc: rearrange data into routes
                      1, ...  %cd: compute and retain relative flight times
                      0, ...  %dd: plot airport metadata
                      0, ...  %ee: plot flight paths used
                      1, ...  %ff: prepare climate indices, plus deseasonalise and/or delinearise if requested below
                      0, ...  %gg: do and plot multilinear regression, one-way
                      0, ...  %gh: do and plot multilinear regression, round-trip
                      0, ...  %hh: do and plot relative histograms, one-way
                      0, ...  %hi: do and plot relative box plots, one-way                      
                      0, ...  %ii: do and plot relative histograms, round-trip                  
                      0, ...  %ij: do and plot relative boxplots, round-trip
                      1, ...  %jj: time series of planes and indices
                      0, ...  %kk: time series of relative time taken
                      0];     %ll: comparisons of flight time to u, v, T

                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%airports to include. Those with too few flights will be discarded later
%These were selected by hand, using maps generated from the output of
%routine aa.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                             
%data cleansing and prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

%exclusion zone around start and end airports
%used to avoid unusual behaviour around takeoff and landing due to e.g. 
%air-traffic control
MasterSettings.DepArrExclusion = 200; %km

%minimum number of flights on a route to use it in our analysis
MasterSettings.MinFlights = 10;

%allowable range of flight times relative to median for route
%this is to exclude unusual flights due to e.g. rerouting
MasterSettings.RelativeTime = [0.85,1.15]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time spacing of downsampled full routes, used to generate maps
MasterSettings.ResampleTime = 1; %minutes

%mesh for geographic maps and u/v/T averaging
MasterSettings.Maps.Lon = -100:1:25;
MasterSettings.Maps.Lat =   20:1:80;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time period we will study
MasterSettings.TimeScale = datenum(1994,1,1):1:datenum(2019,12,31);

%seasons to use
%these don't have to be actual seasons - they could be any arbitrary set of days-of-year
%the programme will use the names of the sub-structures as the "season" names
MasterSettings.Seasons.DJF = date2doy(datenum(2000,12,1):datenum(2001, 3,1)-1);
MasterSettings.Seasons.MAM = date2doy(datenum(2000, 3,1):datenum(2000, 6,1)-1);
MasterSettings.Seasons.JJA = date2doy(datenum(2000, 6,1):datenum(2000, 9,1)-1);
MasterSettings.Seasons.SON = date2doy(datenum(2000, 9,1):datenum(2000,12,1)-1);
MasterSettings.Seasons.All = 1:1:366;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
%other settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path to the IAGOS data tree
%this contains all the individual netCDF flight files from the IAGOS Data Portal
MasterSettings.DataDir = [LocalDataDir,'/IAGOS/Timeseries/'];

%what climate indices to use?
%data will be plotted and multilinear regressed about these
%all will be normalised to a range of -1 to 1 (except HadCRUT)
MasterSettings.Indices = {'ENSO','QBO','HadCRUT','NAM','TSI','Time'};

%what fraction of the data should be used for index-comparison histograms?
MasterSettings.IndexFraction = 0.2; %1 = all data

%how many histogram bins to use, spread across the range of valid relative
%times, and how many bins to smooth output by?
MasterSettings.IndexHistBins   = 40;
MasterSettings.IndexHistSmooth = 3; %must be odd

%do we want to lag the indices for the regression, and if so over how long
%a window and with what step size? 
%this goes combinatorically with the list of indices above - so be careful,
%the runtime can get very large very fast if this is used
MasterSettings.Reg.Lag   = 0; %1 for yes, 0 for no
MasterSettings.Reg.Steps = [-60,-30,-10,-5,-2,0]; %values to try

%what range should we compute ROUND TRIP histograms over
%relative to a single trip
MasterSettings.RTRelativeTime = [1.9,2.1]; 

%should we deseasonalise the data and indices?
%this will not deseasonalise: QBO, ENSO, HadCRUT, TSI, Time
MasterSettings.Deseasonalise = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROGRAMME BEGINS HERE. 
%YOU SHOULDN'T NEED TO MODIFY ANYTHING BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% 4. compute and store relative flight times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(4) == 1;

  %notification
  disp('----------> Computing relative flight times and pairing flights')
  
  %set needed variables
  Settings = struct();
  Settings.NA         = MasterSettings.Airports.NA;
  Settings.Eur        = MasterSettings.Airports.Eur; 
  Settings.Frac       = MasterSettings.IndexFraction;
  Settings.Seasons    = MasterSettings.Seasons;
  
  %call routine
  cd_timecompute(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Relative flight time computation and pairing SKIPPED')
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. plot airport info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(5) == 1;

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
%% 6. plot route info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(6) == 1;

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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. prepare climate indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(7) == 1;

  %notification
  disp('----------> Preparing climate indices ')
  
  %set needed variables
  Settings = struct();
  Settings.Indices   = MasterSettings.Indices;
  Settings.TimeScale = MasterSettings.TimeScale;
  Settings.Reg = MasterSettings.Reg;
  Settings.DS = MasterSettings.Deseasonalise;
  
  
  %call routine
  ff_generateindices(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Climate indices preparation SKIPPED')
  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. regress data against climate indices, one-way
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(8) == 1;

  %notification
  disp('----------> Regressing against climate indices, one-way ')
  
  %set needed variables
  Settings = struct();
  Settings.Indices = MasterSettings.Indices;
  Settings.Seasons = MasterSettings.Seasons;
  Settings.NA      = MasterSettings.Airports.NA;
  Settings.Eur     = MasterSettings.Airports.Eur; 
  Settings.Reg     = MasterSettings.Reg;
  
  %call routine
  gg_regression(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Regression against climate indices, one-way SKIPPED')
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. regress data against climate indices, round-trip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(9) == 1;

  %notification
  disp('----------> Regressing against climate indices, round-trip ')
  
  %set needed variables
  Settings = struct();
  Settings.Indices = MasterSettings.Indices;
  Settings.Seasons = MasterSettings.Seasons;
  Settings.NA      = MasterSettings.Airports.NA;
  Settings.Eur     = MasterSettings.Airports.Eur; 
  Settings.Reg     = MasterSettings.Reg;
  
  %call routine
  gh_regression2(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Regression against climate indices, round-trip SKIPPED')
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10. difference histograms, one-way
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(10) == 1;

  %notification
  disp('----------> Producing one-way distance histograms ')
  
  %set needed variables
  Settings = struct();
  Settings.Indices    = MasterSettings.Indices;
  Settings.Seasons    = MasterSettings.Seasons;
  Settings.NA         = MasterSettings.Airports.NA;
  Settings.Eur        = MasterSettings.Airports.Eur; 
  Settings.Frac       = MasterSettings.IndexFraction;
  Settings.HistBins   = linspace(MasterSettings.RelativeTime(1), ...
                                 MasterSettings.RelativeTime(2), ...
                                 MasterSettings.IndexHistBins+1);
  Settings.HistSmooth = MasterSettings.IndexHistSmooth;
  
  %call routine
  hh_histograms(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> One-way distance histograms SKIPPED')
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11. difference box plots, one-way
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(11) == 1;

  %notification
  disp('----------> Producing one-way distance  box plots ')
  
  %set needed variables
  Settings = struct();
  Settings.Indices    = MasterSettings.Indices;
  Settings.Seasons    = MasterSettings.Seasons;
  Settings.NA         = MasterSettings.Airports.NA;
  Settings.Eur        = MasterSettings.Airports.Eur; 
  Settings.Frac       = MasterSettings.IndexFraction;
  Settings.HistBins   = linspace(MasterSettings.RelativeTime(1), ...
                                 MasterSettings.RelativeTime(2), ...
                                 MasterSettings.IndexHistBins+1);
  Settings.HistSmooth = MasterSettings.IndexHistSmooth;
  
  %call routine
  hi_boxplots(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> One-way distance box plots SKIPPED')
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 12. difference histograms, round-trip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(12) == 1;

  %notification
  disp('----------> Producing round-trip distance histograms ')
  
  %set needed variables
  Settings = struct();
  Settings.Indices    = MasterSettings.Indices;
  Settings.Seasons    = MasterSettings.Seasons;
  Settings.NA         = MasterSettings.Airports.NA;
  Settings.Eur        = MasterSettings.Airports.Eur; 
  Settings.Frac       = MasterSettings.IndexFraction; 
  Settings.HistBins   = linspace(MasterSettings.RTRelativeTime(1), ...
                                 MasterSettings.RTRelativeTime(2), ...
                                 MasterSettings.IndexHistBins+1);
  Settings.HistSmooth = MasterSettings.IndexHistSmooth;
  
  %call routine
  ii_histograms2(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Round-trip distance histograms SKIPPED')
  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 13. difference histograms, round-trip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(13) == 1;

  %notification
  disp('----------> Producing round-trip distance boxplots ')
  
  %set needed variables
  Settings = struct();
  Settings.Indices    = MasterSettings.Indices;
  Settings.Seasons    = MasterSettings.Seasons;
  Settings.NA         = MasterSettings.Airports.NA;
  Settings.Eur        = MasterSettings.Airports.Eur; 
  Settings.Frac       = MasterSettings.IndexFraction; 
  Settings.HistBins   = linspace(MasterSettings.RTRelativeTime(1), ...
                                 MasterSettings.RTRelativeTime(2), ...
                                 MasterSettings.IndexHistBins+1);
  Settings.HistSmooth = MasterSettings.IndexHistSmooth;
  
  %call routine
  ij_boxplots2(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Round-trip distance boxplots SKIPPED')
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 14. metadata about individual aircraft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(14) == 1;

  %notification
  disp('----------> Producing time series of planes and indices ')

  %call routine
  Settings.Frac = MasterSettings.IndexFraction;
  jj_planes_and_indices(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Time series of planes and indices plots SKIPPED')
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 15. time series of relative time taken
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(15) == 1;

  %notification
  disp('----------> Producing raw time series ')
  
  %set needed variables
  Settings = struct();
  Settings.NA         = MasterSettings.Airports.NA;
  Settings.Eur        = MasterSettings.Airports.Eur; 
  Settings.Frac       = MasterSettings.IndexFraction;
  Settings.Seasons    = MasterSettings.Seasons;
  
  %call routine
  kk_timetaken(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> Raw time series SKIPPED')
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 16. compare flight time to u,v,T encountered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MasterSettings.Run(16) == 1;

  %notification
  disp('----------> Comparing to u,v,T ')
  
  %set needed variables
  Settings = struct();
  Settings.Seasons    = MasterSettings.Seasons;
  
  %call routine
  ll_comparison(Settings)
  
  %tidy up
  clearvars -except MasterSettings
  
  else
  
  %notification
  disp('x-x-x-x-x-> u,v,T comparison SKIPPED')
  
end
