clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Impact of climate-dynamical processes on flight times
%
%script to run analyses
%
%Corwin Wright, c.wright@bath.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%what settings are we using?
% SettingsID = 'basic_annual';
% SettingsID = 'basic_noannual';
SettingsID = 'honolulu_annual';

%load them
Settings = load(['data/',SettingsID,'.mat']);
clear SettingsID  %it's contained in the file we loaded already

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% import and generate raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %generate airport geolocation
% %settings choices locked in here: Choices.Airports
% generate_geolocation(Settings);
% 
% 
% %generate IAGOS data
% %settings choices locked in here: Choices.TimeRange, Choices.MinDist
% if sum(ismember(Settings.Choices.DataSets,1)) | sum(ismember(Settings.Choices.DataSets,3)); import_iagos_data(Settings); end
% 
% %generate data from Ed?
% %settings choices locked in here: [none]
% if sum(ismember(Settings.Choices.DataSets,2)); import_ed_data(   Settings); end
% 
% %generate merged dataset
% %settings choices locked in by here: Choices.DataSets
% generate_merged_data(Settings);

% %extract full flight traces
% %settings choices locked in by here: [none]
% generate_track_data(Settings);

% %generate ERA5 wind data
% %settings choices locked in by here: Choices.WindMap
% generate_wind_data(Settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalise and prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %split data into routes
% %settings choices locked in here: [none]
% process_routesplit(Settings);
% 
% %pair routes to produce round trips
% %settings choices locked in here: Choices.MinPairDistance, Choices.Directions
% process_roundtrips(Settings);
% 
% %split into seasons, and normalise flight times
% %settings choices locked in here: Choices.MinFlights, Choices.MaxDeviation, Seasons
% process_seasonsplit(Settings);
% 
% %generate climate indices for the data
% %must be done after all filtering and merging, so arrays match flight info
% %settings choices locked in here: Indices
% process_climate_indices_v3(Settings);
% 
% %compute optimal lags for each index
% %settings choices locked in here: [none]
% process_optimal_lags(Settings);

% %cost of each minute of delay, in both CO2 and USD, and total mber of flights over the Atlantic
% %settings choices locked in here: [none]
% process_scalefactors(Settings);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% carry out and plot analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% data description and assessment
%  (sensitivity tests are handled 
%  separately, outside this script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %climate index independence tests
% %settings choices locked in here: [none]
% analysis_index_stat_tests(Settings)

% % %data coverage and climate indices
% % %settings choices locked in here: [none]
% analysis_coverageandindices(Settings)
% 
% %linear trend analysis by season
% %settings choices locked in here: [none]
% analysis_linear(Settings)

% %linear trend analysis by airline
% %settings choices locked in here: [none]
% analysis_linear_byairline(Settings)
% 
% %number of flights on each route
% %settings choices locked in here: [none], but the programme does contain hardcoded map boundaries
% analysis_routestats(Settings)
% 
% %all flights map
% %settings choices locked in here: Choices.FlightPathMap
% analysis_maps(Settings)

%% data analysis - KDFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%index-split KDFs
%settings choices locked in here: Choices.KDFSplit
% analysis_indexsplit(Settings);


%% data analysis - regressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % linear regression
% %settings choices locked in here: [none]
analysis_regression(Settings);

% %cost implications, in CO2 and USD
% % %settings choices locked in here: [none]
% %requires analysis_regression to have been run to generate necessary coefficients
% analysis_cost_histo(Settings);
% analysis_cost(Settings);


%% data analysis - ERA5 comparisons (not currently working)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % ERA5 comparison of wind speeds
% %settings choices locked in here: [none]
% analysis_winddiff(Settings);

%% data analysis - cruise height impacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cruise height analysis
%settings choices locked in here: Choices.ZScale
% analysis_cruiseheight(Settings);

% % cruise pressure analysis
% %settings choices locked in here: Choices.PScale
% analysis_cruisepressure(Settings);

% % check if the tropopause is correlated with flight level
% %settings choices locked in here: [none]
% analysis_cruise_v_tp(Settings);





% analysis_frac_strat_trop(Settings);