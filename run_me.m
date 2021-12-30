clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IAGOS travel time study, top-level processing script, version 2
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/12/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS 
%all paths and arbitrary choices made in the analysis are set here
%these are the bits you want to change!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%which parts to run?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this is a list of programme components in order - 1 to run, 0 to skip
%this is just to save runtime - if you want a full reset, set all to 1.

%two-letter prefixes refer to the actual functions called - the naming
%convention is related to the order I wrote them and has no further meaning

MainSettings.Run = [0, ...  %aa: generate airport geolocation dataset. Doesn't usually need to be run, and output is handled manually and set in the MainSettings.Airports variables below
                    0, ...  %bb: prep flight and airport data. This is the main data prep routine, but always acts on all the data in Paths.AeolusData - i.e. it only needs rerunning if data is added or removed.
                    0, ...  %cc: prepare climate indices, on the same timescale as the flight data from bb
                    0, ...  %cd: choose which indices will be used raw and which deseasonalised in later analyses
                    0, ...  %dd: split data into routes, normalise flight times, and generate paired routes
                    0, ...  %ee: plot time series of planes used and climate indices
                    0, ...  %ff: plot time taken as a function of time and plane
                    1, ...  %gg: relative duration violin plots
                    0];     


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%file paths to the data and climate indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define paths to data storage directories
Paths.AeolusData = [LocalDataDir,'/IAGOS/TimeSeries'];  %contains individual IAGOS flight netCDF files from iagos.fr
Paths.Indices    = [LocalDataDir,'/Miscellany/'];       %file formats vary - see "generate indices" script for sources
Paths.StoreDir   = ['./data'];                          %storage directory for working files

%define a unique identifier for all source data (i.e. airport metadata, stored flights, and stored climate indices)
%this will be used to SAVE data in routines aa, bb and cc, and to LOAD data in later routines
Paths.SourceIdentifier = 'all'; 

%define a unique identifier to be applied to all postprocessed data (e.g. subsetted data for a particular time period)
%this will be used to SAVE data in routines dd and ee, and to LOAD data in later routines
Paths.PPIdentifier = 'aaa'; 

      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                             
%data filtering and cleaning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

%exclusion zone around start and end airports
%used to avoid unusual behaviour around takeoff and landing due to e.g. 
%air-traffic control
MainSettings.Filtering.DepArrExclusion = 200; %km

%minimum number of flights on a route to use it in our analysis
MainSettings.Filtering.MinFlights = 10;

%allowable range of flight times relative to median for route
%this is to exclude unusual flights due to e.g. rerouting
%these are only applied to SEASONAL data - analyses at the all-data level will include outliers
MainSettings.Filtering.RelativeTime = [0.85,1.15]; 

%deseasonalise flights?
MainSettings.Filtering.DeSeasFlights = 0; %0 no, 1 yes
MainSettings.Filtering.DeSeasPeriod  = 31; %if yes above, how many days should we use as a seasonal averaging window? This will be rounded UP to give whole days at each end of the window

%what is the maximum time difference between two flights paired as a return?
MainSettings.Filtering.MaxDt =2; %days

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gridscales to be used for mapping and comparison to geophysics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time spacing of downsampled full routes, used to generate maps
MainSettings.Maps.ResampleTime = 5; %minutes

%mesh for geographic maps and u/v/T averaging
MainSettings.Maps.Lon = -100:1:25;
MainSettings.Maps.Lat =   20:1:80;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time period we will study. Anything outside this will be discarded.
MainSettings.TimeRange= [datenum(1994,8,3),datenum(2019,12,31)];

%"seasons" to use
%these don't have to be actual seasons - they could be any arbitrary set of days-of-year, and can overlap
%the programme will use the names of the sub-structures as the "season" names
%
%reserved names not to be used: 'tRel','Route','Used','Seasons', 'Indices' (used as variables inside same structs later)
MainSettings.Seasons.DJF = date2doy(datenum(2000,12,1):datenum(2001, 3,1)-1);
MainSettings.Seasons.MAM = date2doy(datenum(2000, 3,1):datenum(2000, 6,1)-1);
MainSettings.Seasons.JJA = date2doy(datenum(2000, 6,1):datenum(2000, 9,1)-1);
MainSettings.Seasons.SON = date2doy(datenum(2000, 9,1):datenum(2000,12,1)-1);
% MainSettings.Seasons.All = 1:1:366;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
%other settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%what climate indices to use?
%data will be plotted and multilinear regressed about these
%all will be normalised to a range of -1 to 1 (except HadCRUT)
MainSettings.Indices = {'ENSO','QBO','HadCRUT','NAM','TSI','SeaIce'};

%what range should we normalise the indices over? (percentiles)
MainSettings.IndexRange = [0,100];

%how far should we smooth "background" data for deseasonalisation?
MainSettings.DSSmooth = 29; %days - must be an odd positive integer

%and which indices should be use deseasonalised, as opposed to raw?
MainSettings.DSIndices = {'SeaIce'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%airports to include. Those with too few flights will be discarded later.
%These were selected by hand, using maps generated from the output of
%routine aa_get_codes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%North America
Airports.NA = {'ATL','BOS','BWI','CDW','CLE','CLT','CVG','DRM','DTW','EWR', ...
               'FOK','IAD','JFK','LUK','MKE','MRB','ORD','PHL','PNE','YMX', ...
               'YQB','YUL','YYR','YYZ','YZD'};
%Europe
Airports.Eur = {'AGA','AGP','AHO','AMM','AMS','ATH','AYT','BCN','BEY','BOD', ...
                'BRE','BRU','BTS','BUD','CAI','CDG','CGN','CIA','CRL','DBV', ...
                'DLM','DME','DRS','DUS','ESB','FCO','FKB','FRA','GHF','GRO', ...
                'HAJ','HAM','HEL','HER','HSK','IST','LCA','LEI','LEJ','LGW', ...
                'LHR','LIS','LNZ','LYS','MAD','MAN','MLA','MRS','MUC','MXP', ...
                'NCE','NUE','ORY','OST','OTP','PMI','PRG','PSA','PUY','RHO', ...
                'RIX','RLG','SDV','SKG','SNN','SPM','STN','SXB','SZG','SZW', ...
                'TLS','TLV','TOJ','TXL','UTC','VIE','ZNV','ZQL','ZRH'};   

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first, save the settings we used so we can refer back to them later if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([Paths.StoreDir,'/settings_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat'],'Paths','MainSettings','Airports')

%ok, let's go...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AA. generate geolocation arrays of airports 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MainSettings.Run(1) == 1;                        
  disp('----------> Geolocating airports')            %notification
  aa_get_codes(Paths)                                 %call routine
  clearvars -except MainSettings Paths Airports       %tidy up
else; disp('x-x-x-x-x-> Airport geolocation SKIPPED'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BB. find and store all flights between the airports 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MainSettings.Run(2) == 1;                        
  disp('----------> Loading and storing flights')     %notification

  %set needed variables
  Settings = struct();
  Settings.MinDist = MainSettings.Filtering.DepArrExclusion;
  Settings.ResampleTime = MainSettings.Maps.ResampleTime;

  bb_prep_data(Paths,Airports,Settings);              %call routine
  clearvars -except MainSettings Paths Airports       %tidy up
else; disp('x-x-x-x-x-> Loading and storing flights SKIPPED'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CC. find and prepare climate indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MainSettings.Run(3) == 1;                        
  disp('----------> Preparing climate indices')     %notification

  %set needed variables
  Settings = struct();
  Indices = MainSettings.Indices;
  IndexRange = MainSettings.IndexRange;
  Settings.DSSmooth = MainSettings.DSSmooth;
  

  cc_generateindices(Paths,Indices,Settings,IndexRange);  %call routine
  clearvars -except MainSettings Paths Airports       %tidy up
else; disp('x-x-x-x-x-> Preparing climate indices SKIPPED'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CD. choose which indices should be raw and which deseasonalised
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MainSettings.Run(4) == 1;                        
  disp('----------> Choosing index time subsetting')  %notification
  Settings.Indices  = MainSettings.Indices;
  Settings.ChosenDS = MainSettings.DSIndices;
  cd_chooseindices(Paths,Settings)         %call routine
  clearvars -except MainSettings Paths Airports                      %tidy up
else; disp('x-x-x-x-x-> Choosing index time subsetting SKIPPED'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DD: split data into routes, normalise flight times, and generate paired routes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MainSettings.Run(5) == 1;          

  disp('----------> Splitting data into routes and seasons')     %notification

  %set needed variables
  Settings = struct();
  Settings.TimeRange     = MainSettings.TimeRange;
  Settings.Seasons       = MainSettings.Seasons;
  Settings.MinFlights    = MainSettings.Filtering.MinFlights;
  Settings.RelativeTime  = MainSettings.Filtering.RelativeTime;
  Settings.DeseasFlights = MainSettings.Filtering.DeSeasFlights;
  Settings.DeseasPeriod  = MainSettings.Filtering.DeSeasPeriod;
  Settings.MaxDt         = MainSettings.Filtering.MaxDt;
  
  dd_routesplit(Paths,Settings,Airports);             %call routine
  clearvars -except MainSettings Paths Airports       %tidy up


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EE. plot time series of planes used and indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MainSettings.Run(6) == 1;                        
  disp('----------> Plotting time series of planes and indices')  %notification
  ee_planes_and_indices(Paths,MainSettings.Indices)                %call routine
  clearvars -except MainSettings Paths Airports                    %tidy up
else; disp('x-x-x-x-x-> Plotting time series of planes and indices SKIPPED'); end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FF. plot time taken as a function of time and plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MainSettings.Run(7) == 1;                        
  disp('----------> Plotting time series of flight times by plane')  %notification
  ff_timetaken(Paths,MainSettings.Filtering.DeSeasFlights)         %call routine
  clearvars -except MainSettings Paths Airports                      %tidy up
else; disp('x-x-x-x-x-> Plotting time series of flight times by plane SKIPPED'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GG. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MainSettings.Run(8) == 1;                        
  disp('----------> Plotting ')  %notification
  disp('xpercent')
  gg_times(Paths,MainSettings.Indices,MainSettings.Seasons, ...
                 MainSettings.Filtering.RelativeTime)                %call routine
  clearvars -except MainSettings Paths Airports                    %tidy up
else; disp('x-x-x-x-x-> Plotting SKIPPED'); end