clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Impact of climate-dynamical processes on flight times
%
%script to generate settings files
%
%Corwin Wright, c.wright@bath.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% name of settings file to generate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%name of this set of analyses
SettingsID = 'basic_noannual';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%file paths to the data and climate indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define paths to data storage directories
Paths.AeolusData = [LocalDataDir,'/IAGOS/TimeSeries'];   %contains individual IAGOS flight netCDF files from iagos.fr
Paths.EdData     = [LocalDataDir,'/corwin/ed_aircraft']; %arr/dep data supplied by Ed Gryspeerdt
Paths.Indices    = [LocalDataDir,'/Miscellany/'];        %file formats vary - see "generate indices" script for sources
Paths.DataDir    = './data/';                          
Paths.Era5Dir    = [LocalDataDir,'/ERA5/'];              %path to ERA5 data    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%which datasets should we use?
%1 is IAGOS
%2 is the data from Ed Gryspeerdt. 
%3 is the IAGOS data subsetted to just days in Ed's dataset
% duplicates are removed when the data is imported, so it's safe to combine multiple similar options
Choices.DataSets = [1]; 

%what time period are we looking over?
Choices.TimeRange = [728295,739340];

%what flight directions do we want to proces results for? E eastwards, W westwards, R round trip
Choices.Directions = {'W','R','E'};

%for data where we have the flight trace, how close to the airport should we discard data?
Choices.MinDist = 10; %km. Selected via sensitivity testing over range 0-1000km

%what is the minimum number of flights (per route-season) to be included in the dataset?
Choices.MinFlights = 10;

%how close in time do flights have to be to be paired as a round-trip?
Choices.MinPairDistance = 1; %days

%what is the maximum distance from the (route-season) median travel time before which the flight will be discarded?
Choices.MaxDeviation = 20; %percent. Set to 100 or greater to not filter

%how should we define the airports used
Choices.Airports = 'list'; %'list' for whitelist definition, or 'geo' for geographic definition

%how large a gap in an IAGOS data record is acceptable before we discard it?
Choices.Maxdt = 15*60; %seconds
Choices.MaxdLat = 10; %degrees. Need to be generous as some flights get near the pole.
Choices.MaxdLon = 10; %degrees

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%how should we split the data into "seasons"?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%"seasons" to use
%these don't have to be actual seasons - they could be any arbitrary set of days-of-year, and can overlap
%the programme will use the names of the sub-structures as the "season" names
Seasons.All = 1:1:366;
Seasons.DJF = date2doy(datenum(2000,12,1):datenum(2001, 3,1)-1);
Seasons.MAM = date2doy(datenum(2000, 3,1):datenum(2000, 6,1)-1);
Seasons.JJA = date2doy(datenum(2000, 6,1):datenum(2000, 9,1)-1);
Seasons.SON = date2doy(datenum(2000, 9,1):datenum(2000,12,1)-1);

Seasons.List = fieldnames(Seasons);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%climate index handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%do we want to lag the data? 
%(optimal lags will be computed anyway, this determines if they are applied)
Choices.ApplyLags = 1;

%what is the maximum lag to test out to?
Choices.MaxLag = 365; % days

%what indices should we use?
% Indices.List = {'Annual','ENSO','NAO','QBO','Time','TSI'};
Indices.List = {'ENSO','NAO','QBO','Time','TSI'};

%what range should we normalise the indices over? (percentiles)
Indices.IndexRange = [0,100];

%which indices should be used deseasonalised, as opposed to raw?
Indices.DS = {'SSTs','SeaIce'}; %any not specified in Indices.List will be ignored

%which indices should be delinearised?
Indices.DL = {'SSTs','SeaIce','AMO'};  %any not specified in Indices.List will be ignored

%how far should we smooth "background" data for index deseasonalisation, if requested.
Indices.DSSmooth = 61; %days - must be an odd positive integer

%how many days should we smooth the climate indices by? Set to 0 to not smooth
%this is different to the above:
% - the above is just a window for deaseasonalisation and is REMOVED from the data 
% - this is an overall smoothing APPLIED TO the data 
Indices.SmoothLength = 7;

%which indices should NOT be smoothed?
Indices.DoNotSmooth = {'NAO'};

%colours to use for the indices in plots
Indices.Colours.ENSO    = [ 57,159,228]./255;
Indices.Colours.Fuel    = [  0,  0,  0]./255;
Indices.Colours.HadCRUT = [255,178,102]./255;
Indices.Colours.NAO     = [ 46,148,130]./255;
Indices.Colours.NAM     = [ 69,174, 98]./255;
Indices.Colours.QBO     = [255,209,107]./255;
Indices.Colours.SeaIce  = [152, 51, 91]./255;
Indices.Colours.u1060   = [255,102,178]./255;
Indices.Colours.TSI     = [196, 66, 79]./255;
Indices.Colours.SSTs    = [204,204,0]./255;
Indices.Colours.Time    = [1,1,1].*0.6;
Indices.Colours.Annual  = [113, 69,168]./255;
Indices.Colours.AMO     = [255,102,178]./255;

%marker symbols to use for the indices in plots
Indices.Symbols.ENSO    = 's';
Indices.Symbols.Fuel    = 'o';
Indices.Symbols.HadCRUT = 'o';
Indices.Symbols.NAO     = 'v';
Indices.Symbols.NAM     = 'o';
Indices.Symbols.QBO     = '^';
Indices.Symbols.SeaIce  = 'o';
Indices.Symbols.TSI     = 's';
Indices.Symbols.SSTs    = 'o';
Indices.Symbols.Time    = 'd';
Indices.Symbols.Annual  = 'o';
Indices.Symbols.AMO     = 'o';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%airports to include. We'll provide two definitions, one a whitelist of
%specific airports and the other a geographic definition, which can  be 
%switched around with a flag above in the main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%whitelist approach. These were selected by hand, using maps generated from 
%the output of routine aa_get_codes for the IAGOS set and then manually from 
%a list generated from the files for the additional airports added with Ed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%North America
Airports.List.NA = {'ATL','BOS','BWI','CDW','CLE','CLT','CVG','DTW','EWR','FOK', ...
                    'IAD','JFK','LUK','MKE','MRB','ORD','PHL','PNE','YMX','YQB', ...
                    'YUL','YYR','YYZ', ...
                    ... %below added after additional data supplied by Ed Gryspeerdt, not present in IAGOS
                    'DFW','IAH','MCO','MIA','PBI','RSW','SAV', ...
                    'SFB','TPA','YHZ','YYT','CMH','IND','MCI','MDW','MEM','MSP', ...
                    'MSY','PVD','RDU','RIC','SDF','TUL','YOW','BNA','BUF','DCA', ...
                    'JAX','PIT','STL','SYR','YWG','LGA'};
%Europe
Airports.List.Eur = {'AGA','AGP','AHO','AMM','AMS','ATH','AYT','BCN','BEY','BOD', ...
                     'BRE','BRU','BTS','BUD','CAI','CDG','CGN','CIA','CRL','DBV', ...
                     'DLM','DME','DRS','DUS','ESB','FCO','FKB','FRA','GHF','GRO', ...
                     'HAJ','HAM','HEL','HER','HSK','IST','LCA','LEI','LEJ','LGW', ...
                     'LHR','LIS','LNZ','LYS','MAD','MAN','MLA','MRS','MUC','MXP', ...
                     'NCE','NUE','ORY','OST','OTP','PMI','PRG','PSA','PUY','RHO', ...
                     'RIX','RLG','SDV','SKG','SNN','SPM','STN','SXB','SZG','SZW', ...
                     'TLS','TLV','TOJ','TXL','UTC','VIE','ZRH', ...
                     ... %below added after additional data supplied by Ed Gryspeerdt, not present in IAGOS
                     'BLL','BLQ','CPH','DUB','EDI','FAO','GLA','GVA','LTN','LUX', ...
                     'OPO','SCQ','STR','VRN','BFS','EIN','NAP','TRN','VCE', ...
                     'BDS','BGY','CTA','IBZ','PMO','BRI','CAG','BER','ALC','BHX'};  

                 
%geographic approach.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Airports.Bounds.NA  = [-95,-65,30,50]; %lonmin,lonmax,latmin,latmax
Airports.Bounds.Eur = [-10, 20,35,60]; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analysis-specific options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ERA5 wind map resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Choices.WindMap.LonScale   = -120:5:10;
Choices.WindMap.LatScale   = 30:5:70;
Choices.WindMap.TimeScale  = min(Choices.TimeRange):1:max(Choices.TimeRange);


%flight path map choices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%map gridding
Choices.FlightPathMap.LonScale = -98:0.2:18;
Choices.FlightPathMap.LatScale = 24:0.2:72;

%to save time, only regenerate the data if a file containing them doesn't exist
%regeneration can be forced by deleting the store file, or overriding this flag manually in aaa_run_me
Choices.FlightPathMap.Regenerate = 0;


%index-split KDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%what percentage of top and bottom data should we use for the KDF split analysis?
Choices.KDFSplit.CutOff = 20;

%what bins should we use (in minutes delay)
Choices.KDFSplit.Bins = -60:1:60;

%statistical significance threshold
Choices.KDFSplit.Alpha = 0.05;

%summary percentile bands
Choices.KDFPc.Percentiles = [2.5,18,50,82,97.5];

%lagged regression possible lags (days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Choices.LaggedRegression.LagScale = 0:1:60;


%height bins for cruise height/pressure analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%note that, since we take a MODE, the results are sensitive to this choice
Choices.ZScale = 0:0.2:14;
Choices.PScale = 150:10:350;


%fixed date to use in cost computstion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Choices.PriceFixDate = datenum(2023,5,15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save!
ID = SettingsID;
save([Paths.DataDir,SettingsID,'.mat'],'Paths','Airports','Choices','ID','Indices','Seasons')
