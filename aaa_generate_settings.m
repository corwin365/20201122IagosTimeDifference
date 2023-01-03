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
SettingsID = 'test567'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%file paths to the data and climate indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define paths to data storage directories
Paths.AeolusData = [LocalDataDir,'/IAGOS/TimeSeries'];   %contains individual IAGOS flight netCDF files from iagos.fr
Paths.EdData     = [LocalDataDir,'/corwin/ed_aircraft']; %arr/dep data supplied by Ed Gryspeerdt
Paths.Indices    = [LocalDataDir,'/Miscellany/'];        %file formats vary - see "generate indices" script for sources
Paths.DataDir    = './data/';                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%which datasets should we use?
Choices.DataSets = [1,2]; %1 is IAGOS, 2 is the data from Ed Gryspeerdt

%what time period are we looking over?
Choices.TimeRange = [datenum(2005,1,1),datenum(2007,12,31)];

%for data where we have the flight trace, how close to the airport should we discard data?
Choices.MinDist = 100; %km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%climate index handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%what indices should we use?
Indices.List = {'ENSO','NAO','QBO','SSTs','TSI'};

%how far should we smooth "background" data for index deseasonalisation?
Indices.DSSmooth = 61; %days - must be an odd positive integer

%what range should we normalise the indices over? (percentiles)
Indices.IndexRange = [0,100];

%which indices should be use deseasonalised, as opposed to raw?
Indices.DS = {'SSTs','SeaIce'};

%which indices should be delinearised?
Indices.DL = {'SSTs','SeaIce'};

%how many days should we smooth the climate indices by? Set to 0 to not smooth
%this is different to the above - the above is just a window for deaseasonalisation and is REMOVED from the data 
%this is an overall smoothing APPLIED TO the data 
Indices.SmoothLength = 7;

%which indices should NOT be smoothed?
Indices.DoNotSmooth = {'NAO'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%airports to include. Those with too few flights will be discarded later.
%These were selected by hand, using maps generated from the output of
%routine aa_get_codes for the IAGOS set and then manually from a list
%generated from the files for the additional airports added with Ed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%North America
Airports.NA = {'ATL','BOS','BWI','CDW','CLE','CLT','CVG','DTW','EWR','FOK', ...
               'IAD','JFK','LUK','MKE','MRB','ORD','PHL','PNE','YMX','YQB', ...
               'YUL','YYR','YYZ', ...
               ... %below added after additional data supplied by Ed Gryspeerdt, not present in IAGOS
               'ALC','BER','BHX','DFW','IAH','MCO','MIA','PBI','RSW','SAV', ...
               'SFB','TPA','YHZ','YYT','CMH','IND','MCI','MDW','MEM','MSP', ...
               'MSY','PVD','RDU','RIC','SDF','TUL','YOW','BNA','BUF','DCA', ...
               'JAX','PIT','STL','SYR','YWG'};
%Europe
Airports.Eur = {'AGA','AGP','AHO','AMM','AMS','ATH','AYT','BCN','BEY','BOD', ...
                'BRE','BRU','BTS','BUD','CAI','CDG','CGN','CIA','CRL','DBV', ...
                'DLM','DME','DRS','DUS','ESB','FCO','FKB','FRA','GHF','GRO', ...
                'HAJ','HAM','HEL','HER','HSK','IST','LCA','LEI','LEJ','LGW', ...
                'LHR','LIS','LNZ','LYS','MAD','MAN','MLA','MRS','MUC','MXP', ...
                'NCE','NUE','ORY','OST','OTP','PMI','PRG','PSA','PUY','RHO', ...
                'RIX','RLG','SDV','SKG','SNN','SPM','STN','SXB','SZG','SZW', ...
                'TLS','TLV','TOJ','TXL','UTC','VIE','ZRH', ...
                ... %below added after additional data supplied by Ed Gryspeerdt, not present in IAGOS
                'BLL','BLQ','CPH','DUB','EDI','FAO','GLA','GVA','LTN','LUX', ...
                'OPO','SCQ','STR','VRN','BFS','EIN','LGA','NAP','TRN','VCE', ...
                'BDS','BGY','CTA','IBZ','PMO','BRI','CAG'};  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save!
ID = SettingsID;
save([Paths.DataDir,SettingsID,'.mat'],'Paths','Airports','Choices','ID','Indices')

