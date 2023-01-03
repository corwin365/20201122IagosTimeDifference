function [] = bb_prep_data(Paths,Airports,Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find all flights between the chosen airports in Europe and North America
%then retain their metadata and some basic info (e.g. path, u, v, T)
%
%this routine is almost definitely the dominant runtime sink, and will 
%only need rerunning if new data is added. 
%
%
%additional data checks:
%A. remove any data within 100km of start and end of flight, to avoid 
%looping around the airport before or after departure
%B. modified to retain unique identifiers for individual flights, to 
%identify any issues due to equipment change
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/12/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flight filenames
FlightFiles = wildcardsearch(Paths.AeolusData,'*.nc');


%possible airports, and their locations
Airports.IDs  = [Airports.Eur,Airports.NA];
Airports.Lons = NaN(numel(Airports.IDs),1);
Airports.Lats = Airports.Lons;

%metainfo about flights
Flights.Dep     = cell(numel(FlightFiles),1);  %flight start
Flights.Arr     = cell(numel(FlightFiles),1);  %flight end
Flights.t       = NaN( numel(FlightFiles),1);  %flight time
Flights.Date    = NaN( numel(FlightFiles),1);  %date of flight 
Flights.PlaneID = cell(numel(FlightFiles),1);  %unique aircraft identifier
Flights.InstID  = cell(numel(FlightFiles),1);  %unique instrument identifier

%downsampled flight tracks, for later analysis
Flights.Paths   = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% where are the airports?
%this information isn't used here, but will be needed to generate maps later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Codes = load([Paths.StoreDir,'/airport_codes_',Paths.SourceIdentifier,'.mat']);
for iAirport=1:1:numel(Airports.IDs)
  
  ID = Airports.IDs{iAirport};
  ID = find(contains(Codes.Codes,ID));
  if numel(ID) == 0; continue; end %airport not found in data
  Coords = strsplit(Codes.Coords{ID},' ');
  
  Airports.Lons(iAirport) = str2num(Coords{1});
  Airports.Lats(iAirport) = str2num(Coords{2});  
end
clear Codes iAirport ID Coords

%drop all airports that haven't been found
Good = find(~isnan(Airports.Lons+Airports.Lats));
Airports.IDs  = Airports.IDs( Good);
Airports.Lons = Airports.Lons(Good);
Airports.Lats = Airports.Lats(Good);
clear Good

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, find all the flights between these east-west airport-pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textprogressbar('Processing flight files ')
Count = 0; %count of flights which meet criteria, for indexing purposes
for iFlight = 1:1:numel(FlightFiles)
try
 if mod(iFlight,100) == 0; textprogressbar(iFlight./numel(FlightFiles).*100); end
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % process metadata
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %load the flight 
  Data = rCDF(FlightFiles{iFlight});
  
  %get the file metadata
  MetaData = Data.MetaData.Attributes.Global;  
  
  %hence, identify start and end location, date, and plane/instrument ID
  NFields = numel(MetaData);
  for iField=1:1:NFields
    if strcmp(MetaData(iField).Name,'departure_airport')
      Dep = MetaData(iField).Value; Dep = Dep(1:3);
    elseif strcmp(MetaData(iField).Name,'arrival_airport')
      Arr = MetaData(iField).Value; Arr = Arr(1:3);
    elseif strcmp(MetaData(iField).Name,'departure_UTC_time')
      Date = MetaData(iField).Value;
    elseif strcmp(MetaData(iField).Name,'title')
      ID = MetaData(iField).Value;
    end
  end
  clear MetaData NFields iField
  
  %convert the ID string to an *aircraft* ID and an *instrument* ID
  ID = split(ID,',');
  Plane = ID{4}; Instrument = ID{1};
  Plane(Plane == ' ') = []; Instrument(Instrument == ' ') = []; %remove excess spaces
  clear ID
  
  %convert date to Matlab time
  Date = datenum(Date,'YYYY-mm-DDTHH:MM:ss');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % is this a valid flight?
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %is each airport on one of our lists?
  if ~ismember(Arr,Airports.IDs); continue; end
  if ~ismember(Dep,Airports.IDs); continue; end    
  
  %find which continent the arrival and depature are in
  Start.NA  = ismember(Dep,Airports.NA);
  Start.Eur = ismember(Dep,Airports.Eur);
  End.NA    = ismember(Arr,Airports.NA);
  End.Eur   = ismember(Arr,Airports.Eur);    
    
  %are the two ends on different continents?
  if Start.NA  + End.NA  == 2; continue; end
  if Start.Eur + End.Eur == 2; continue; end
  clear Start End
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % let's clean up the data a bit, for fairness
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  
  %1. remove any points within Settings.MinDist of DEP or ARR
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  distanceFromDep = nph_haversine([ones(numel(Data.lat),1).*Data.lat(1),ones(numel(Data.lat),1).*Data.lon(1)], ...
                                  [Data.lat,Data.lon]);
  distanceFromArr = nph_haversine([ones(numel(Data.lat),1).*Data.lat(end),ones(numel(Data.lat),1).*Data.lon(end)], ...
                                  [Data.lat,Data.lon]);
  Bad = find(distanceFromArr < Settings.MinDist | distanceFromDep < Settings.MinDist);
  
  Data.UTC_time(Bad) = NaN;
  clear Bad distanceFromArr distanceFromDep
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % store flight metadata
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
  Count = Count+1; %valid flight
  
  %store the information we want
  Flights.Dep{    Count} = Dep;
  Flights.Arr{    Count} = Arr;
  Flights.t(      Count) = nanmax(Data.UTC_time) - nanmin(Data.UTC_time); %seconds
  Flights.Date(   Count) = Date;
  Flights.PlaneID{Count} = Plane;
  Flights.InstID{ Count} = Instrument;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % store flight path
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  
  %store the full flight path, downsampled to ~Settings.MinTime
  IAGOS = prep_iagos(FlightFiles{iFlight},'SamplingRate',Settings.ResampleTime./24./60);
  Flights.Paths.Lat{Count} = IAGOS.lat;
  Flights.Paths.Lon{Count} = IAGOS.lon;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % store some geophysical information, if it exists
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if isfield(IAGOS,'zon_wind_AC'); Flights.Paths.U{Count} = IAGOS.zon_wind_AC; else Flights.Paths.U{Count} = NaN; end
  if isfield(IAGOS,'mer_wind_AC'); Flights.Paths.V{Count} = IAGOS.mer_wind_AC; else Flights.Paths.V{Count} = NaN; end
  if isfield(IAGOS,'air_temp_AC'); Flights.Paths.T{Count} = IAGOS.air_temp_AC; else Flights.Paths.T{Count} = NaN; end




  clear IAGOS

catch; end
end
textprogressbar(100)  ; textprogressbar('!')

clear Arr Dep Date Plane Instrument Data FlightFiles iFlight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% drop empty rows of data 
%caused by diff between valid flights and all flights in array declaration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Flights.Dep     = Flights.Dep(    1:Count);
Flights.Arr     = Flights.Arr(    1:Count);
Flights.t       = Flights.t(      1:Count);
Flights.Date    = Flights.Date(   1:Count);
Flights.PlaneID = Flights.PlaneID(1:Count);
Flights.InstID  = Flights.InstID( 1:Count);
clear Count



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%store the airport list used as well as the choices made about data handling, to check later
Settings.Airports = Airports;

%and save
save([Paths.StoreDir,'/iagos_flight_data_',Paths.SourceIdentifier,'.mat'],'Flights','Settings')

%and return
return
