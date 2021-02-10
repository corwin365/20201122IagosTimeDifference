function [] = bb_prep_data(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find all flights between the chosen airports in Europe and North America
%
%additional data checks:
%A. remove any data within 100km of start and end of flight, to avoid 
%looping around the airport before or after departure
%B. modified to retain unique identifiers for individual flights, to 
%identify any issues due to equipment change
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flight filenames
Flights = wildcardsearch(Settings.DataDir,'*.nc');

%possible airports, and their locations
Airports.IDs  = [Settings.Eur,Settings.NA];
Airports.Lons = NaN(numel(Airports.IDs),1);
Airports.Lats = Airports.Lons;

%metainfo about flights
Results.Dep     = cell(numel(Flights),1);  %flight start
Results.Arr     = cell(numel(Flights),1);  %flight end
Results.t       = NaN( numel(Flights),1);   %flight time
Results.Date    = NaN( numel(Flights),1);   %date of flight 
Results.PlaneID = cell(numel(Flights),1);  %unique aircraft identifier
Results.InstID  = cell(numel(Flights),1);  %unique instrument identifier

%downsampled flight tracks, for map generation
Results.Paths   = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% where are the airports?
%this information isn't used here, but will be needed to generate maps later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Codes = load('data/airport_codes.mat');
for iAirport=1:1:numel(Airports.IDs)
  
  ID = Airports.IDs{iAirport};
  ID = find(contains(Codes.Codes,ID));
  Coords = strsplit(Codes.Coords{ID},' ');
  
  Airports.Lons(iAirport) = str2num(Coords{1});
  Airports.Lats(iAirport) = str2num(Coords{2});  
end
clear Codes iAirport ID Coords

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, find all the flights between these east-west airport-pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textprogressbar('Processing flights ')
Count = 0; %count of flights which meet criteria, for indexing purposes
for iFlight = 1:1:numel(Flights)
try
 if mod(iFlight,100) == 0; textprogressbar(iFlight./numel(Flights).*100); end
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % process metadata
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %load the flight 
  Data = rCDF(Flights{iFlight});
  
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
  Date =datenum(Date,'YYYY-mm-DDTHH:MM:ss');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % is this a valid flight?
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %is each airport on one of our lists?
  if ~ismember(Arr,Airports.IDs); continue; end
  if ~ismember(Dep,Airports.IDs); continue; end    
  
  %find which continent the arrival and depature are in
  Start.NA  = ismember(Dep,Settings.NA);
  Start.Eur = ismember(Dep,Settings.Eur);
  End.NA    = ismember(Arr,Settings.NA);
  End.Eur   = ismember(Arr,Settings.Eur);    
    
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
  Results.Dep{    Count} = Dep;
  Results.Arr{    Count} = Arr;
  Results.t(      Count) = nanmax(Data.UTC_time) - nanmin(Data.UTC_time); %seconds
  Results.Date(   Count) = Date;
  Results.PlaneID{Count} = Plane;
  Results.InstID{ Count} = Instrument;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % store flight path
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  
  %also store the full flight path, downsampled to ~Settings.MinTime
  IAGOS = prep_iagos(Flights{iFlight},'SamplingRate',Settings.ResampleTime./24./60);
  Results.Paths.Lat{Count} = IAGOS.lat;
  Results.Paths.Lon{Count} = IAGOS.lon;
  
  clear IAGOS



catch; end
end
textprogressbar(100)  ; textprogressbar('!')

clear Arr Dep Date Plane Instrument Data Flights iFlight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% drop empty rows of data 
%caused by diff between valid flights and all flights in array declaration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Results.Dep     = Results.Dep(    1:Count);
Results.Arr     = Results.Arr(    1:Count);
Results.t       = Results.t(      1:Count);
Results.Date    = Results.Date(   1:Count);
Results.PlaneID = Results.PlaneID(1:Count);
Results.InstID  = Results.InstID( 1:Count);
clear Count



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('data/flight_data.mat','Results','Airports','Settings')
