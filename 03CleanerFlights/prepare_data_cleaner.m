clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyse IAGOS data for effect of ENSO, NAO, QBO and HadCRUT on
%trans-atlantic flight times
%
%
%additional data checks:
%1. remove any data within 100km of start and end of flight, to avoid looping around the airport before or after departure
%2. modified to retain unique identifiers for individual flights, to identify any issues due to equipment change
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time period
Settings.TimeScale = datenum(1994,1,1):1:datenum(2020,12,31);

%datadir
Settings.DataDir = [LocalDataDir,'/IAGOS/Timeseries/'];

%definition of trans-Atlantic flight:
%between one of these 
Settings.NA = {'ATL','BOS','BWI','CDW','CLE','CLT','CVG','DRM','DTW','EWR','FOK','IAD','JFK','LUK','MKE','MRB','ORD','PHL','PNE','YMX','YQB','YUL','YYR','YYZ','YZD'};
%and one of these
Settings.Eur = {'AGA','AGP','AHO','AMM','AMS','ATH','AYT','BCN','BEY','BOD','BRE','BRU','BTS','BUD','CAI','CDG','CGN','CIA','CRL','DBV','DLM','DME','DRS','DUS','ESB','FCO','FKB','FRA','GHF','GRO','HAJ','HAM','HEL','HER','HSK','IST','LCA','LEI','LEJ','LGW','LHR','LIS','LNZ','LYS','MAD','MAN','MLA','MRS','MUC','MXP','NCE','NUE','ORY','OST','OTP','PMI','PRG','PSA','PUY','RHO','RIX','RLG','SDV','SKG','SNN','SPM','STN','SXB','SZG','SZW','TLS','TLV','TOJ','TXL','UTC','VIE','ZNV','ZQL','ZRH'};
%in either direction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flight filenames
Flights = wildcardsearch(Settings.DataDir,'*.nc');

%possible airports
Airports = [Settings.Eur,Settings.NA];

%results
Results.Dep     = NaN(numel(Flights),1); %index in airports array
Results.Arr     = Results.Dep;           %index in airports array
Results.t       = Results.Dep;           %flight time
Results.Date    = Results.Dep;           %date of flight 
Results.PlaneID = cell(numel(Flights));  %unique aircraft identifier
Results.InstID  = Results.PlaneID;       %unique instrument identifier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textprogressbar('Processing flights ')
for iFlight = 1:1:numel(Flights)

 if mod(iFlight,100) == 0; textprogressbar(iFlight./numel(Flights).*100); end
 
 try
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % process metadata
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %load the flight 
  Data = rCDF(Flights{iFlight});
  
  %get the dep and arr airports
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
  clear ID
  
  %convert date to Matlab time
  Date =datenum(Date,'YYYY-mm-DDTHH:MM:ss');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % is this a valid flight?
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %is it worth checking for a valid pair?
  if ~ismember(Arr,Airports); continue; end
  if ~ismember(Dep,Airports); continue; end    
  
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
  
  %1. remove any points within 100km of DEP or ARR
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  CutOff = 100;
  distanceFromDep = nph_haversine([ones(numel(Data.lat),1).*Data.lat(1),ones(numel(Data.lat),1).*Data.lon(1)], ...
                                  [Data.lat,Data.lon]);
  distanceFromArr = nph_haversine([ones(numel(Data.lat),1).*Data.lat(end),ones(numel(Data.lat),1).*Data.lon(end)], ...
                                  [Data.lat,Data.lon]);
  Bad = find(distanceFromArr < CutOff | distanceFromDep < CutOff);
  
  Data.UTC_time(Bad) = NaN;
  clear Bad distanceFromArr distanceFromDep
  
  %2. ... TBD
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % store
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  
  %store the information we want
  Results.Dep(    iFlight) = find(contains(Airports,Dep));
  Results.Arr(    iFlight) = find(contains(Airports,Arr));
  Results.t(      iFlight) = nanmax(Data.UTC_time) - nanmin(Data.UTC_time); %seconds
  Results.Date(   iFlight) = Date;
  Results.PlaneID{iFlight} = Plane;
  Results.InstID{ iFlight} = Instrument;
 
  
  clear Arr Dep Date Plane Instrument
 catch;end
  
end
textprogressbar(100)  ; textprogressbar('!')
  

save('flightpairs_identified.mat','Results','Airports','Settings')
