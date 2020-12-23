clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analyse IAGOS data for effect of ENSO, NAO, QBO and HadCRUT on
%trans-atlantic flight times
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/11/22
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
Results.Dep  = NaN(numel(Flights),1); %index in airports array
Results.Arr  = Results.Dep;           %index in airports array
Results.t    = Results.Dep;           %flight time
Results.Date = Results.Dep;           %date of flight 

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
  %hence, identify start and end location, and date
  NFields = numel(MetaData);
  for iField=1:1:NFields
    if strcmp(MetaData(iField).Name,'departure_airport')
      Dep = MetaData(iField).Value; Dep = Dep(1:3);
    elseif strcmp(MetaData(iField).Name,'arrival_airport')
      Arr = MetaData(iField).Value; Arr = Arr(1:3);
    elseif strcmp(MetaData(iField).Name,'departure_UTC_time')
      Date = MetaData(iField).Value;
    end
  end
  clear MetaData NFields iField
  
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
  
  %store the start and end locations, and the flight time
  Results.Dep( iFlight) = find(contains(Airports,Dep));
  Results.Arr( iFlight) = find(contains(Airports,Arr));
  Results.t(   iFlight) = nanmax(Data.UTC_time) - nanmin(Data.UTC_time); %seconds
  Results.Date(iFlight) = Date;
catch;end
  
end
textprogressbar(100)  ; textprogressbar('!')
  

save('flightpairs.mat','Results','Airports','Settings')
