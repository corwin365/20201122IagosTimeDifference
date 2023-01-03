function [] = bc_format_extra_data(Paths,Airports,Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find all flights between the chosen airports in Europe and North America
%then retain their metadata and some basic info (e.g. path, u, v, T)
%
%this loads data supplied by Ed Gryspeerdt in CSV format, rather than the
%highly structured IAGOS data. This is harder to clean, but much more 
%voluminous, which helps with the stats...
%
%Corwin Wright, c.wright@bath.ac.uk, 2022/07/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load IATA<->ICAO airport conversion
AirportConversion = [Paths.EdData,'/airports.csv'];

%flight filenames
% FlightFiles = wildcardsearch(Paths.EdData,'*.txt');
FlightFiles = wildcardsearch(Paths.EdData,'*tfms*.csv');

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

%get conversion table between airport identifiers
AirportConversion = readtable(AirportConversion);
AirportConversion = AirportConversion(:,[13,14,3,9]); %used only when finding unique airports, below

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find all new airports in the dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % % % % 
% % % % % % % % 
% % % % % % % % Store = {''};
% % % % % % % % for iFile=1:1:numel(FlightFiles)
% % % % % % % %   [iFile,numel(FlightFiles)]
% % % % % % % % 
% % % % % % % %   %get data
% % % % % % % %   Data = import_ed_flights_v2(FlightFiles{iFile});
% % % % % % % % 
% % % % % % % %   %trim down the flights so they are unique
% % % % % % % %   [~,idx] = unique(Data.UNIQUE_FLIGHT_ID);
% % % % % % % %   Data = Data(idx,:);
% % % % % % % %   clear idx
% % % % % % % % 
% % % % % % % %   %strip out bad characters introduced by (I assume) Ed's upstream data pipeline
% % % % % % % %   for iFlight=1:1:size(Data,1)
% % % % % % % %     a = strrep(Data.AIRCRAFT_ID(iFlight),"b'","");
% % % % % % % %     Data.AIRCRAFT_ID(iFlight) = strrep(a,"'","");
% % % % % % % % 
% % % % % % % %     a = strrep(Data.DEPT_APRT(iFlight),"b'","");
% % % % % % % %     Data.DEPT_APRT(iFlight) = strrep(a,"'","");
% % % % % % % % 
% % % % % % % %     a = strrep(Data.ARR_APRT(iFlight),"b'","");
% % % % % % % %     Data.ARR_APRT(iFlight) = strrep(a,"'","");
% % % % % % % %   end
% % % % % % % %   clear a iFlight
% % % % % % % % 
% % % % % % % %   %store airports
% % % % % % % %   Store = [Store;table2cell(Data(:,4))];
% % % % % % % %   Store = [Store;table2cell(Data(:,7))]; 
% % % % % % % %   Store = unique(string(Store(2:end)));
% % % % % % % % end
% % % % % % % % 
% % % % % % % % 
% % % % % % % % %convert airports codes to IATA
% % % % % % % % [DepList,~,idxD] = unique(Store);
% % % % % % % % for iDep=1:1:numel(DepList)
% % % % % % % %   [~,idx] = ismember(DepList(iDep),AirportConversion{:,1});
% % % % % % % %   if idx ~= 0; DepList(iDep) = AirportConversion{idx,2}; end
% % % % % % % % end
% % % % % % % % Store = DepList(idxD);
% % % % % % % % 
% % % % % % % % %discard all small airports, those not in relevant countries, and those we already have
% % % % % % % % for iAP=1:1:numel(Store)
% % % % % % % % 
% % % % % % % %   %size
% % % % % % % %   try
% % % % % % % %     [~,idx] = ismember(Store(iAP),AirportConversion{:,2});
% % % % % % % %     if ~strcmp(AirportConversion{idx,3},'large_airport');
% % % % % % % %       Store(iAP) = '';
% % % % % % % %     end
% % % % % % % %   catch; Store(iAP) = ''; %if it's not on the list it's probably not big...
% % % % % % % %   end
% % % % % % % % 
% % % % % % % %   %country - first check if it passed the above
% % % % % % % %   if strlength(Store(iAP)) ~= 0
% % % % % % % %     Country = AirportConversion{idx,4};
% % % % % % % %     ValidCountries = {'AT','BE','CH','CZ','DE','DK','ES','GB','IE','IT','LU','MC','NL','PT','US','CA'};
% % % % % % % %     if sum(strcmp(Country,ValidCountries)) == 0;
% % % % % % % %       Store(iAP) = '';
% % % % % % % %     end
% % % % % % % %   end
% % % % % % % % 
% % % % % % % %   %is it new?
% % % % % % % %   if strlength(Store(iAP)) ~= 0
% % % % % % % %     if sum(strcmp(Store(iAP),Airports.NA)) + sum(strcmp(Store(iAP),Airports.Eur)) ~= 0;
% % % % % % % %       Store(iAP) = '';
% % % % % % % %     end
% % % % % % % %   end
% % % % % % % %   
% % % % % % % % end
% % % % % % % % Store = unique(Store); Store = Store(2:end);
% % % % % % % % stop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, find all the flights between these east-west airport-pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%drop all non-major airports from the list. Searching this list is a dominant runtime sink otherwise.
Good = find(strcmp(AirportConversion{:,3},'large_airport') == 1);
AirportConversion = AirportConversion(Good,:);
clear Good

textprogressbar('Processing secondary flight files ')
Count = 0; %count of flights which meet criteria, for indexing purposes
for iFile = 1:1:numel(FlightFiles)
%   try
    textprogressbar(iFile./numel(FlightFiles).*100);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load data for this month
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %load file
    Data = import_ed_flights_v2(FlightFiles{iFile});

    %trim down the flights so they are unique 
    [~,idx] = unique(Data.UNIQUE_FLIGHT_ID);
    Data = Data(idx,:);
    clear idx

    %strip out bad characters introduced by (I assume) Ed's upstream data pipeline
    for iFlight=1:1:size(Data,1)
      a = strrep(Data.AIRCRAFT_ID(iFlight),"b'","");
      Data.AIRCRAFT_ID(iFlight) = strrep(a,"'","");

      a = strrep(Data.DEPT_APRT(iFlight),"b'","");
      Data.DEPT_APRT(iFlight) = strrep(a,"'","");

      a = strrep(Data.ARR_APRT(iFlight),"b'","");
      Data.ARR_APRT(iFlight) = strrep(a,"'","");      
    end
    clear a iFlight


    %convert airports codes to IATA
    [DepList,~,idxD] = unique(string(table2cell(Data(:,4))));
    [ArrList,~,idxA] = unique(string(table2cell(Data(:,7))));
    for iDep=1:1:numel(DepList)
      [~,idx] = ismember(DepList(iDep),AirportConversion{:,1});
      if idx ~= 0; DepList(iDep) = AirportConversion{idx,2}; end
    end
    for iArr=1:1:numel(ArrList)
      [~,idx] = ismember(ArrList(iArr),AirportConversion{:,1});
      if idx ~= 0; ArrList(iArr) = AirportConversion{idx,2}; end
    end
    Data{:,4} = DepList(idxD);
    Data{:,7} = ArrList(idxA);
    clear iDep iArr idx DepList idxD ArrList idxA 

    %remove airports where AT LEAST ONE of the terminals is not in our list
    idxStore = zeros(size(Data,1),1);
    for iAirport=1:1:numel(Airports.IDs)
      Dep = find(strcmp(Airports.IDs{iAirport},Data{:,4}));
      Arr = find(strcmp(Airports.IDs{iAirport},Data{:,7}));
      if numel(Dep) +numel(Arr) == 0;  
        continue;
    
      end
      idxStore(Dep) = idxStore(Dep) + 1;
      idxStore(Arr) = idxStore(Arr) + 1;
    end
    Data = Data(idxStore == 2,:);
    clear idxStore iAirport Dep Arr


    %ok, now the data are manageable! select the array and start formatting them to match the IAGOS data
    Dep = Data{:,4};
    Arr = Data{:,7};
    Plane  = table2cell(Data(:,1)); %not the same, but we don't have the aircraft tail IDs
    Instrument = cell(numel(Plane),1);
    for iCell=1:1:numel(Instrument); Instrument{iCell} = 'NOT_IAGOS'; end%use this to filter out the above when tracking individual planes
    clear iCell

    %merge timestamp parts
    DepTime= datenum(string(table2array(Data(:,5))),'yyyymmdd') + table2array(Data(:,6))./24;
    ArrTime= datenum(string(table2array(Data(:,8))),'yyyymmdd') + table2array(Data(:,9))./24;

    %check the times agree with those Ed computed, to within 30 minutes. Sometimes they don't.
    %this is usually one of the dates being bad (19000101 is a common one)
    %for now just delete - on first few tests seems only a small effect. if this hits the data hard, come back and upgrade.
    for iFlight=1:1:numel(Dep)
      if abs((ArrTime(iFlight)-DepTime(iFlight)) - Data.FLIGHT_TIME(iFlight)./24) > 1./48;
        ArrTime(iFlight) = NaN;
        DepTime(iFlight) = NaN;
      end
    end; clear iFlight

    t = ArrTime-DepTime; t = t.*24.*60.*60; %seconds
    Date = floor(DepTime);

    %store
    for iFlight=1:1:numel(t)
      %store the information we want
      Count = Count+1;
      Flights.Dep{    Count} = Dep{iFlight};
      Flights.Arr{    Count} = Arr{iFlight};
      Flights.t(      Count) = t(iFlight);
      Flights.Date(   Count) = Date(iFlight);
      Flights.PlaneID{Count} = Plane{iFlight};
      Flights.InstID{ Count} = Instrument{iFlight};
    end

    %tidy up
    clearvars -except Flights FlightFiles Count Settings Paths Airports  AirportConversion

%   catch; end
end
textprogressbar(100)  ; textprogressbar('!')
clear AirportConversion FlightFiles


% % % %transpose, to match other dataset
% % % Fields= fieldnames(Flights);
% % % for iF=1:1:numel(Fields); Flights.(Fields{iF}) = Flights.(Fields{iF})'; end
% % % clear iF Fields

%create empty flight traces, to protect any logic using them from IAGOS
Flights.Paths = struct();
Fields = {'Lat','Lon','U','V','T'};
for iF=1:1:numel(Fields);
  Flights.Paths.(Fields{iF}) = cell(1,numel(Flights.Dep));
end
clear iF Fields



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%store the airport list used as well as the choices made about data handling, to check later
Settings.Airports = Airports;

%and save
save([Paths.StoreDir,'/extra_flight_data_',Paths.SourceIdentifier,'.mat'],'Flights','Settings')

%and return
return


