function import_ed_data(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%standardise Ed's data
%
%%finds all flights between the chosen airports in Europe and North America
%then retain their metadata, including the name of the original file
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Formatting data from Ed G')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load IATA<->ICAO airport conversion
AirportConversion = [Settings.Paths.DataDir,'/airports.csv'];

%flight filenames
FlightFiles = wildcardsearch(Settings.Paths.EdData,'*tfms*.csv');

%possible airports, and their locations
Airports.IDs  = [Settings.Airports.Eur,Settings.Airports.NA];
Airports.Lons = NaN(numel(Airports.IDs),1);
Airports.Lats = Airports.Lons;

%get conversion table between airport identifiers
AirportConversion = readtable(AirportConversion);
AirportConversion = AirportConversion(:,[13,14,3,9]); %used only when finding unique airports, below

%drop all non-major airports from the list. Searching this list is a dominant runtime sink otherwise.
Good = find(strcmp(AirportConversion{:,3},'large_airport') == 1);
AirportConversion = AirportConversion(Good,:);
clear Good

%find a date for each file
idx1 = strfind(FlightFiles,'tfms_');
Dates = NaN(numel(idx1),1);
for iFile=1:1:numel(FlightFiles)
  thetext = FlightFiles{iFile};
  datestring = thetext(idx1{iFile}+[1:1:7]+4);
  Dates(iFile) = datenum(str2num(datestring(1:4)),1,str2num(datestring(5:7)));
end
clear iFile idx1 text datestring

%hence, select only those date we want
InTimeRange = find(Dates >= min(Settings.Choices.TimeRange) ...
                 & Dates <= max(Settings.Choices.TimeRange));

FlightFiles = FlightFiles(InTimeRange);
clear Dates InTimeRange



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, find all the flights between the selected east-west airport-pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a data table to keep the values in
VarNames = [   "Dep",   "Arr",     "t",  "Date", "Plane","Instrument","DataSource","FilePath"];
VarTypes = ["string","string","double","double","string",    "string",    "string",  "string"];
FlightData = table('Size',[numel(FlightFiles).*100,numel(VarNames)], ... %placeholder length
                   'VariableTypes',VarTypes,                         ...
                   'VariableNames',VarNames                          );
clear VarNames VarTypes

textprogressbar('Processing flight files from Ed Gryspeerdt ')
Count = 0; %count of flights which meet criteria, for indexing purposes
for iFile = 1:1:numel(FlightFiles)
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
      Count = Count+1;
      FlightData{Count,1} = string(Dep{iFlight});
      FlightData{Count,2} = string(Arr{iFlight});
      FlightData{Count,3} = t(iFlight);
      FlightData{Count,4} = Date(iFlight);
      FlightData{Count,5} = string(Plane{iFlight});
      FlightData{Count,6} = string(Instrument{iFlight});
      FlightData{Count,7} = "FromEd";
      FlightData{Count,8} = string(FlightFiles{iFile});
    end

    %tidy up
    clearvars -except Flights FlightFiles Count Settings Paths Airports  AirportConversion FlightData

end
textprogressbar(100)  ; textprogressbar('!')
clear AirportConversion FlightFiles

%drop any empty rows
FlightData = rmmissing(FlightData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_edG.mat'],'FlightData')

disp('--------------------------')
disp('Data from Ed G formatted')
disp('--------------------------')

end
