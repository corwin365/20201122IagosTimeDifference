function import_iagos_data(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%standardise IAGOS data
%
%finds all flights between the chosen airports in Europe and North America
%then retain their metadata, including the name of the original file
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Formatting data from IAGOS')
disp('+++++++++++++++++++++++++++')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find files for the time period selected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flight filenames
FlightFiles = wildcardsearch(Settings.Paths.AeolusData,'*.nc');

%find a date for each file
idx1 = strfind(FlightFiles,'IAGOS_timeseries_');
Dates = NaN(numel(idx1),1);
for iFile=1:1:numel(FlightFiles)
  text = FlightFiles{iFile};
  datestring = text(idx1{iFile}+[1:1:8]+16);
  Dates(iFile) = datenum(datestring,'yyyymmdd');
end
clear iFile idx1 text datestring

%hence, select only those date we want
InTimeRange = find(Dates >= min(Settings.Choices.TimeRange) ...
                 & Dates <= max(Settings.Choices.TimeRange));

FlightFiles = FlightFiles(InTimeRange);
clear Dates InTimeRange

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the information we need from each file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a data table to keep the values in
VarNames = [   "Dep",   "Arr",     "t",  "Date", "Plane","Instrument","DataSource","FilePath"];
VarTypes = ["string","string","double","double","string",    "string",    "string",  "string"];
FlightData = table('Size',[numel(FlightFiles),numel(VarNames)], ...
                   'VariableTypes',VarTypes,                    ...
                   'VariableNames',VarNames                     );
clear VarNames VarTypes
Filled = zeros([numel(FlightFiles),1]);


textprogressbar('Importing IAGOS flights ')
for iFlight=1:1:numel(FlightFiles)

  %get flight information
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


  %load the flight
  Data = rCDF(FlightFiles{iFlight});

  %get the file metadata
  MetaData = Data.MetaData.Attributes.Global;  
  
  %hence, identify start and end location, date, and plane/instrument ID
  NFields = numel(MetaData);
  for iField=1:1:NFields
    if strcmp(MetaData(iField).Name,'departure_airport')
      Dep = MetaData(iField).Value(1:3);
    elseif strcmp(MetaData(iField).Name,'arrival_airport')
      Arr = MetaData(iField).Value(1:3);
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


  % is this a valid flight?
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  %is each airport on one of our lists?
  if ~ismember(Arr,[Settings.Airports.NA,Settings.Airports.Eur]); continue; end
  if ~ismember(Dep,[Settings.Airports.NA,Settings.Airports.Eur]); continue; end    

  %find which continent the arrival and depature are in
  Start.NA  = ismember(Dep,Settings.Airports.NA);
  Start.Eur = ismember(Dep,Settings.Airports.Eur);
  End.NA    = ismember(Arr,Settings.Airports.NA);
  End.Eur   = ismember(Arr,Settings.Airports.Eur);    
    
  %are the two ends on different continents?
  if Start.NA  + End.NA  == 2; continue; end
  if Start.Eur + End.Eur == 2; continue; end
  clear Start End
  
  %remove any points within Settings.MinDist of DEP or ARR
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  distanceFromDep = nph_haversine([ones(numel(Data.lat),1).*Data.lat(1),ones(numel(Data.lat),1).*Data.lon(1)], ...
                                  [Data.lat,Data.lon]);
  distanceFromArr = nph_haversine([ones(numel(Data.lat),1).*Data.lat(end),ones(numel(Data.lat),1).*Data.lon(end)], ...
                                  [Data.lat,Data.lon]);
  Bad = find(distanceFromArr < Settings.Choices.MinDist | distanceFromDep < Settings.Choices.MinDist);
  
  Data.UTC_time(Bad) = NaN;
  clear Bad distanceFromArr distanceFromDep

  % store!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
  Filled(iFlight) = 1;
  
  FlightData{iFlight,1} = string(Dep);
  FlightData{iFlight,2} = string(Arr);
  FlightData{iFlight,3} = nanmax(Data.UTC_time) - nanmin(Data.UTC_time); %seconds
  FlightData{iFlight,4} = Date;
  FlightData{iFlight,5} = string(Plane);
  FlightData{iFlight,6} = string(Instrument);
  FlightData{iFlight,7} = "IAGOS";
  FlightData{iFlight,8} = string(FlightFiles{iFlight});


  if mod(iFlight,100) == 0;textprogressbar(iFlight./numel(FlightFiles).*100); end
end; clear iFile
textprogressbar(100);textprogressbar('!')

FlightData = FlightData(Filled == 1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_iagos.mat'],'FlightData')

disp('--------------------------')
disp('Data from IAGOS formatted')
disp('--------------------------')

end

