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
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get valid airport list
Airports = load([Settings.Paths.DataDir,'/',Settings.ID,'_airportinfo.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find files for the time period selected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flight filenames
FlightFiles = wildcardsearch(Settings.Paths.AeolusData,'*.nc*');


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
Reasons = {'list','region','dt','dLat','dLon'};
Rejected = zeros(numel(Reasons),1);
for iFlight=1:1:numel(FlightFiles)

  %get flight information
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


  %load the flight
  Data = rCDF(FlightFiles{iFlight});

  %get the file metadata
  MetaData = Data.MetaData.Attributes.Global;
  NFields = numel(MetaData);  

  %are we using the old or the new format?
  if strcmp(class(MetaData(1).Value),'string'); Format = 2; else Format = 1; end %newer format (post roughly 2022?) is 2, original flavour is 1

  %if we're using the newer data format, we need to convert the string arrays to char arrays
  if Format == 2;  for iField=1:1:NFields; MetaData(iField).Value = convertStringsToChars(MetaData(iField).Value ); end; end
  
  %hence, identify start and end location, date, and plane/instrument ID
  for iField=1:1:NFields
    if strcmp(MetaData(iField).Name,'departure_airport')
      Dep = MetaData(iField).Value(1:3);
    elseif strcmp(MetaData(iField).Name,'arrival_airport')
      Arr = MetaData(iField).Value(1:3);
    elseif strcmp(MetaData(iField).Name,'departure_UTC_time')
      Date = MetaData(iField).Value;
    elseif Format == 1 && strcmp(MetaData(iField).Name,'title');
      ID = MetaData(iField).Value;
    elseif Format == 2 && strcmp(MetaData(iField).Name,'platform');
      ID = MetaData(iField).Value;
    end
  end
  clear MetaData NFields iField

  %convert date to Matlab time
  Date = datenum(Date,'YYYY-mm-DDTHH:MM:ss');  

  %convert the ID string to an *aircraft* ID and an *instrument* ID
  ID = split(ID,',');
  Plane = ID{4}; Instrument = ID{1};
  Plane(Plane == ' ') = []; Instrument(Instrument == ' ') = []; %remove excess spaces
  clear ID


  % is this a valid flight?
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  %is each airport on one of our lists?
  if ~ismember(Arr,Airports.AirportData.Code); Rejected(1) = Rejected(1)+1;  continue; end
  if ~ismember(Dep,Airports.AirportData.Code); Rejected(1) = Rejected(1)+1; continue; end    

  %find which continent the arrival and depature are in
  Start.NA  = ismember(Dep,Airports.AirportData.Code(Airports.AirportData.Continent == 1));
  Start.Eur = ismember(Dep,Airports.AirportData.Code(Airports.AirportData.Continent == 2));
  End.NA    = ismember(Arr,Airports.AirportData.Code(Airports.AirportData.Continent == 1));
  End.Eur   = ismember(Arr,Airports.AirportData.Code(Airports.AirportData.Continent == 2));
    
  %are the two ends on different continents?
  if Start.NA  + End.NA  == 2; Rejected(2) = Rejected(2)+1; continue; end
  if Start.Eur + End.Eur == 2; Rejected(2) = Rejected(2)+1; continue; end
  clear Start End

  %qc: check maximum gap in data
  if Settings.Choices.Maxdt   < max(diff(Data.UTC_time)); Rejected(3) = Rejected(3)+1; continue; end
  if Settings.Choices.MaxdLat < max(diff(Data.lon));      Rejected(4) = Rejected(4)+1; continue; end
  if Settings.Choices.MaxdLon < max(diff(Data.lat));      Rejected(5) = Rejected(5)+1; continue; end

  

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

disp('Selection/QC rejections:')
for iReason=1:1:numel(Reasons)
  disp([Reasons{iReason},': ',num2str(Rejected(iReason))])
end; 
clear iReason Reasons Rejected


disp('--------------------------')
disp('Data from IAGOS formatted')
disp('--------------------------')



end

