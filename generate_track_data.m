function [] = generate_track_data(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate flight track data used in subsequent analyses

%Corwin Wright, c.wright@bath.ac.uk, 2024/09/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Generating track data')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_merged.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%edit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get wind data based on this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate a list of dates that will include every flight in the dataset
FullDateRange = (floor(min(FlightData.Date))-1):1:(ceil(max(FlightData.Date))+1);

%create a store for the flight traces
FlightStore = {};
textprogressbar('Importing tracks ')
for iFlight=1:1:size(FlightData,1)

  if FlightData.DataSource(iFlight) ~= "IAGOS"; continue; end

  %I generated the flight list on my laptop but want to also run it on 0184,
  %so this line just switches filepaths if needed
  FilePath = FlightData.FilePath(iFlight);
  FilePath = strrep(FilePath,'D:\Data\',LocalDataDir);
  FilePath = strrep(FilePath,'\','/');

  %load data and store
  Data =  rCDF(FilePath);

  try; %a very small number of flights are missing some variables, this stops them crashing the routine
  Store.U{iFlight}    = Data.zon_wind_AC;
  Store.V{iFlight}    = Data.mer_wind_AC;
  Store.Lat{iFlight}  = Data.lat;
  Store.Lon{iFlight}  = Data.lon;
  Store.Time{iFlight} = FlightData.Date(iFlight)+Data.UTC_time/86400;
  Store.Prs{iFlight}  = Data.air_press_AC;
  Store.Z{iFlight}    = Data.baro_alt_AC;
  catch; continue; end

  textprogressbar(iFlight./size(FlightData,1).*100)

end
textprogressbar('!')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([Settings.Paths.DataDir,'/',Settings.ID,'_flighttracks.mat'],'Store','-v7.3')



disp('--------------------------')
disp('Track data generated')
disp('--------------------------')
