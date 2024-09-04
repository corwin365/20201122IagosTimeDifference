function process_emissions(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%transfer information to and from MEC emissions calculator
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/09/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++++++')
disp('Emissions calculator interface')
disp('+++++++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each flight, compute flight levels and distances of each cruise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iFlight=1:1:size(FlightData,1);

  %load data
  Flight = rCDF(FlightData.FilePath(iFlight));

  %convert altitudes to feet, then round off to the nearest 10 000
  %this gives us flight levels
  Altitude = round((Flight.baro_alt_AC.*3.28./100),-2);

  %split the data into cruises at a given flight level
  CruiseStarts = [1;find(diff(Altitude) ~=0)];

  %hence, find the length and height of each cruise
  Cruises = NaN(numel(idx)-1,2);
  for iCruise=1:1:numel(idx)-1;
    
    Cruises(iCruise,1) = Altitude(CruiseStarts(iCruise)); %flight level of each cruise

    startlat = Flight.lat(CruiseStarts(iCruise));
    startlon = Flight.lon(CruiseStarts(iCruise));    
    endlat   = Flight.lat(CruiseStarts(iCruise+1)-1);
    endlon   = Flight.lon(CruiseStarts(iCruise+1)-1);      
    Cruises(iCruise,2) = nph_haversine([startlat,startlon],[endlat,endlon])./0.54; %cruise length in nautical miles
    clear startlat startlon endlat endlon

  end





end