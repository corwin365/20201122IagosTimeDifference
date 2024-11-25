function process_scalefactors(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute estimates of CO2 emissions, USD fuel costs, and total flights
% per month
%
%also saves the regression results
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/11/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++++++')
disp('Scale Factor Computations')
disp('+++++++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flight data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])

%load fuel price index
Fuel = load([Settings.Paths.Indices,'jet_fuel_price.mat']);

%load aviation data
AllFlights = load([LocalDataDir,'/Miscellany/flight_stats_us.mat']);
AllFlights = AllFlights.USBTSconcatflights;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate fuel use per km for each plane in the dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plane type and properties for each tail ID in the dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%list of plane types, so we can refer to it later
PlaneTypeList = {'A340-300','A340-200','A330-300','A340-600'};


%assign a plane type to each flight (plane type for each tail numnber found from planespotters.net)
PlaneType = NaN(size(FlightData,1),1);
for iFlight=1:1:numel(PlaneType)

  if ismissing(FlightData.Plane(iFlight)); continue; end %round trips don't have individual planes assigned

  switch FlightData.Plane(iFlight)
    case 'D-AIGI'; PlaneType(iFlight) = 1; %A340-300
    case 'D-AIGF'; PlaneType(iFlight) = 1; %A340-300
    case 'V5-NME'; PlaneType(iFlight) = 1; %A340-300
    case 'D-AIGT'; PlaneType(iFlight) = 1; %A340-300
    case 'F-GLZG'; PlaneType(iFlight) = 1; %A340-300
    case 'F-GLZU'; PlaneType(iFlight) = 1; %A340-300
    case 'C-GEFA'; PlaneType(iFlight) = 1; %A340-300
    case 'EC-GUQ'; PlaneType(iFlight) = 1; %A340-300      

    case 'OE-LAG'; PlaneType(iFlight) = 2; %A340-200      

    case 'D-AIKO'; PlaneType(iFlight) = 3; %A330-300
    case 'D-AIKE'; PlaneType(iFlight) = 3; %A330-300
    case 'F-GZCO'; PlaneType(iFlight) = 3; %A330-300

     
    case 'D-AIHE'; PlaneType(iFlight) = 4; %A340-600
      
    otherwise;  %this shouldn't happen
      disp('Warning: plane type not included in fuel cost calculator')
      stop
      PlaneType(iFlight) = 0;

  end
end

%a litre of kerosene weighs:
litreweight = 0.819; %kg

% % % %the fuel calculation uses the spreadsheet obtained here:
% % % % https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/2HMEHB
% % % 
% % % %these are the properties needed for this fuel use calculation, all in kg:
% % %  % (type 1) a340-300: https://www.airbus.com/en/who-we-are/company-history/commercial-aircraft-history/previous-generation-aircraft/a340-family/a340-300
% % %  % (type 2)  a340-200: https://www.airbus.com/en/who-we-are/company-history/commercial-aircraft-history/previous-generation-aircraft/a340-family/a340-200
% % %  % (type 3)  a330-300: https://aircraft.airbus.com/en/aircraft/a330-advanced-to-boost-profitability/a330-300
% % %  % (type 4)  a340-600: https://www.airbus.com/en/who-we-are/company-history/commercial-aircraft-history/previous-generation-aircraft/a340-family/a340-600
% % % 
% % % %maximum take-off weight:
% % % MTOW(1) = 276.50.*1000;
% % % MTOW(2) = 275.*1000;
% % % MTOW(3) = 242.*1000; %kg
% % % MTOW(4) = 380.*1000;
% % % 
% % % %maximum zero fuel weight
% % % MZFW(1) = 183.*1000;
% % % MZFW(2) = 169.*1000;
% % % MZFW(3) = 175.*1000; %lb->kg
% % % MZFW(4) = 251.*1000;
% % % 
% % % %maximum fuel weight
% % % MFW(1) = 147850 .* litreweight;
% % % MFW(2) = 155040 .* litreweight;
% % % MFW(3) = 139090 .*litreweight;
% % % MFW(4) = 204500 .* litreweight;
% % % 
% % % 
% % % %other assumptions: 
% % % %passenger weight: 95km (inc bags!)
% % % %flight distance: 7500 km (fra to atl)
% % % %distance to alternate airport: 370km (standard value)
% % % %cruise mach number: 0.72

%typical number of passengers
Passengers(1) = 270; %1: A340-300
Passengers(2) = 230; %2: A340-200
Passengers(3) = 300; %3: A330-300
Passengers(4) = 350; %5: A340-600

%taking a mid-cruise value, putting these values into the spreadsheets gives fuel estimates in kg per 100 passenger-km of:
FuelUse(1) = 3.03;
FuelUse(2) = 4.07;
FuelUse(3) = 2.00;
FuelUse(4) = 3.27;

%which converts to the following in gallons per km
FuelUse = FuelUse ./ litreweight ./ 3.785 ./100 .* Passengers ; %3.785 is the conversion from litres to US gallons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert to cost per minute of delay per plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tracks = load([Settings.Paths.DataDir,'/',Settings.ID,'_flighttracks.mat']);
Tracks = Tracks.Store;

UniquePlanes = unique(FlightData.Plane(~ismissing(FlightData.Plane)));
SpeedStore = NaN(numel(UniquePlanes),numel(Tracks.Lat));

for iFlight=1:1:numel(Tracks.Lat)

  %get plane index
  PlaneIdx = find(UniquePlanes == FlightData.Plane(iFlight));

  %take middle third of flight and compute time and distance
  Lon  = Tracks.Lon{ iFlight};
  Lat  = Tracks.Lat{ iFlight};
  Time = Tracks.Time{iFlight};
  a = round(numel(Lat)./3); b = 2.*a; %i.e. from the point 1/3 of the way along to the point2/3 along
  if a == 0; continue; end
  dx = nph_haversine([Lat(a),Lon(a)],[Lat(b),Lon(b)]).*1000;
  dt = (Time(b)-Time(a)).*60.*60.*24;
  SpeedStore(PlaneIdx ,iFlight) = dx./dt;
end
clear iFlight PlaneIdx Lon Lat Time a b dx dt
SpeedPerPlane = nanmean(SpeedStore,2)';
%these are all pretty close (standard deviation < 3m/s), so just take the mean to simply everything else
PlaneSpeed = nanmean(SpeedPerPlane(:)); %m/s

%hence, how much fuel do we use per minute, in gallons?
GallonsPerMinute = FuelUse .* (PlaneSpeed./1000.*60);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate CO2 emissions per minute delay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CO2 emission factor for Jet and Turboprop aircraft 	
Co2PerKgFuel = 3.15; %from Eurocontrol

%convert from per kg of fuel to per gallon of fuel
Co2PerGallon = Co2PerKgFuel .* litreweight .*3.785; %co2 per gallon of fuel. this is consistent with the US EIA estimate of 9.75 kg-co2 per gallon of jet fuel.

%hence, Co2 per minute of delay
Co2PerMinute = GallonsPerMinute .* Co2PerGallon;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate total number of flights per day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%create time axis
AllFlights.Time = datenum(AllFlights.year,AllFlights.month,1);
TimeScale = unique(AllFlights.Time(~isnan(AllFlights.Time)));

%how many seats per plane on this route?
AllFlights.SeatsPerPlane = AllFlights.seats ./ AllFlights.departures_performed;
AllFlights.SeatsPerPlane(~isfinite(AllFlights.SeatsPerPlane)) = NaN;

%deduplicate the table 
Tbl = AllFlights(:,[1,2,3,7]); %year, month, origin, dest
[~,indexToUniqueRows,~] = unique(Tbl);
AllFlights = AllFlights(indexToUniqueRows,:);
clear indexToUniqueRows

%filtering
%%%%%%%%%%%%%%

% Good = find(Data.Time < datenum(2022,1,1));
% Data = Data(Good,:);


%require a minimum travel distance of 4000km
Good = find(AllFlights.distance > 4000./1.6);
AllFlights = AllFlights(Good,:);

%require jet aircraft (aircraft_group >=6);
Good = find(AllFlights.Aircraft_Group >= 6);
AllFlights = AllFlights(Good,:);

%require passenger aircraft (aircraft_config = 1)
Good = find(AllFlights.Aircraft_Config == 1);
AllFlights = AllFlights(Good,:);

%require at least 150 seats per plane
Good = find(AllFlights.SeatsPerPlane >= 150);
AllFlights = AllFlights(Good,:);

clear Good

%find flights between North America and Europe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%us locations have WACs less than 100. 
%    Exclude all outside contiguous US: these are WAC < 10
%canadian locations have WACS in the 900 range
%european locations have WACs in the 400 range

Origin_NA      = find((AllFlights.origin_wac >  10 & AllFlights.origin_wac <  100) | (AllFlights.origin_wac >= 900 & AllFlights.origin_wac < 1000));
Origin_Eur     = find(AllFlights.origin_wac >= 400 & AllFlights.origin_wac < 500);

Dest_NA      = find((AllFlights.dest_wac >  10 & AllFlights.dest_wac <  100) | (AllFlights.dest_wac >= 900 & AllFlights.dest_wac < 1000));
Dest_Eur     = find(AllFlights.dest_wac >= 400 & AllFlights.dest_wac < 500);

Data_Westwards = AllFlights(intersect(Dest_NA, Origin_Eur),:);
Data_Eastwards = AllFlights(intersect(Dest_Eur, Origin_NA),:);
clear Origin_NA Origin_Eur Dest_NA Dest_Eur AllFlights

%hence, estimate total flights per day in each direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get monthly total
NFlights = NaN(2,numel(TimeScale)); %1 is westwards, 2 is eastwards
for iMonth=1:1:numel(TimeScale)

  idx = find(Data_Westwards.Time == TimeScale(iMonth));
  NFlights(1,iMonth) = nansum(Data_Westwards.departures_performed(idx));

  idx = find(Data_Eastwards.Time == TimeScale(iMonth));
  NFlights(2,iMonth) = nansum(Data_Eastwards.departures_performed(idx));
end

%scale to days and put on our primary daily data scale
NFlightsb = NFlights.*12./365; clear NFlights

NFlights.Time =Settings.Choices.TimeRange(1):1:Settings.Choices.TimeRange(2);
NFlights.Val = interp1(TimeScale,NFlightsb',NFlights.Time,'nearest');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([Settings.Paths.DataDir,'/',Settings.ID,'_fuel_and_emissions.mat'],'GallonsPerMinute','Co2PerMinute','Fuel','PlaneType','PlaneTypeList','NFlights')
