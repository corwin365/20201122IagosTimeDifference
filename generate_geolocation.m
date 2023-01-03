function [] = generate_geolocation(Settings)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%geolocate airport codes
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Generating geolocation data')
disp('+++++++++++++++++++++++++++')

%load big file of airport metadata
InfoArray = readtable([Settings.Paths.DataDir,'/airports.csv']);
InfoArray = InfoArray(:,[5,6,7,8,9,14,]); %lat,lon,elev,continent,country,iata code

%drop any that don't have an IATA code
[~,idx] = rmmissing(InfoArray.iata_code);
InfoArray(idx,:) = [];
clear idx

%create merged list of airports
Codes = [Settings.Airports.NA,Settings.Airports.Eur]';

%which continent is each one on?
%1 for NA, 2 for Eur
Continent = [ones(size(Settings.Airports.NA)).*1,ones(size(Settings.Airports.Eur)).*2]'; 

%create storage array for geolocation
Coords = NaN(numel(Codes),3); %lon,lat,z

%fill this cell array
IATA = InfoArray.iata_code;
for iAirport=1:1:numel(Codes)

  %find code in InfoArray
  idx = find(strcmp(Codes{iAirport},IATA));
  
  %pull out lat and lon as string (see above - backwards compatability)
  Coords(iAirport,:) = [InfoArray.longitude_deg(idx), ...
                        InfoArray.latitude_deg( idx), ...
                        InfoArray.elevation_ft( idx)];
end
clear IATA idx iAirport

%put into table
AirportData = table(Codes,Coords(:,1),Coords(:,2),Coords(:,3),Continent);
AirportData = renamevars(AirportData,["Codes","Var2","Var3","Var4"],["Code","Lon","Lat","Elev"]);

%save results
save([Settings.Paths.DataDir,'/',Settings.ID,'_airportinfo.mat'],'AirportData')

disp('--------------------------')
disp('Geolocation data generated')
disp('--------------------------')
