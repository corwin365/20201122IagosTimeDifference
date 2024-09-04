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

%drop any that don't have an IATA code
[~,idx] = rmmissing(InfoArray.iata_code);
InfoArray(idx,:) = [];
clear idx

%drop anything except a 'large' or 'medium' airport
idx = find(contains(InfoArray.type,'large') | contains(InfoArray.type,'medium'));
InfoArray = InfoArray(idx,:);
clear idx


InfoArray = InfoArray(:,[5,6,7,8,9,14,]); %lat,lon,elev,continent,country,iata code



if strcmpi(Settings.Choices.Airports,'list');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% choose airports - whitelist approach
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %create merged list of airports
  Codes = [Settings.Airports.List.NA,Settings.Airports.List.Eur]';

  %which continent is each one on?
  %1 for NA, 2 for Eur
  Continent = [ones(size(Settings.Airports.List.NA)).*1,ones(size(Settings.Airports.List.Eur)).*2]';


elseif strcmpi(Settings.Choices.Airports,'geo');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% choose airports -  geographic approach
  %impose a minimum airport size, as otherwise we end up with thousands
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %North America
  idx.NA  = find(InfoArray.longitude_deg >= Settings.Airports.Bounds.NA(1) ...
               & InfoArray.longitude_deg <= Settings.Airports.Bounds.NA(2) ...
               & InfoArray.latitude_deg  >= Settings.Airports.Bounds.NA(3) ...
               & InfoArray.latitude_deg  <= Settings.Airports.Bounds.NA(4));
  %Europe
  idx.Eur = find(InfoArray.longitude_deg >= Settings.Airports.Bounds.Eur(1) ...
               & InfoArray.longitude_deg <= Settings.Airports.Bounds.Eur(2) ...
               & InfoArray.latitude_deg  >= Settings.Airports.Bounds.Eur(3) ...
               & InfoArray.latitude_deg  <= Settings.Airports.Bounds.Eur(4));

  Codes = InfoArray.iata_code([idx.NA;idx.Eur]);
  Continent = [ones(size(idx.NA)).*1;ones(size(idx.Eur)).*2];

else
  disp('Error - incorrect approach specified for identifying airports')
  stop
end


%create storage array for geolocation
Coords = NaN(numel(Codes),3); %lon,lat,z

%fill this cell array
IATA = InfoArray.iata_code;
for iAirport=1:1:numel(Codes)

  %find code in InfoArray
  idx = find(strcmp(Codes{iAirport},IATA));
  
  %pull out lat and lon as string (see above - backwards compatability)
  try
    Coords(iAirport,:) = [InfoArray.longitude_deg(idx), ...
                          InfoArray.latitude_deg( idx), ...
                          InfoArray.elevation_ft( idx)];
  catch; end
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
