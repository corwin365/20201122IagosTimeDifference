function [] = ab_geolocate(Airports,Paths)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%geolocate airport codes
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load big file of airport metadata
InfoArray = readtable([Paths.StoreDir,'/airports.csv']);
InfoArray = InfoArray(:,[5,6,7,8,9,14,]); %lat,lon,elev,continent,country,iata code

%drop any that don't have an IATA code
[~,idx] = rmmissing(InfoArray.iata_code);
InfoArray(idx,:) = [];
clear idx

%create merged list of airports
Codes = [Airports.NA,Airports.Eur];

%create storage array for geolocation
%weird storage format is for backwards compatability with older code
Coords = cell(1,numel(Codes));

%fill this cell array
IATA = InfoArray.iata_code;
for iAirport=1:1:numel(Codes)

  %find code in InfoArray
  idx = find(strcmp(Codes{iAirport},IATA));
  
  %pull out lat and lon as string (see above - backwards compatability)
  Coords{iAirport} = [num2str(InfoArray.longitude_deg(idx)),' ',num2str(InfoArray.latitude_deg(idx))];
  
end
clear IATA idx iAirport


%save results
save([Paths.StoreDir,'/airport_codes_',Paths.SourceIdentifier,'.mat'],'Codes','Coords')
