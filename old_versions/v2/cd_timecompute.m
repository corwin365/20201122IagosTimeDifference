function [] = cd_timecompute(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot time series of time taken over the instrument record
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/12/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load route and flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlightData = load('data/flight_data.mat');
RouteData  = load('data/routes.mat');
Indices    = load('data/indices.mat');

%number of airports
NPorts = numel(RouteData.Airports);

%list of seasons
Seasons = fieldnames(Settings.Seasons);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute relative flight times, and split into Eward and Wward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlightData.Results.tRel  = NaN.*FlightData.Results.t;
FlightData.Results.Eward = NaN.*FlightData.Results.t;

for iDep=1:1:NPorts
  for iArr=1:1:NPorts
   for iSeason=1:1:numel(Seasons)
    
     %find all flights meeting these criteria
     Flights = squeeze(RouteData.Flights(iSeason,iDep,iArr,:));
     Flights = Flights(~isnan(Flights));
     if numel(Flights) == 0; continue; end

     %what is the relative time taken by each flights?
     FlightData.Results.tRel(Flights) = FlightData.Results.t(Flights)./nanmedian(FlightData.Results.t(Flights));
     
     %are we heading east or west?
     if ismember(RouteData.Airports{iDep},Settings.NA); E = 1; else E = 2; end  %1 eastward, 2 westward
     FlightData.Results.Eward(Flights) = E;
     
   end
  end
end
clear iDep iArr iSeason Flights NPorts E

%also compute the mean overall flight time, so we can convert our
%coefficients from relative time to minutes
MedianFlightTime = nanmedian(FlightData.Results.t(:))./60; %MINUTES


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pair the flights to produce round trips
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%let's take westward flights as normal (arbitrarily, but chosen
%because the planes used are mostly European), then try
%and pair each of them with the nearest eastward flight between
%the same pair of airports
%impose a maximum separation time to make sure we're not averaging 
%completely different climate index states
Settings.maxdt = 2; %days

%create an array fr the flight pairs
%deliberately too big - trimmed at end
FlightPairs = NaN(numel(Seasons), ...
                  numel(RouteData.Flights(1,:,:,:)), ...
                  2);


for iSeason=1:1:numel(Seasons)

  count = 0;
  for iArr=1:1:numel(RouteData.Airports)

    %is this airport in Europe?
    if ismember(RouteData.Airports{iArr},Settings.Eur)
      continue
    end
    
    %get all flights arriving at this airport
    Flights = flatten(RouteData.Flights(iSeason,:,iArr,:));
    Flights = Flights(~isnan(Flights));
    
    %ok. for each flight:
    for iFlight=1:1:numel(Flights)
     
      %find its origin
      Origin = FlightData.Results.Dep(Flights(iFlight));
      
      %and date
      Date = FlightData.Results.Date(Flights(iFlight));
      
      %now, find all flights going from here TO the origin
      FromHere = find(ismember(FlightData.Results.Dep,RouteData.Airports{iArr}));
      ToThere  = find(ismember(FlightData.Results.Arr(FromHere),Origin));

      %which flight is the nearest in time, and is it close enough?
      [dt,idx] = min(abs(FlightData.Results.Date(FromHere(ToThere)) - Date));
      if dt > Settings.maxdt; continue; end
      
      %ok, the flight pair is...
      FlightPair = [Flights(iFlight),FromHere(ToThere(idx))];

      %store it
      count = count+1;     
      FlightPairs(iSeason,count,:) = FlightPair;
      
    end
    
  end
end

%tidy up
clear count Date dt FlightPair Flights FromHere iArr idx
clear iFlight iSeason Origin ToThere

used = find(squeeze(nansum(FlightPairs,[1,3])) > 0);
FlightPairs = FlightPairs(:,used,:);
clear used




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('data/relative_times.mat','MedianFlightTime','FlightPairs')

Airports = FlightData.Airports; Results = FlightData.Results; Settings = FlightData.Settings;
save('data/flight_data.mat','Airports','Results','Settings')


end
