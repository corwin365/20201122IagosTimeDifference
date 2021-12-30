function [] = dd_routesplit(Paths,Settings,Airports)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%split the data into routes and "seasons", to facilitate later analysis
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/12/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load dataset of individual flights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flights
FlightData = load([Paths.StoreDir,'/flight_data_',Paths.SourceIdentifier,'.mat'],'Flights');
FlightData = FlightData.Flights;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each season (plus all-data), create a set of arrays corresponding to FlightData, called 'Working'
% field 'Used' will say if the data is used in that season
% field 'tRel' will tlel the relative flight time
% field 'Route' will associate it to a specific route, handled in a cell array called 'Routes'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cell array handling routes
RouteInfo = cell(numel(FlightData.t),7);  %this is way more than the number of possible routes expected,
                                          % but the number of routes can't be MORE than the number of flights, 
                                          % and we can trim it to size later


%top-level of struct (all-dataset properties)
Working = struct();
Working.Used  = ones(numel(FlightData.t),1);  %assume included by default
Working.Route = Working.Used.*NaN;            %unique route number
Working.tRel  = Working.Used.*NaN;            %relative flight time to other flights on route
Working.Pair  = Working.Used.*NaN;            %return flight paired with this outgoing flight

%make data outside our time range as universally bad
OutOfRange = find(FlightData.Date <= Settings.TimeRange(1) ...
                | FlightData.Date >= Settings.TimeRange(2));
Working.Used(OutOfRange) = 0;
disp('------')
disp([num2str(numel(OutOfRange)),' flights discarded as out of time range'])
clear OutOfRange

%then define season-level structures
Seasons = fieldnames(Settings.Seasons);
for iSeason=1:1:numel(Seasons);
  Working.(Seasons{iSeason}).Used  = Working.Used.*0;  %assume flights are NOT in the season until proved otherwise
  Working.(Seasons{iSeason}).Route = Working.Used.*NaN; 
  Working.(Seasons{iSeason}).tRel  = Working.Used.*NaN;   
end; clear iSeason

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% assign each flight to a route
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%identify all unique DEP and ARR airports
UniqueAirports = unique(cat(1,FlightData.Dep,FlightData.Arr)); 


%assign a unique numbers to each route
RouteCount = 0;
for iArr=1:1:numel(UniqueAirports)
  for iDep=1:1:numel(UniqueAirports)

    %find all flights with this ARR and DEP
    RouteCount = RouteCount+1;
    ThisArr  = find(ismember(FlightData.Arr,UniqueAirports{iArr}));
    ThisDep  = find(ismember(FlightData.Dep,UniqueAirports{iDep})); 
    ThisPair = intersect(ThisArr,ThisDep); 

    %identify the endpoints
    RouteInfo{RouteCount,1} = RouteCount;    
    RouteInfo{RouteCount,2} = UniqueAirports{iDep};
    RouteInfo{RouteCount,3} = UniqueAirports{iArr};
    RouteInfo{RouteCount,4} = numel(ThisPair);

    %is the flight eastward or westward?
    if ismember(UniqueAirports{iDep},Airports.NA); E = 1; else E = 2; end  %1 eastward, 2 westward
    RouteInfo{RouteCount,5} = E;

    %and tie the flights on this route back to the orginal list
    if numel(ThisPair) ~= 0; Working.Route(ThisPair) = RouteCount; end

  end
end
clear iArr iDep RouteCount ThisArr ThisDep ThisPair UniqueAirports E


%drop insufficiently used routes
RouteN = cell2mat(RouteInfo(:,4));
Good = find(RouteN >= Settings.MinFlights);
Bad  = find(RouteN <  Settings.MinFlights);
RouteInfo = RouteInfo(Good,:);
for iBad=1:1:numel(Bad); 
  Working.Used( Working.Route == Bad(iBad)) = 0; 
end
disp([num2str(sum(RouteN(Bad))),' flights dropped due to insufficiently busy routes at whole-dataset level'])
clear Good Bad RouteN iBad iSeason
disp('------')

%duplicate route identifiers down to seasons
for iSeason=1:1:numel(Seasons); Working.(Seasons{iSeason}).Route = Working.Route; end; clear iSeason

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optionally, deseasonalise travel times on each route
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Settings.DeseasFlights == 1;

  for iRoute=1:1:size(RouteInfo,1);

    %find flights on this route
    OnThisRoute = find(Working.Route == RouteInfo{iRoute,1} & Working.Used == 1);

    %find their day-of-year and flight duration
    DoY = floor(date2doy(FlightData.Date(OnThisRoute)));
    FT  = FlightData.t(OnThisRoute);

    %now, produce a rolling mean flight time for each day-of-the-year:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %1. duplicate data so we can get the end-of-year wraparound right
    DoY = [DoY-366;DoY;DoY+366];
    FT  = [FT;FT;FT];

    %2. sort by DoY
    [DoY,idx] = sort(DoY,'asc'); 
    FT = FT(idx);

    %3. average with rolling n-day window
    RollingAverage = NaN(366,1);
    for iDay=1:1:366;
      InRange = 366 + iDay+[-1,1].*0.5.*Settings.DeseasPeriod;
      DataInPeriod = find(DoY >= floor(InRange(1)) & DoY <= ceil(InRange(2)));
      RollingAverage(iDay) = nanmean(FT(DataInPeriod));
    end

    %finally, interpolate it to the data and remove it
    %note that this will NOT alter the disk-stored flight times, just the memory copy used in this file
    %as we do not save these modified data - it will affect anything computed in this routine though, i.e.
    %RELATIVE flight times overall and per-season

    %remove the seasonal variation from the signal, ADDING BACK THE ANNUAL MEAN to keep relative time stats used later sane.
    DoY = floor(date2doy(FlightData.Date(OnThisRoute)));
    FlightData.t(OnThisRoute) = FlightData.t(OnThisRoute) - interp1(1:1:366,RollingAverage,DoY) + nanmean(RollingAverage);


  end



end
clear iRoute AvPeriod DataInPeriod DoY FT iDay idx InRange OnThisRoute RollingAverage 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the median flight time for each route, and normalise flight times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Routes = unique(cell2mat(RouteInfo(:,1)));
for iRoute=1:1:numel(Routes)

  %get flights on route
  %%%%%%%%%%%%%%%%%%%%%  
  FlightsOnThisRoute = find(Working.Route == Routes(iRoute));

  %all data
  %%%%%%%%%%%%%%%%%%%

  %median over all flights
  RouteInfo{iRoute,6} = nanmedian(FlightData.t(FlightsOnThisRoute));

  %normalised individual flight times
  Working.tRel(FlightsOnThisRoute) = FlightData.t(FlightsOnThisRoute)./RouteInfo{iRoute,6};


  %by season
  %%%%%%%%%%%%%%%%%%%%%%

  DoY = floor(date2doy(FlightData.Date));
  SeasonMedians = NaN(numel(Seasons),1);
  for iSeason=1:1:numel(Seasons)

    %get indices of flights on this route in this system
    InThisSeason = find(ismember(DoY(FlightsOnThisRoute),Settings.Seasons.(Seasons{iSeason})));
    idx = FlightsOnThisRoute(InThisSeason); %this is much less typing

    %find median flight time
    SeasonMedians(iSeason) = nanmedian(FlightData.t(idx));

    %normalise flight times
    Working.(Seasons{iSeason}).tRel(idx) = FlightData.t(idx)./SeasonMedians(iSeason);

    %and flag as used
    Working.(Seasons{iSeason}).Used(idx) = 1;

  end


  %store
  RouteInfo{iRoute,7} = SeasonMedians;

end

%also remove the flights we marked as bad overall from the seasonal sets
for iSeason=1:1:numel(Seasons); Working.(Seasons{iSeason}).Used(Working.Used == 0) = 0; end
clear iRoute Routes FlightsOnThisRoute SeasonMedians DoY iSeason InThisSeason idx

%compute an overall median flight time for the whole dataset
OverallMedian = nanmedian(FlightData.t(Working.Used == 1));
disp(['Overall median flight time in retained dataset is ',num2str(OverallMedian./60),' minutes'])
disp('------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% drop any routes with:
% a. insufficient flights
% b. flights which are too long or too short 
% in a given season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%unpopular routes
%%%%%%%%%%%%%%%%%%%%

for iSeason=1:1:numel(Seasons)
  RtS = Working.(Seasons{iSeason}).Route; %routes this season
  RtS(Working.(Seasons{iSeason}).Used == 0) = NaN;
  uRtS = unique(RtS); uRtS = uRtS(~isnan(uRtS));

  Sigma = 0;
  for iRoute=1:1:numel(uRtS)
    ThisSeasonThisRoute = find(RtS == uRtS(iRoute));
    if numel(ThisSeasonThisRoute) < Settings.MinFlights;
      Working.(Seasons{iSeason}).Used(ThisSeasonThisRoute) = 0;
      Sigma = Sigma+numel(ThisSeasonThisRoute);
    end
  end
  disp([num2str(Sigma),' flights dropped due to insufficiently busy routes in ',Seasons{iSeason}])
end

disp('------')

%too long or too short
%%%%%%%%%%%%%%%%%%%%%%


for iSeason=1:1:numel(Seasons)

  %find flights that are either oddly long or oddly short relative to others in the season
  BadTime = find(Working.(Seasons{iSeason}).tRel < min(Settings.RelativeTime) ...
               | Working.(Seasons{iSeason}).tRel > max(Settings.RelativeTime));
  
  %remove them
  Working.(Seasons{iSeason}).Used(BadTime) = 0;

  %and tell the user about it
  disp([num2str(numel(BadTime)),' flights dropped due to being outside acceptable length variation in ',Seasons{iSeason}])
end
clear iSeason BadTime
disp('------')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pair the flights to produce round trips
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%let's take westward flights as normal (arbitrarily, but chosen
%because the planes used are mostly European), then try
%and pair each of them with the nearest eastward flight between
%the same pair of airports
%
%impose a maximum separation time to make sure we're not averaging 
%completely different climate index states


for iRoute=1:1:size(RouteInfo,1)

  %if the route is eastwards (i.e. towards europe), skip it
  if RouteInfo{iRoute,5} == 1; continue; end

  %find the return route between these two airports
  Origin      = RouteInfo{iRoute,2};  
  Destination = RouteInfo{iRoute,3};


  %now, find the return route (which may not exist, if it has too few flights and was discarded above)
  ReturnRoute= intersect(find(contains(RouteInfo(:,2),Destination)), ...
                         find(contains(RouteInfo(:,3),Origin)));
  if numel(ReturnRoute) == 0; continue; end %no return route found

  %now, find the times of these lists of flights
  OutList    = find(Working.Used == 1 & Working.Route == RouteInfo{     iRoute,1});
  ReturnList = find(Working.Used == 1 & Working.Route == RouteInfo{ReturnRoute,1});

  %ok. work through the outgoing list and find the closest flight temporally
  for iFlight=1:1:numel(OutList)
    [dt,idx] = min(abs(FlightData.Date(OutList(iFlight)) - FlightData.Date(ReturnList)));

    %is the time permissible?
    if dt > Settings.MaxDt; continue; end

    %store this pair
    Working.Pair(OutList(iFlight)) = ReturnList(idx);
  end
end
clear iRoute Origin Destination ReturnRoute OutList ReturnList iFlight dt idx

%duplicate data into seasonal arrays
for iSeason=1:1:numel(Seasons); Working.(Seasons{iSeason}).Pair = Working.Pair; end; clear iSeason



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([Paths.StoreDir,'/routes_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat'],'Working','RouteInfo','OverallMedian')


%and return
return