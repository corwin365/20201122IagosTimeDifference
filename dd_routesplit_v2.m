function [] = dd_routesplit_v2(Paths,Settings,Airports)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%split the data into routes and "seasons", to facilitate later analysis
%
%version 2, as I was having trouble finding a bug in the original

%
%Corwin Wright, c.wright@bath.ac.uk, 2022/02/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load dataset of individual flights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flights
FlightData = load([Paths.StoreDir,'/flight_data_',Paths.SourceIdentifier,'.mat'],'Flights');
FlightData = FlightData.Flights;

%load list of "seasons"
SeasonNames = fieldnames(Settings.Seasons);

%add 'All' to this list - this will include all data
SeasonNames{end+1} = 'All';
Settings.Seasons.All = 1:1:366; %i.e. all days of the year

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create a series of arrays that map onto the flights
% and contain information we will use later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create struct
Working = struct();

%for each "season"...
%%%%%%%%%%%%%%%%%%%%%%%
for iSeason=1:1:numel(SeasonNames); 

  %create flag arrays for if the flight is included in this subset
  Working.Used.(SeasonNames{iSeason}) = true(size(FlightData.Dep));  %assume a flight is INCLUDED until removed

  %create arrays for relative flight time within the season, compared to all flights on that route in that direction
  Working.tRel.(SeasonNames{iSeason}) = NaN(size(FlightData.Dep));

end; clear iSeason

%and over all flights...
%%%%%%%%%%%%%%%%%%%%%%%%

%what route number is the flight on?
Working.RouteNumber = NaN(size(FlightData.Dep));

%is the flight eastward or westward?
Working.Eastward = false(size(FlightData.Dep));

%what is the eastward 'pair' of a given westward flight?
Working.Pair = NaN(size(FlightData.Dep));

%for paired flights, what is the difference between Eward and Wward flight times?
%this is based on the ideas of 10.1038/NCLIMATE2715
Working.PairDelta = NaN(size(FlightData.Dep));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% discard any flights outside our overall time range
% and divide data into 'seasons'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%overall out-of-range
OutOfRange = find(FlightData.Date <= Settings.TimeRange(1) ...
                | FlightData.Date >= Settings.TimeRange(2));

dd = floor(date2doy(FlightData.Date));
for iSeason=1:1:numel(SeasonNames)

  %find all dates *not* in this season
  DaysInSeason = Settings.Seasons.(SeasonNames{iSeason});
  OutofSeason = find(~ismember(dd,DaysInSeason));

  %add on the overall out-of-range days, and flag to ignore
  Bad = unique([OutOfRange;OutofSeason]);
  ToFlag = Working.Used.(SeasonNames{iSeason});
  ToFlag(Bad) = false;
  Working.Used.(SeasonNames{iSeason}) = ToFlag;

end; 

disp('--------')
disp(['Data divided into seasons; ',num2str(numel(OutOfRange)),' flights discarded as out of time range'])
clear dd OutOfRange iSeason DaysInSeason OutofSeason Bad ToFlag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% assign each flight to a route
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%identify all unique DEP and ARR airports
UniqueAirports = unique(cat(1,FlightData.Dep,FlightData.Arr)); 

%create an array to store route metadata
RouteInfo = cell(1,5); %The 1 will grow. The 6 will be: [RouteID,Origin,Destination,NFlights,Eastward,FlightMedians]

%assign a unique number to each route
RouteCount = 0;
for iArr=1:1:numel(UniqueAirports)
  for iDep=1:1:numel(UniqueAirports)
    if iArr == iDep; continue; end

    %find all flights with this ARR and DEP
    ThisArr  = find(ismember(FlightData.Arr,UniqueAirports{iArr}));
    ThisDep  = find(ismember(FlightData.Dep,UniqueAirports{iDep})); 
    ThisPair = intersect(ThisArr,ThisDep); 
    if numel(ThisPair) == 0; continue; end
    RouteCount = RouteCount+1;    

    %identify and store some route metadata
    RouteInfo{RouteCount,1} = RouteCount;    
    RouteInfo{RouteCount,2} = UniqueAirports{iDep};
    RouteInfo{RouteCount,3} = UniqueAirports{iArr};
    RouteInfo{RouteCount,4} = numel(ThisPair);

    %is the flight eastward or westward?
    if ismember(UniqueAirports{iDep},Airports.NA); 
      RouteInfo{RouteCount,5} = true;
    else                                           
      RouteInfo{RouteCount,5} = false;
    end

    %mark all individual flights on this route with their route number and direction
    Working.RouteNumber(ThisPair) = RouteCount;
    Working.Eastward(   ThisPair) = RouteInfo{RouteCount,5};

  end
end
clear iArr iDep  ThisArr ThisDep ThisPair UniqueAirports

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% drop insufficiently busy routes within each season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------')
for iSeason=1:1:numel(SeasonNames)
  Sigma = 0;

  for iRoute=1:1:RouteCount;

    %find all flights on this route in this season
    RouteFlights = find(Working.RouteNumber                 == iRoute ...
                      & Working.Used.(SeasonNames{iSeason}) == 1      );

    %if too few flights left on this route in this season...
    if numel(RouteFlights) < Settings.MinFlights;
    
      %... remove these flights from the valid dataset
      zz = Working.Used.(SeasonNames{iSeason});
      zz(RouteFlights) = false;
      Working.Used.(SeasonNames{iSeason}) = zz;
      %and increment count, for information
      Sigma = Sigma+numel(RouteFlights);
      
    end

  end;
  disp(['Flights dropped due to insufficiently busy routes in ',SeasonNames{iSeason},': ',num2str(Sigma)])


end;
clear iRoute iSeason zz RouteFlights Sigma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the median flight time for each season-route combination, and normalise flight times to this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iRoute=1:1:RouteCount;
  SeasonMedians = NaN(numel(SeasonNames),1);
  for iSeason=numel(SeasonNames):-1:1 %the negative ordering is to allow certain tests - see comment below

    %find all flights on this route in this season
    RouteFlights = find(Working.RouteNumber                 == iRoute ...
                      & Working.Used.(SeasonNames{iSeason}) == 1      );

    %find the median flight time
    SeasonMedians(iSeason) = median(FlightData.t(RouteFlights),'omitnan');

    %scale all relevant flights to this
    zz = Working.tRel.(SeasonNames{iSeason});
    zz(RouteFlights) = FlightData.t(RouteFlights) ./ SeasonMedians(iSeason);
    Working.tRel.(SeasonNames{iSeason}) = zz;
    
% %     %use this for tests where tRel seasonal and tRel all need to be the same - should be commented out otherwise
% %     zz = Working.tRel.All;
% %     zz(RouteFlights) = FlightData.t(RouteFlights) ./ SeasonMedians(iSeason);
% %     Working.tRel.All = zz; 

  end

  %store these medians as part of the route info
  RouteInfo{iRoute,6} = SeasonMedians;

end
clear iRoute SeasonMedians iSeason RouteFlights zz 

%also find a set of representative overall medians
for iSeason=1:1:numel(SeasonNames);
  OverallMedians.(SeasonNames{iSeason}).East = median(FlightData.t(Working.Eastward ==  true & Working.Used.(SeasonNames{iSeason}) == true),'omitnan');
  OverallMedians.(SeasonNames{iSeason}).West = median(FlightData.t(Working.Eastward == false & Working.Used.(SeasonNames{iSeason}) == true),'omitnan');
  OverallMedians.(SeasonNames{iSeason}).Both = median(FlightData.t(                            Working.Used.(SeasonNames{iSeason}) == true),'omitnan');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% drop any flights outside the defined useful range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------')
for iSeason=1:1:numel(SeasonNames)

  Times = Working.tRel.(SeasonNames{iSeason}); 
  Used  = Working.Used.(SeasonNames{iSeason});

  Bad = find(Times < min(Settings.RelativeTime) | Times > max(Settings.RelativeTime));
  Used(Bad) = false;

  Working.Used.(SeasonNames{iSeason}) = Used;

  disp(['Flights dropped as too fast or too slow in ',SeasonNames{iSeason},': ',num2str(numel(Bad))])
  
end
clear iSeason Times Used Bad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pair flights to produce composite round trips
%require both members to be valid members of the 'All' season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%loop by flights. It would be more efficient to loop by routes, but easier to introduce a bug.
for iFlight=1:1:numel(FlightData.Dep)

  %require outgoing flight to be a valid member of 'All'
  if ~Working.Used.All(iFlight); continue; end

  %require outgoing flight to be westbound
  if Working.Eastward(iFlight); continue; end

  %find all flights with flipped DEP and ARR
  OppositeRoute = find(strcmp(RouteInfo(:,3),FlightData.Dep{iFlight}) ...
                     & strcmp(RouteInfo(:,2),FlightData.Arr{iFlight}));
  OnOppositeRoute = find(Working.RouteNumber == OppositeRoute);

  %drop out if opposite route has no flights
  if numel(OnOppositeRoute) < 1; continue; end

  %drop any opposite-route flights which are not valid mmbers of 'All'
  OnOppositeRoute = OnOppositeRoute(Working.Used.All(OnOppositeRoute) == true);

  %now, find the closest in time to the incoming flight
  [dt,idx] = min(abs(FlightData.Date(iFlight)-FlightData.Date(OnOppositeRoute)));
  if numel(idx) ==0; continue; end %this can happen when using a time period reduced from the full set, due to the way we select OnOppositeRoute in time

  %if it's too far away in time, skip out
  if dt > Settings.MaxDt; continue; end

  %otherwise, store it as a pair
  Working.Pair(iFlight) = OnOppositeRoute(idx);

  %and store the travel-time difference between the pair
  Working.PairDelta(iFlight) = FlightData.t(iFlight) - FlightData.t(idx);

end
clear iFlight OnOppositeRoute dt idx OppositeRoute
disp('--------')
disp('Flights paired')
disp('--------')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% produce a version of 'All' which has the annual cycle removed
% this can be optionally swapped in using the DeseasFlights flag
% in place of the main dataset, but will always be made.
%
%needs tom be done separately for east and west
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iEast=1:2;

  if iEast == 1; EastLogical = true; else EastLogical = false; end

  %first, find the day-of-year and relative flight time of each used flight
  dd   = floor(date2doy(FlightData.Date(Working.Used.All == 1 & Working.Eastward == EastLogical)))';
  tRel = Working.tRel.All(              Working.Used.All == 1 & Working.Eastward == EastLogical)'-1; %the -1 is to remove the overall median

  %duplicate out over three years, to facilitate accurate smoothing at ends
  dd = [dd-366,dd,dd+366];
  tRel = [tRel,tRel,tRel];

  %now, generate a rolling time-mean
  TimeMean = NaN.*dd;
  for iDay=367:1:366*2
    InWindow = find(dd > iDay-Settings.DeseasPeriod./2 ...
                  & dd < iDay+Settings.DeseasPeriod./2);
    TimeMean(iDay) = mean(tRel(InWindow),'omitnan');
  end; clear iDay tRel InWindow dd


  %and cut down to one year again
  dd = 1:1:366;
  TimeMean = TimeMean(367:2*366);

  %interpolate to the times of the original flights
  ToRemove = interp1(dd,TimeMean,date2doy(FlightData.Date(Working.Eastward == EastLogical)));

  %and produce the DS timeseries
  Working.tRel.DS(Working.Eastward == EastLogical) = Working.tRel.All(Working.Eastward == EastLogical) - ToRemove;
  Working.Used.DS = Working.Used.All; %this will happen twice, but is a cheap operation so fine.

  %finally, store the seasonal cycle - may be independently useful
  SeasonalCycle.Cycle(iEast,:) = TimeMean;
  SeasonalCycle.DoY(  iEast,:) = 1:1:366;
  SeasonalCycle.East( iEast)   = EastLogical;

  clear dd TimeMean ToRemove EastLogical

end; clear iEast
Working.tRel.DS = Working.tRel.DS'; %undo earlier flip

%finally, do we want to replace the primary data with this?
if Settings.DeseasFlights == 1;
  disp('Switching in deseasonalised data for raw')
  Working.tRel.Raw = Working.tRel.All;
  Working.Used.Raw = Working.Used.All;
  Working.tRel.All = Working.tRel.DS;
  Working.tRel = rmfield(Working.tRel,'DS');
  Working.Used = rmfield(Working.Used,'DS');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rename 'Used' to 'InSeason' (I changed my mind at the end and don't want to go throguh to fix every case)
Working.InSeason = Working.Used; Working = rmfield(Working,'Used');

save([Paths.StoreDir,'/routes_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat'], ...
      'Working','RouteInfo','OverallMedians','SeasonalCycle')


%and return
return