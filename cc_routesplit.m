function [] = cc_routesplit(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%split the data into routes, to facilitate later analysis
%data will also be split by "season" here
%
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load dataset of individual flights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = load('data/flight_data.mat');

%compute day-of-year of each flight, for "seasonal" splitting
DoY = floor(date2doy(Data.Results.Date));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% identify all unique DEP and ARR airports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Airports = unique(cat(1,Data.Results.Dep,Data.Results.Arr)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each pair of airports in each "season", identify all flights and store them as a group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get list of seasons
Seasons = fieldnames(Settings.Seasons);

%create storage arrays
Routes = NaN(numel(Seasons),numel(Airports),numel(Airports),numel(Data.Results.t)); %last dim will be massively reduced at end
MaxFlights = 0; %used for said reduction


% textprogressbar('Pairing routes ')
for iArr=1:1:numel(Airports)
  for iDep=1:1:numel(Airports)

    %find all flights with this ARR and DEP
    ThisArr     = find(ismember(Data.Results.Arr,Airports{iArr}));
    ThisDep     = find(ismember(Data.Results.Dep,Airports{iDep})); 
    ThisPairAll = intersect(ThisArr,ThisDep); %the 'all' is as opposed to "seasonal"
    
    %do we have enough flights to be useful?
    if numel(ThisPairAll) < Settings.MinFlights; continue; end
    
    %now go over "seasons"
    for iSeason = 1:1:numel(Seasons);
      
      %find all days in this "season"
      ThisPairSeason = ThisPairAll(find(ismember(DoY(ThisPairAll),Settings.Seasons.(Seasons{iSeason}))));

    
      %do we have enough flights to be useful?
      if numel(ThisPairSeason) < Settings.MinFlights; continue; end
      
      %what is the relative time taken by each flights?
      RelativeTime = Data.Results.t(ThisPairSeason)./nanmedian(Data.Results.t(ThisPairSeason));
      
      %hence, which should we exclude?
      Bad = find(RelativeTime < Settings.RelativeTime(1) ...
               | RelativeTime > Settings.RelativeTime(2));
      
      %if any were too fast or too slow, remove them
      %then check again that we have enough flights
      if numel(Bad) > 0
        ThisPairSeason(Bad) = [];
        if numel(ThisPairSeason) < Settings.MinFlights; continue; end
      end
      
      %we have enough. store them!
      Routes(iSeason,iDep,iArr,1:numel(ThisPairSeason)) = ThisPairSeason;
      
      %is this the largest set?
      if numel(ThisPairSeason) > MaxFlights; MaxFlights = numel(ThisPairSeason); end
    end
  end;
%   textprogressbar(iDep./numel(Airports).*100)
end;
% textprogressbar('!')
clear ThisArr ThisDep iArr iDep ThisPairAll RelativeTime Bad ThisPairSeason

%drop empty elements
Routes = Routes(:,:,:,1:MaxFlights);


%we now have an array of DEP x ARR, where the third dimension lists
%all flights in the dataset between these airports

%drop any combinations where there are no flights each way
%these are the airports which never meet our min-flights criterion

Total = squeeze(sum(~isnan(Routes),[1,4]));
UsedDep = find(sum(Total,2) > 0);
UsedArr = find(sum(Total,1) > 0);
UsedAirports = unique([UsedDep;UsedArr']);

clear Total UsedDep UsedArr MaxFlights Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finally, retain:
%1. a list of valid ARR and DEP airports
%2. the number of flights between each pair each way
%3. a list of the individual flight indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Airports = Airports(UsedAirports);
Flights  = Routes(:,UsedAirports,UsedAirports,:);
NFlights = sum(~isnan(Flights),4);

clear Routes UsedAirports


save('data/routes.mat')
