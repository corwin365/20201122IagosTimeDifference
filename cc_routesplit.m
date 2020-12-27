function [] = cc_routesplit(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%split the data into routes, to facilitate later analysis
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load dataset of individual flights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = load('data/flight_data.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% identify all unique DEP and ARR airports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Airports = unique(cat(1,Data.Results.Dep,Data.Results.Arr)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each pair, identify all flights and store them as a group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Routes = NaN(numel(Airports),numel(Airports),numel(Data.Results.t)); %last dim will be massively reduced at end
MaxFlights = 0; %used for said reduction


textprogressbar('Pairing routes ')
for iArr=1:1:numel(Airports)
  for iDep=1:1:numel(Airports)
    
    %find all flights with this ARR and DEP
    ThisArr = find(ismember(Data.Results.Arr,Airports{iArr}));
    ThisDep = find(ismember(Data.Results.Dep,Airports{iDep})); 
    ThisPair = intersect(ThisArr,ThisDep);
    
    %do we have enough flights to be useful?
    if numel(ThisPair) < Settings.MinFlights; continue; end
    
    %we have enough. store them!
    Routes(iDep,iArr,1:numel(ThisPair)) = ThisPair;
    
    %is this the largest set?
    if numel(ThisPair) > MaxFlights; MaxFlights = numel(ThisPair); end

  end;
  textprogressbar(iDep./numel(Airports).*100)
end;
textprogressbar('!')
clear ThisArr ThisDep iArr iDep ThisPair

%drop empty elements
Routes = Routes(:,:,1:MaxFlights);

%we now have an array of DEp x ARR, where the third dimension lists
%all flights in the dataset between these airports

%drop any combinations where there are no flights each way
%these are the airports which never meet our min-flights criterion
Total = sum(~isnan(Routes),3);
UsedDep = find(sum(Total,2) > 0);
UsedArr = find(sum(Total,2) > 0);
UsedAirports = unique([UsedDep,UsedArr]);

clear Total UsedDep UsedArr MaxFlights Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finally, retain:
%1. a list of valid ARR and DEP airports
%2. the number of flights between each pair each way
%3. a list of the individual flight indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Airports = Airports(UsedAirports);
Flights  = Routes(UsedAirports,UsedAirports,:);
NFlights = sum(~isnan(Flights),3);

clear Routes UsedAirports


save('data/routes.mat')
