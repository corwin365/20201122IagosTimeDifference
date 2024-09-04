function process_routesplit(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify and quantify routes
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Identifying flight routes ')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% import data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flight data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_merged.mat'])

%get valid airport list
Airports = load([Settings.Paths.DataDir,'/',Settings.ID,'_airportinfo.mat']);
Airports.NA  = Airports.AirportData.Code(Airports.AirportData.Continent == 1);
Airports.Eur = Airports.AirportData.Code(Airports.AirportData.Continent == 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% produce a list of routes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%first, create a table to store route info in
VarNames = ["RouteID",   "Dep",   "Arr", "RoutePair",   "Direction","NFlights", "Flights"];
VarTypes = [ "double","string","string",    "double", "categorical",  "double",    "cell"];
RouteData = table('Size',[numel(Airports.NA)*numel(Airports.Eur).*2,numel(VarNames)], ...
                   'VariableTypes',VarTypes,                                          ...
                   'VariableNames',VarNames                                           );
clear VarNames VarTypes

%now, associate routes with IDs
Count = 0;
textprogressbar('Generating route table ')
for iX=1:1:numel(Airports.NA)
  for iY=1:1:numel(Airports.Eur)
    for Direction=[0,1]; %one away from Europe, one towards
      Count = Count+1;
      RouteData{Count,1} = Count;
      if Direction == 0; %eastward
        RouteData{Count,2} = string(Airports.NA{iX});
        RouteData{Count,3} = string(Airports.Eur{iY});
        RouteData{Count,4} = Count+1;
        RouteData{Count,5} = 'E';
      else %westward
        RouteData{Count,2} = string(Airports.Eur{iY});
        RouteData{Count,3} = string(Airports.NA{iX});
        RouteData{Count,4} = Count-1;
        RouteData{Count,5} = 'W';
      end
    end
  end
  textprogressbar(iX./numel(Airports.NA).*100)
end
clear Count iX iY
textprogressbar('!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% associate flights with routes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get lists of airports
a = table2array(RouteData(:,2)); A = table2array(FlightData(:,1));
b = table2array(RouteData(:,3)); B = table2array(FlightData(:,2));

textprogressbar('Associating flights with routes ')
for iRoute=1:1:size(RouteData,1)

  %get list of flights on this route
  ThisArr = find(strcmp(a(iRoute),A));
  ThisDep = find(strcmp(b(iRoute),B));
  ThisRoute = intersect(ThisDep,ThisArr);


  %store information
  RouteData.NFlights(iRoute) = numel(ThisRoute);
  RouteData.Flights{ iRoute} = ThisRoute;

  textprogressbar(iRoute./size(RouteData,1).*100)
end; 
clear clear a b A B iRoute ThisArr ThisDep ThisRoute
textprogressbar('!')

%drop completely empty routes
Good = find(RouteData.NFlights ~= 0);
RouteData = RouteData(Good,:);
clear Good

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([Settings.Paths.DataDir,'/',Settings.ID,'_routes.mat'],'RouteData')


disp('--------------------------')
disp('Flight routes identified ')
disp('--------------------------')

end

