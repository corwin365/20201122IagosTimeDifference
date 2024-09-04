function process_roundtrips(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalise flight durations 
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Producing round trips ')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% import data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flight data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_merged.mat'])

%load route information
load([Settings.Paths.DataDir,'/',Settings.ID,'_routes.mat'])

%join the datasets to provide the information we need in a handier form
MergedData = outerjoin(FlightData,RouteData,'Keys',{'Dep','Arr'}, ...
                                            'LeftVariables', {'FlightIndex','Dep','Arr','Date','t'}, ...
                                            'RightVariables',{'RouteID','RoutePair','Direction'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pair flights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%create an array to store one possible pair for each flight.
FlightPairs = NaN(size(FlightData,1),2);

%find the best valid pair for each flight. Most won't get one.
textprogressbar('Pairing flights ')
for iFlight=1:1:size(FlightData,1)

  %update progress
  if mod(iFlight,100) == 1;textprogressbar(iFlight./size(FlightData,1).*100); end  

  %use WESTWARD flights as our basis to avoid duplicates (arbitrary choice)
  if MergedData.Direction(iFlight) ~= 'W'; continue; end

  %find all flights on opposing route
  PossiblePairs = find(MergedData.RouteID == MergedData.RoutePair(iFlight));
  if numel(PossiblePairs) == 0;continue; end

  %find closest difference in time to all possible pairs (in days)
  [dt,idx] = min(abs(MergedData.Date(iFlight) - MergedData.Date(PossiblePairs)));
  
  %is the pair close enough in time to be useful?
  if dt > Settings.Choices.MinPairDistance; continue; end

  %ok, we have a valid pair. Store it.
  FlightPairs(iFlight,:) = [iFlight,PossiblePairs(idx)];

end; 
textprogressbar(100);textprogressbar('!')
clear iFlight dt idx PossiblePairs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ok. now we're going to create a new category of route ('R', for 
%'round trip') and create flights using this route trip from the 
%merged pairs. Fiddly now, but means we can just treat them as 
%normal flights in all downstream code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%add 'return flights' to the master flight list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first, identify all flights in the paired dataset
idxPaired = find(~isnan(sum(FlightPairs,2)));
FlightPairs = FlightPairs(idxPaired,:);

%now extract the output (Eur -> NA) and inbound (NA -> Eur pairs). Remember these can have negative time seperation,
%so the distinction between the two types is purely for deduplicative purposes
PairedFlights.Out = MergedData(FlightPairs(:,1),:);
PairedFlights.In  = MergedData(FlightPairs(:,2),:);

%replace the time in the 'Out' array with the sum of the round trip time
PairedFlights.Out.t = PairedFlights.Out.t + PairedFlights.In.t;

%replace the type of trip with 'R'
PairedFlights.Out.Direction(:) = 'R';

%store the original separate flight data, useful for debugging
PairedFlights.Out.OriginalW = PairedFlights.Out.FlightIndex;
PairedFlights.Out.OriginalE = PairedFlights.In.FlightIndex;

%replace the indices with new numbers
PairedFlights.Out.FlightIndex = transpose(1:1:size(PairedFlights.Out,1))+max(FlightData.FlightIndex);

%wipe the route IDs, as we're about to replace these
PairedFlights.Out.RouteID(  :) = NaN;
PairedFlights.Out.RoutePair(:) = NaN;


clear idxPaired FlightPairs

%generate a list of the "new routes" in the same format as the real routes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get list of possible airports
Arr = unique(PairedFlights.Out.Arr);
Dep = unique(PairedFlights.Out.Dep);

%generate storage table, same format as the routes table
VarNames = ["RouteID",   "Dep",   "Arr", "RoutePair",   "Direction","NFlights", "Flights"];
VarTypes = [ "double","string","string",    "double", "categorical",  "double",    "cell"];
NewRoutes = table('Size',[numel(Dep)*numel(Arr),numel(VarNames)], ...
                   'VariableTypes',VarTypes,                         ...
                   'VariableNames',VarNames                          );
clear VarNames VarTypes

%fill it
Count  = 0;
for iDep=1:1:numel(Dep)
  for iArr=1:1:numel(Arr)

    %simple information
    Count = Count+1;
    NewRoutes.RouteID(  Count) = max(RouteData.RouteID)+Count;
    NewRoutes.Dep(      Count) = Dep(iDep);
    NewRoutes.Arr(      Count) = Arr(iArr);
    NewRoutes.RoutePair(Count) = NaN;
    NewRoutes.Direction(Count) = 'R';

    %number and list of flights
    ListOfFlights = find(strcmp(PairedFlights.Out.Dep,Dep(iDep)) & strcmp(PairedFlights.Out.Arr,Arr(iArr)));
    NewRoutes.NFlights(Count) = numel(ListOfFlights);
    NewRoutes.Flights{Count} = PairedFlights.Out.FlightIndex(ListOfFlights);

  end
end

%drop empty routes
Good = find(NewRoutes.NFlights > 0);
NewRoutes = NewRoutes(Good,:);

%add it to the main routes table
RouteData = [RouteData;NewRoutes];


clear iDep iArr Arr Dep Count Good ListOfFlights


%now, apply the "new route" IDs to the "new flight" records, and store
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iFlight=1:1:size(PairedFlights.Out,1)

  idx = find(NewRoutes.Dep == PairedFlights.Out.Dep(iFlight) & NewRoutes.Arr == PairedFlights.Out.Arr(iFlight));
  PairedFlights.Out.RouteID(   iFlight) = NewRoutes.RouteID(idx);
  PairedFlights.Out. RoutePair(iFlight) = NewRoutes.RouteID(idx);  

end


%we'll need to add some empty columns to the original data so we can add on 'OriginalW' and 'originalE' for the merged data
VarNames = ["OriginalW","OriginalE"];
VarTypes = [ "double","double"];
Add = table('Size',[size(MergedData,1),numel(VarNames)], ...
            'VariableTypes',VarTypes,                    ...
            'VariableNames',VarNames                     );
Add = standardizeMissing(Add,0);
MergedData = [MergedData,Add];
clear VarNames VarTypes Add

%now, add the new round trips to the master list of flights
MergedData = [MergedData;PairedFlights.Out;];
clear idx iFlight  PairedFlights NewRoutes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% merge everything together into one result table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%merge the MergedData and FlightData arrays and tidy up. The goal here is 
%to copy over the variables present in single but not paired data, so it will
%also take some tidying
CombinedFlights = outerjoin(MergedData,FlightData);

%drop unneeded duplicate columns. These arise because the roundtrip flights don't have all the same information as the one-ways, e.g. filepath etc
CombinedFlights = removevars(CombinedFlights,...
                             {'FlightIndex_FlightData','Date_FlightData','t_FlightData','Dep_FlightData','Arr_FlightData'});

%remove the round trip info from the raw data, as we've now made the roundtrip data
CombinedFlights = removevars(CombinedFlights,{'RoutePair'});

%and rename the remaining ones to their original names
CombinedFlights = renamevars(CombinedFlights, ...
                            {'FlightIndex_MergedData','Dep_MergedData','Arr_MergedData','Date_MergedData','t_MergedData'}, ...
                            {'FlightIndex','Dep','Arr','Date','t'});


disp('Roundtrips integrated')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlightData = CombinedFlights;

save([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_inc_roundtrips.mat'],'RouteData','FlightData')

disp('--------------------------')
disp('Round trips produced      ')
disp('--------------------------')

end

