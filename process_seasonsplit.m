function process_seasonsplit(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalise flight durations 
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Splitting data into seasons and normalising ')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% import data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flight and routing data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_inc_roundtrips.mat'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% split data by season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find day-of-year for each flight
DoY = floor(date2doy(FlightData.Date));

%find median flight time for each route in each season
MedianTimes        = NaN(size(RouteData,1),numel(Settings.Seasons.List));
OverallMedianTimes = struct();

for iSeason=1:1:numel(Settings.Seasons.List)

  %find all flights in each season
  InThisSeason = find(ismember(DoY,Settings.Seasons.(Settings.Seasons.List{iSeason})));

  %find an overall median for the season first
  OverallMedianTimes.(Settings.Seasons.List{iSeason}) = median(FlightData.t(InThisSeason));

  %pull out the relevant rows of the flight table for this season so we can work on them
  SeasonFlightData = FlightData(InThisSeason,:);

  %find the median flight time for each route in this season
  for iRoute=1:1:size(RouteData,1)

    %find flights on this route this season
    ThisRouteAndSeason = find(SeasonFlightData.RouteID == RouteData.RouteID(iRoute));

    %skip if we don't have enough
    if numel(ThisRouteAndSeason) < Settings.Choices.MinFlights; continue; end

    %find median flight time for this route and season
    MedianTimes(iRoute,iSeason) = median(SeasonFlightData.t(ThisRouteAndSeason));

  end



 
end;
clear iSeason InThisSeason SeasonFlightData iRoute ThisRouteAndSeason

%convert to table
MedianTimes = array2table(MedianTimes);
a = fieldnames(MedianTimes); a = a(1:numel(Settings.Seasons.List));
MedianTimes = renamevars(MedianTimes,a,string(Settings.Seasons.List));
clear a

%and store in the master route table
RouteData = [RouteData,mergevars(MedianTimes,string(Settings.Seasons.List),'MergeAsTable',true)];
RouteData = renamevars(RouteData,["Var1"],"MedianTravelTime");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalise flight times to season median
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NormalisedTimes = NaN(numel(DoY),numel(Settings.Seasons.List));

for iSeason=1:1:numel(Settings.Seasons.List)
  for iRoute=1:1:size(RouteData,1)


    %get median, and skip if NaN (must have been filtered out above)
    TheMedian = table2array(MedianTimes(iRoute,iSeason));
    if isnan(TheMedian); continue; end

    %get flights on this route in this season
    ThisRouteAndSeason = find(FlightData.RouteID == RouteData.RouteID(iRoute) ...
                            & ismember(DoY,Settings.Seasons.(Settings.Seasons.List{iSeason})));


    %normalise flight times
    NormalisedTimes(ThisRouteAndSeason,iSeason) = FlightData.t(ThisRouteAndSeason)./TheMedian;

  end
end
clear iSeason iRoute

%discard flights which are too long or short for their season
if Settings.Choices.MaxDeviation < 100;
  Min =    (1-Settings.Choices.MaxDeviation./100);
  Max = 1./(1-Settings.Choices.MaxDeviation./100);
  Bad = find(NormalisedTimes < Min | NormalisedTimes > Max);
  NormalisedTimes(Bad) = NaN;
end


%add normalised times to master flight table
aa = NormalisedTimes; %keep for later
NormalisedTimes = array2table(NormalisedTimes);
a = fieldnames(NormalisedTimes); a = a(1:numel(Settings.Seasons.List));
NormalisedTimes = renamevars(NormalisedTimes,a,string(Settings.Seasons.List));

FlightData = [FlightData,mergevars(NormalisedTimes,string(Settings.Seasons.List),'MergeAsTable',true)];
FlightData = renamevars(FlightData,["Var1"],"RelativeTravelTime");

%also store which season(s) each flight falls into
InSeasons = array2table(~isnan(aa));
a = fieldnames(InSeasons); a = a(1:numel(Settings.Seasons.List));
InSeasons = renamevars(InSeasons,a,string(Settings.Seasons.List));
FlightData = [FlightData,mergevars(InSeasons,string(Settings.Seasons.List),'MergeAsTable',true)];
FlightData = renamevars(FlightData,["Var1"],"InSeasons");

clear aaa NormalisedTimes InSeasons a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% convert each flight time to a notional delay in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%convert
RTT = FlightData.RelativeTravelTime;
for iSeason=1:1:numel(Settings.Seasons.List)
  RTT.(Settings.Seasons.List{iSeason}) = (RTT.(Settings.Seasons.List{iSeason})-1) .* OverallMedianTimes.(Settings.Seasons.List{iSeason});
end; clear iSeason

%store
FlightData.Delay = RTT;
clear RTT



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'],'RouteData','FlightData','OverallMedianTimes')

disp('--------------------------')
disp('Data split into seasons and normalised')
disp('--------------------------')

end

