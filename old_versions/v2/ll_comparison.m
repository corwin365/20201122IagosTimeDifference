function [] = ll_comparison(Settings)

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

%list of seasons
Seasons = fieldnames(Settings.Seasons);

%median time and flight pairs
A =  load('data/relative_times.mat');
Pairs = A.FlightPairs;
MedianTime = A.MedianFlightTime;
clear A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find some interesting percentiles of each var in the data
%hardcode for now, but should factor up to master if this code 
%ends up useful
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Percentiles = [2.5,18,50,82,97.5];
Variables   = {'U','V','T'};
NFlights    = numel(FlightData.Results.tRel);
Store = NaN(numel(Variables), ...
            NFlights,         ...
            numel(Percentiles));

for iFlight=1:1:NFlights
  for iVar=1:1:numel(Variables)
    TS = FlightData.Results.Paths.(Variables{iVar}); TS = TS{iFlight};
    Store(iVar,iFlight,:) = prctile(TS,Percentiles);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scatter each variable against flight times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf,'color','w')

nbins_x = 25; nbins_y = 30;
tRel = (FlightData.Results.tRel.*MedianTime)-MedianTime;


for iVar=1:1:numel(Variables)
  for iPrc=1:1:numel(Percentiles)

    %generate subplot
    subplot(numel(Variables),numel(Percentiles),iPrc+numel(Percentiles).*(iVar-1))

    %generate density plot
    x = linspace(min(tRel),max(tRel),nbins_x);
    y = linspace(prctile(squeeze(Store(iVar,:,iPrc)),1),prctile(squeeze(Store(iVar,:,iPrc)),99),nbins_y);
    [xi,yi] = meshgrid(x,y);
    zz = bin2mat(tRel,squeeze(Store(iVar,:,iPrc)),ones(numel(tRel),1),xi,yi,'@nansum');
    zz = zz./nansum(zz(:)).*100;

    %plot data
    pcolor(xi,yi,zz); shading flat
    colormap(cbrewer('seq','YlOrRd',32))
    caxis([0 2])
    axis square

    if iPrc == 1; ylabel(Variables{iVar}); end
    if iVar == 1; title(Percentiles(iPrc)); end

    drawnow
  end
end