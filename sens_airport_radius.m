clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sensitivity test of aircport exclusion radius
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/09/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Prefix = 'sens_mindist_';
Settings.Values = [0 0.5 1	1.5	2.4	3.7	5.6	8.7	13.3	20.5	31.6	48.7	75	115.5	177.8	273.8	421.7	500 649.4	700 1000];
Settings.Percentiles = [0:0.5:100]; 
Settings.NBoots = 10000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BootStore = NaN(numel(Settings.Values),numel(Settings.Percentiles));

for iFile=1:1:numel(Settings.Values)

  %load the flight data file
  InFile = ['data/sensitivity_AER/',Settings.Prefix,num2str(Settings.Values(iFile)),'_flightinfo_merged.mat'];
  if ~exist(InFile,'file'); clear InFile; continue; end
  load(InFile);  
  clear InFile


  %load the route file
  InFile = ['data/sensitivity_AER/',Settings.Prefix,num2str(Settings.Values(iFile)),'_routes.mat'];
  if ~exist(InFile,'file'); clear InFile FlightData; continue; end
  load(InFile);
  clear InFile

  %for each route, find the median travel time and normalise all flights to it, then store
  tStore = [];
  for iRoute=1:1:size(RouteData,1)
    t = FlightData.t(cell2mat(RouteData.Flights(iRoute,:)));
    tStore = [tStore;t./median(t)];
  end; clear iRoute RouteData FlightData t

  tStore = tStore*100; % convert to %

  %bootstrap the variance to get a confidence estimate range
  idx = randi(numel(tStore),[Settings.NBoots,numel(tStore)]);
  variance = var(tStore(idx),[],2);
  BootStore(iFile,:) = prctile(variance,Settings.Percentiles);
  clear idx variance tStore

end; clear iFile


%remove empty points (for testing)
Good = find(~isnan(BootStore(:,1)));
BootStore = BootStore(Good,:);
Settings.Values = Settings.Values(Good);
clear Good

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
hold on
axis([1e-1 1500 10 30])
patch([1.5e-1,0.5,0.5,1.5e-1],[10.04,10.04,29.96,29.96],'w','edgecolor','none')
% set(gca,'Layer','top')


%set first value to be just slight off the edge of the data
Settings.Values(1) = 0.1;


%plot data
for iPC=1:1:(numel(Settings.Percentiles)-1)./2
  x = [Settings.Values,Settings.Values(end:-1:1)];
  y = [BootStore(:,iPC);BootStore(end:-1:1,end-iPC+1)]';
  patch(x,y,[255,153,51]./255,'facealpha',0.04,'edgecolor','none')

  if any(Settings.Percentiles(iPC) == [0,2.5,18,50,82,97.5,100]);
    plot(Settings.Values,BootStore(:,iPC),      'k:','linewi',0.5)
    text(Settings.Values(end)+1,BootStore(end,iPC),[num2str(Settings.Percentiles(iPC)),'%'],'fontsize',6)

    plot(Settings.Values,BootStore(:,end-iPC+1),'k:','linewi',0.5)
    text(Settings.Values(end)+1,BootStore(end,end-iPC+1),[num2str(Settings.Percentiles(end-iPC+1)),'%'],'fontsize',6)
  end

end; clear iPc

plot(Settings.Values,BootStore(:,closest(Settings.Percentiles,50)),'k-o','linewi',1,'markerfacecolor','k','markersize',4)

%tidy up
set(gca,'xscale','log')
box on
ylabel('Variance in flight time [%]')
xlabel('Airport exclusion radius [km]')
title('(a) Minimum Distance')
set(gca,'xtick',10.^(-1:1:3),'xticklabel',[0,10.^(0:1:3)])
