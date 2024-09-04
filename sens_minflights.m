clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sensitivity test of minimum flights per route
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/09/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Prefix = 'sens_minflights';
Settings.Values = 1:1:50;
Settings.Percentiles = [0:0.5:100]; 
Settings.NBoots = 10000;

Settings.BasisFile = 'all';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate files needed for main routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % for iValue=1:1:numel(Settings.Values)
% % % % % 
% % % % %   %define ID of output
% % % % %   OutString = [Settings.Prefix,num2str(Settings.Values(iValue))]
% % % % %   %load the master settings file
% % % % %   A = load(['data/',Settings.BasisFile,'.mat']);
% % % % % 
% % % % %   %change value
% % % % %   A.Choices.MinFlights = Settings.Values(iValue);
% % % % %   A.ID = OutString;
% % % % % 
% % % % %   %save to a new settings file
% % % % %   OutFile = ['data/',Settings.Prefix,num2str(Settings.Values(iValue)),'.mat'];
% % % % %   save(OutFile,'-struct','A');
% % % % % 
% % % % %   %also make copies of the data files
% % % % % 
% % % % %   copyfile(['data/',Settings.BasisFile,'_routes.mat'],           ['data/',OutString,'_routes.mat'] )
% % % % %   copyfile(['data/',Settings.BasisFile,'_flightinfo_merged.mat'],['data/',OutString,'_flightinfo_merged.mat'] )
% % % % %   copyfile(['data/',Settings.BasisFile,'_airportinfo.mat'],      ['data/',OutString,'_airportinfo.mat'] )
% % % % % end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BootStore = NaN(numel(Settings.Values),numel(Settings.Percentiles));

for iFile=1:1:numel(Settings.Values)

  %load the flight data file
  InFile = ['data/',Settings.Prefix,num2str(Settings.Values(iFile)),'_flightinfo_inc_roundtrips.mat'];
  if ~exist(InFile,'file'); clear InFile; continue; end
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
hold on

%plot data
for iPC=1:1:(numel(Settings.Percentiles)-1)./2
  x = [Settings.Values,Settings.Values(end:-1:1)];
  y = [BootStore(:,iPC);BootStore(end:-1:1,end-iPC+1)]';
  patch(x,y,[51,153,255]./255,'facealpha',0.04,'edgecolor','none')

  if any(Settings.Percentiles(iPC) == [0,2.5,18,50,82,97.5,100]);
    plot(Settings.Values,BootStore(:,iPC),      'k:','linewi',0.5)
    text(Settings.Values(end)+1,BootStore(end,iPC),[num2str(Settings.Percentiles(iPC)),'%'],'fontsize',6)

    plot(Settings.Values,BootStore(:,end-iPC+1),'k:','linewi',0.5)
    text(Settings.Values(end)+1,BootStore(end,end-iPC+1),[num2str(Settings.Percentiles(end-iPC+1)),'%'],'fontsize',6)
  end

end; clear iPc

plot(Settings.Values,BootStore(:,closest(Settings.Percentiles,50)),'k-o','linewi',1,'markerfacecolor','k','markersize',4)

%tidy up

box on
ylabel('Variance in flight time [%]')
xlabel('Minimum flights per route')
title('(b) Minimum Flights')
