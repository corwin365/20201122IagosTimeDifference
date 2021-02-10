function [] = ee_routeplot(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the routes used in each direction
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create new figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
clf
set(gcf,'color','w','position',[101.4 347.8 1422.8 442])
subplot = @(m,n,p) subtightplot (m, n, p, 0.02, 0.02, [0.02 0.06]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load route and flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlightData = load('data/flight_data.mat');
RouteData  = load('data/routes.mat');
NPorts = numel(RouteData.Airports);

%get list of seasons
Seasons = fieldnames(Settings.Seasons);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% produce a map of flight density in each direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create arrays needed
[xi,yi] = meshgrid(Settings.Maps.Lon,Settings.Maps.Lat);
Maps = zeros([numel(Seasons),size(xi),2]);


%loop over  airports then seasons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textprogressbar('Mapping flights ')
for iDep=1:1:NPorts
  
  for iSeason=1:1:numel(Seasons)
  
    %find all flights DEPARTING from this airport in this season
    DepFlights = unique(RouteData.Flights(iSeason,iDep,:,:));
    
    %is this a westbound or eastbound flight?
    if   ismember(RouteData.Airports{iDep},Settings.NA); Dir = 1;
    else                                                 Dir = 2;
    end
    
    %for each flight, bin the route geographically then add to our array
    LonStore = [];
    LatStore = [];
    
    for iFlight=1:1:numel(DepFlights)
      if isnan(DepFlights(iFlight)); continue; end %no more flights from here
      
      
      LonStore = [LonStore,FlightData.Results.Paths.Lon{DepFlights(iFlight)}];
      LatStore = [LatStore,FlightData.Results.Paths.Lat{DepFlights(iFlight)}];
       
    end
    zz = bin2mat(LonStore,LatStore,ones(size(LonStore)),xi,yi,'@nansum');
    Maps(iSeason,:,:,Dir) = squeeze(Maps(iSeason,:,:,Dir)) + zz;
    
  end
  clear iFlight Lon Lat zz DepFlights iSeason LoNStore LatStore
    
  textprogressbar(iDep./NPorts.*100)
end; clear iDep
textprogressbar('!')



%remove 0s, and take log
Maps(Maps == 0) = NaN;
Maps = log10(Maps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load surface imagery
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,Map,~,Coasts] = topo_etc([min(xi(:)) max(xi(:))], ...
                               [min(yi(:)) max(yi(:))], ...
                               0,0,0,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% and plot the two maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
k = 0;

for iDir=1:2

  for iSeason=1:1:numel(Seasons)-1%omit all-seasons - redundant
    
    k = k+1;
    subplot(2,4,k)
    
    %produce map
    m_proj('albers','lat',[30 75],'lon',[-92 20]);
    m_image(Map.LonScale,Map.LatScale,Map.Map);
    hold on
    
    %add data     
    m_pcolor(xi,yi,squeeze(Maps(iSeason,:,:,iDir)));
    m_grid('fontsize',12);
    colormap(cbrewer('seq','YlOrBr','16'));
    
    %overlay coastlines
    for iCoast=1:1:numel(Coasts)
      m_plot(Coasts(iCoast).Lon,Coasts(iCoast).Lat,'-','linewi',0.5,'color',[1,1,1].*0)
      hold on
    end
    
    %set colour range
    caxis([0 3.5])
    

    drawnow

    
    % % %     m_coast('patch',[1,1,1].*0.9,'edgecolor','none');
% % %     hold on
% % %     
% % %     if iDir == 1; title([Seasons{iSeason},' Eward']);
% % %     else;         title([Seasons{iSeason},' Wward']);
% % %     end
% % %     
% % %    
    
    
  end
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% colourbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
cb = colorbar('eastoutside','position',[0.95 0.26 0.01 0.40]);
cb.Label.String = 'log_{10} (Flight-Minutes)';
