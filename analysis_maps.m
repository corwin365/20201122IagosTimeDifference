function analysis_coverageandindices(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%meta information - map of routes
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/06/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Map of flight routes by season')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep flight maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RouteMapFile = [Settings.Paths.DataDir,'/',Settings.ID,'_routemap.mat'];

if ~exist(RouteMapFile,'file') | Settings.Choices.FlightPathMap.Regenerate == 1;

  disp('No map file exists or regeneration has been requested, generating')

  %load data
  load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])
  Data = FlightData;
  clear OverallMedianTimes RouteData FlightData


  %prep map
  Map.Lon = Settings.Choices.FlightPathMap.LonScale;
  Map.Lat = Settings.Choices.FlightPathMap.LatScale;
  [xi,yi] = meshgrid(Map.Lon,Map.Lat);

  FlightMaps = NaN(size(Data,1),numel(Map.Lon),numel(Map.Lat));

  parfor iFlight = 1:1:size(Data,1)
    if ~ismissing(Data.FilePath(iFlight))
      if exist(Data.FilePath(iFlight),'file')
        
        %load data
        ThisFlight = rCDF(Data.FilePath(iFlight));

        %interpolate to constant timescale
        OldTime = ThisFlight.UTC_time;
        NewTime = min(OldTime):1:max(OldTime);
        [~,idx] = unique(OldTime); %there are some duplicate points in the data
        ThisFlight.lat = interp1(OldTime(idx),ThisFlight.lat(idx),NewTime,'nearest');
        ThisFlight.lon = interp1(OldTime(idx),ThisFlight.lon(idx),NewTime,'nearest');

        %bin onto map
        FlightMaps(iFlight,:,:) = bin2mat(ThisFlight.lon,ThisFlight.lat,ones(size(ThisFlight.lat)),xi,yi,'@nansum')';
      end
    end
  end


  %split by season and direction
  NSeasons = numel(Settings.Seasons.List);
  Store = NaN(numel(Settings.Choices.Directions),NSeasons,numel(Map.Lon),numel(Map.Lat));
  for iDir=1:1:numel(Settings.Choices.Directions)
    for iSeason=1:1:NSeasons
      idx = find(table2array(Data.InSeasons(:,iSeason)) == 1 & Data.Direction == Settings.Choices.Directions{iDir});
      Store(iDir,iSeason,:,:) = nansum(FlightMaps(idx,:,:),1);
    end
  end
  Map.Map = Store;


  %save
  save(RouteMapFile,'Map')
  clear xi yi iFlight ThisFlight FlightMaps iSeason idx Store Data NSeasons iDir OldTime NewTime

else
  load(RouteMapFile);
end
clear RouteMapFile

%compute a min,max and median path in terms of latitude at each longitude
Limits = NaN(numel(Settings.Choices.Directions),numel(Settings.Seasons.List),size(Map.Lon,2),3); %the 3 is min,median,max lat

for iDir=1:1:numel(Settings.Choices.Directions)
  for iSeason=1:1:numel(Settings.Seasons.List)

    tracks = squeeze(Map.Map(iDir,iSeason,:,:));
    for iLon=1:1:size(tracks,1)
      if Map.Lon(iLon) < -70; continue; end
      if Map.Lon(iLon) > -5; continue; end
      idx = find(tracks(iLon,:) ~= 0);
      if numel(idx) == 0; continue; end
      
      Limits(iDir,iSeason,iLon,1) = Map.Lat(idx(1));
      Limits(iDir,iSeason,iLon,2) = Map.Lat(closest(cumsum(tracks(iLon,:)),sum(tracks(iLon,:))./2));
      Limits(iDir,iSeason,iLon,3) = Map.Lat(idx(end));
    end
  end
end
clear iDir iSeason tracks idx iLon 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep mean winds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WindMapFile = [Settings.Paths.DataDir,'/',Settings.ID,'_windmap.mat'];
if ~exist(WindMapFile,'file') | Settings.Choices.FlightPathMap.Regenerate == 1;

  disp('No wind file exists or regeneration has been requested, generating')

  %get climatological winds
  x = -180:1.5:180; y = -90:1.5:90; t = datenum(2020,1,1:1:366); %year will be ignored
  [xi,yi,ti] = meshgrid(x,y,t);
  stop
  BG = get_context(xi,yi,'TimePoints',ti,'Wind',true,'Era5_Clim',true,'Pressure',250); 
  clear xi yi ti

  %average into seasons
  WindMaps = NaN(numel(Settings.Seasons.List),numel(x),numel(y));
  for iSeason = 1:1:numel(Settings.Seasons.List)
    idx = Settings.Seasons.(Settings.Seasons.List{iSeason});
    WindMaps(iSeason,:,:) = mean(BG.U(:,:,idx),3)';
  end

  Wind.Lon = x;
  Wind.Lat = y; 
  Wind.U   = WindMaps;

  %save
  save(WindMapFile,'Wind');
  clear x y BG WindMaps iSeason idx



else
  load(WindMapFile)
end

%all_windmap.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot map for each season and direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get background map
x = -180:0.25:180; y = -90:0.25:90;
[xi,yi] = meshgrid(x,y);
BG = get_context(xi,yi,'SurfaceImage',true);%,'SurfaceImage_Image','land_ocean_ice');
clear xi yi

Map.Map(Map.Map == 0) = NaN;

% figure('Position',[-1919,-951,1920,970])
clf
subplot = @(m,n,p) subtightplot (m, n, p, 0.01, 0.06,[0.01 0.1]);

k = 0;
for iDir=[1,3] %omitting roundtrips
  for iSeason=1:1:numel(Settings.Seasons.List)


    %generate panel
    k = k+1;
    h = subplot(3,numel(Settings.Seasons.List),k);
    colormap(h,cbrew('YlOrRd',15));

    % %prep basemap
    % m_proj('ortho','lat',55','long',-50,'radius',45);
    m_proj('satellite','lat',50,'long',-50,'radius',1);
    % m_proj('albers','lat',[30 75],'lon',[-92 20]);
    m_image(x,y,BG.SurfaceImage);
    hold on

    %plot flights for season
    m_pcolor(Map.Lon,Map.Lat,log10(squeeze(Map.Map(iDir,iSeason,:,:)./60)'));

    %plot limits
    Lims = squeeze(Limits(iDir,iSeason,:,:));
    Good = find(~isnan(Lims(:,2)));
    if iDir == 1; Colour = 'k'; else; Colour = 'b'; end
    Lims = squeeze(Limits(iDir,iSeason,:,:));
    Good = find(~isnan(Lims(:,2)));
    % m_plot(Map.Lon(Good),smoothn(Lims(Good,1)',31),'-','linewi',1,'color',Colour)
    m_plot(Map.Lon(Good),smoothn(Lims(Good,2)',31),'-','linewi',2,'color',Colour)
    % m_plot(Map.Lon(Good),smoothn(Lims(Good,3)',31),'-','linewi',1,'color',Colour)

    %tidy up
    m_grid('xtick',-180:10:180,'ytick',-90:10:90,'xticklabel',{},'yticklabel',{});
    m_coast('color','k');
    caxis([0 3.5])

    %title
    if iDir == 1; title(Settings.Seasons.List{iSeason},'fontsize',24); end

    %done!
    drawnow

  end
end



%colours


cb = colorbar('eastoutside','position',[0.90 0.50 0.02 0.33],'ytick',[0:1:3],'yticklabel',10.^(0:1:3));
cb.Label.String = 'Total Flight-minutes';


LineLon = Map.Lon;
clear BG cb h iDir iSeason Map



clear WindMapFile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot wind maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSeason=1:1:numel(Settings.Seasons.List)


  %generate panel
  k = k+1;
  if strcmpi(Settings.Seasons.List{iSeason},'All'); continue; end %don't plot annual-mean wind
  h = subplot(3,numel(Settings.Seasons.List),k);
  colormap(h,cbrew('RdBu',60));

  % %prep basemap
  % m_proj('ortho','lat',55','long',-50,'radius',45);
  m_proj('satellite','lat',50,'long',-50,'radius',1);
  hold on

  %plot flights for season
  m_contourf(Wind.Lon,Wind.Lat,squeeze(Wind.U(iSeason,:,:))',-50:1:50,'edgecolor','none');

  %plot limits
  for iDir=[1,3]
    if iDir == 1; Colour = 'k'; else; Colour = 'b'; end
    Lims = squeeze(Limits(iDir,iSeason,:,:));
    Good = find(~isnan(Lims(:,2)));
    % m_plot(LineLon(Good),smoothn(Lims(Good,1)',31),'-','linewi',1,'color',Colour)
    m_plot(LineLon(Good),smoothn(Lims(Good,2)',31),'-','linewi',2,'color',Colour)
    % m_plot(LineLon(Good),smoothn(Lims(Good,3)',31),'-','linewi',1,'color',Colour)

  end

  %tidy up
  m_grid('xtick',-180:10:180,'ytick',-90:10:90,'xticklabel',{},'yticklabel',{});
  m_coast('color','k');
  caxis([-1,1].*40)

  %done!
  drawnow

end

cb = colorbar('eastoutside','position',[0.90 0.1 0.02 0.23]);
cb.Label.String = 'Climatological Zonal Wind [m/s]';

