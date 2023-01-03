function [] = clim_cc(Paths,Settings)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%maps of routes taken, for the HAdCRUT-only analyses
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2022/01/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load route and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%all flight data
Flights = load([Paths.StoreDir,'/flight_data_',Paths.SourceIdentifier,'.mat']);
Flights = Flights.Flights;

%which flights are actually used?
load([Paths.StoreDir,'/routes_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat']);


%get season names, and add 'all' 
%list of seasons
SeasonNames = fieldnames(Settings.Seasons);
SeasonNames{end+1} = 'All';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% produce geographic flight-density data for each season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Maps.Lon = -100:.15:30;
Maps.Lat = 20:.15:90;
% Maps.Lon = -100:5:30;
% Maps.Lat = 20:5:90;
[xi,yi] = meshgrid(Maps.Lon,Maps.Lat);

Maps.Data    = NaN(numel(SeasonNames),2,numel(Maps.Lon),numel(Maps.Lat)); %2 is east and west
Maps.Medians = NaN(numel(SeasonNames),2,numel(Maps.Lon));  %median latitude at each longitude


for iSeason=1:1:numel(SeasonNames)
  for iDir=1:2

    if     iDir == 1; FlightsToUse = find(Working.InSeason.(SeasonNames{iSeason}) == 1 & Working.Eastward == true);
    elseif iDir == 2; FlightsToUse = find(Working.InSeason.(SeasonNames{iSeason}) == 1 & Working.Eastward == false);
    end;
      
    LatList = []; LonList = [];


    for iFlight=1:1:numel(FlightsToUse)

      LatList = [LatList,Flights.Paths.Lat{FlightsToUse(iFlight)}];
      LonList = [LonList,Flights.Paths.Lon{FlightsToUse(iFlight)}];      

    end; clear iFlight 

    Maps.Data(   iSeason,iDir,:,:) =  bin2mat(LonList,LatList,ones(size(LonList)),xi,yi,'@nansum')';

    %also find median lat at each lon
    for iLon=1:1:numel(Maps.Lon)-1
      Maps.Medians(iSeason,iDir,iLon) = median(LatList(inrange(LonList,Maps.Lon([iLon,iLon+1]))),'omitnan');
    end; clear iLon

    clear FlightsToUse LatList LonList


  end; clear iDir
end; clear iSeason

%smooth out the medians a bit, as they can jump around between tracks in the md-Atlantic
Maps.Medians = smoothn(Maps.Medians,[1,1,11]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load topography, wind climatology, and prepare colour table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load surface imagery 
[~,~,Image] = topo_v2([-180 180],[-90,90],'Image','pale');

%load wind clima
ClimU = load([Paths.Indices,'/flightclima.mat']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare figure
% figure
set(gcf,'color','w','position',[259 20 1493 976])
subplot = @(m,n,p) subtightplot (m, n, p, 0.01, 0.05, [0.02,0.08]);
clf
k = 0;
ColourRange = [0 floor(log10(max(Maps.Data(:))))];
ColourMap = cbrewer('seq','YlOrRd',3.*(1+range(ColourRange))); 
ColourMap = ColourMap((3*1+1):end,:);

for iDir=[2,1,3] %in this case, '3' is the wind climatology
  for iSeason=1:1:numel(SeasonNames)

    %prepare panel
    k = k+1;
    h = subplot(3,numel(SeasonNames),k);
    cla

    m_proj('satellite','lon',-40,'lat',58,'rad',0.6,'rot',5);


    %plot surface image and coasts
    m_image(Image.Lon,Image.Lat,Image.Map);
    hold on




    %get map data, and plot
    if iDir < 3;
      Data = squeeze(Maps.Data(iSeason,iDir,:,:))';
      Data(Data == 0) = NaN;
      Data = log10(Data);
      m_pcolor(Maps.Lon,Maps.Lat,Data);
      clear Data


      %overplot median routing
      %'satellite' projection of m_map is buggy for this, so lot individual points
      for iPoint=inrange(Maps.Lon,[-70,-10])
        m_plot(Maps.Lon(iPoint),squeeze(Maps.Medians(iSeason,iDir,iPoint))','b.','markersize',3)
      end

      colormap(h,ColourMap)
      caxis(ColourRange)

    else

      
      Data = squeeze(ClimU.Clima.u(iSeason,:,:))';
      m_contourf(ClimU.Lon,ClimU.Lat,Data,-50:1:50,'edgecolor','none');
      colormap(h,flipud(cbrewer('div','RdBu',32)))
      caxis([-1,1].*50)

      %overplot median routing
      %'satellite' projection of m_map is buggy for this, so lot individual points
      for iPoint=inrange(Maps.Lon,[-70,-10])
        for jDir=1:2;
        m_plot(Maps.Lon(iPoint),squeeze(Maps.Medians(iSeason,jDir,iPoint))','b.','markersize',3)
        end
      end
    end

    %tidy up panel
    m_coast('color',[1,1,1].*0.3);
    m_grid('xticklabel',[],'yticklabel',[],'linestyle',':','color',[1,1,1].*0.4);


    %contextual additions
    if iDir == 2; title(SeasonNames{iSeason},'fontsize',24); end

    if iDir == 2 && iSeason == 1
      cb = colorbar('westoutside','position',[0.94 0.44 0.015 0.3]);
      set(cb,'ytick',0:1:5,'yticklabel',{'1','10','10^2','10^3','10^4','10^5'})
      ylabel(cb,'Flight-Minutes')
      set(cb,'TickLength',0.075)
    elseif iDir ==3 & iSeason == 1;
      cb = colorbar('westoutside','position',[0.94 0.15 0.015 0.3]);
      set(cb,'ytick',-50:20:50)
      ylabel(cb,'Climatological Zonal Wind [ms^{-1}]')
      set(cb,'TickLength',0.075)
    end

    drawnow


  end
end



