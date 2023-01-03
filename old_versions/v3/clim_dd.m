function [] = clim_dd(Paths,Settings)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%maps of routes taken, for the HAdCRUT-only analyses
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2022/01/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Seasons  = Settings.Seasons;


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


%get season names, and add 'all' 
SeasonNames = fieldnames(Seasons);
SeasonNames{end+1} = 'All';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load comparator records
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load HadCRUT
HadCRUT = rCDF([Paths.Indices,'/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc']);
HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);


%interpolate to same dates as the flights
sz = size(HadCRUT.tas_mean);
ts = reshape(HadCRUT.tas_mean,sz(1)*sz(2),sz(3))';
ts = interp1(HadCRUT.MatlabTime,ts,Flights.Date);
HadCRUT.T = reshape(ts',sz(1),sz(2),numel(Flights.Date));

%and remove unwanted fields (remember we have the time in the flights struct)
H.lon = HadCRUT.longitude; H.lat = HadCRUT.latitude;
H.T = HadCRUT.T;
HadCRUT = H; clear H;


disp('ERA5 to add later')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, split data into directions and seasons, and correlate each box temporally
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Corrs.HadCRUT = NaN(3,numel(SeasonNames),sz(1),sz(2));
Corrs.ERA5    = Corrs.HadCRUT; %to do

Corrs.PHadCRUT = Corrs.HadCRUT;



for iDir=1:1:3 %1 is eastward, 2 is westward, 3 is round-trip


  for iSeason=1:1:numel(SeasonNames)

    %select flights in this season and direction
    %also load their flight times


    if iDir == 1;
      ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                     & Working.Eastward                        == true);
      TotalFlightTime = Working.tRel.(SeasonNames{iSeason}); 
      TotalFlightTime = TotalFlightTime(ThisGroup);
      OverallMedian = OverallMedians.(SeasonNames{iSeason}).East;
    elseif iDir == 2;
      ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                     & Working.Eastward                        == false);
      TotalFlightTime = Working.tRel.(SeasonNames{iSeason}); 
      TotalFlightTime = TotalFlightTime(ThisGroup);
      OverallMedian = OverallMedians.(SeasonNames{iSeason}).West;
    elseif iDir == 3;
      ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                     & ~isnan(Working.Pair));
      TotalFlightTime = Working.tRel.(SeasonNames{iSeason}); 
      TotalFlightTime = TotalFlightTime(ThisGroup)+TotalFlightTime(Working.Pair(ThisGroup));
      OverallMedian = (OverallMedians.(SeasonNames{iSeason}).East+OverallMedians.(SeasonNames{iSeason}).West)/2;
    end

    %convert flight time to delay in minutes
    TotalFlightTime = TotalFlightTime.*OverallMedian;
    if iDir == 3; TotalFlightTime = TotalFlightTime-OverallMedian.*2;
    else          TotalFlightTime = TotalFlightTime-OverallMedian;
    end
    TotalFlightTime = TotalFlightTime./60;

    %produce reduced datasets on the same timescale
    H = HadCRUT.T(:,:,ThisGroup);

    %may get some NaNs for the paired routes - remove to be safe
    Good = find(~isnan(TotalFlightTime));
    H = H(:,:,Good); TotalFlightTime = TotalFlightTime(Good);

    %and correlate each gridbox
    for iX=1:1:sz(1)
      for iY=1:1:sz(2)

        [r,p] = corrcoef(H(iX,iY,:),TotalFlightTime);
        if numel(r) < 4; continue; end %not enough data to correlate
        Corrs.HadCRUT( iDir,iSeason,iX,iY) = r(2);
        Corrs.PHadCRUT(iDir,iSeason,iX,iY) = p(2);
        


      end
    end
    clear iX iY r H TotalFlightTime ThisGroup OverallMedian p Goos



  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% and map them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%prepare figure
% figure
% % set(gcf,'color','w','position',[250 265 1493 646])
subplot = @(m,n,p) subtightplot (m, n, p, 0.01, 0.05, [0.02,0.08]);
clf
k = 0;

ColourMap = flipud(cbrewer('div','RdYlBu',15)); 


for iDir=[1,3,2] %eastward,round trip, westward
  for iSeason=1:1:numel(SeasonNames)

    %prepare panel
    k = k+1;
    subplot(3,numel(SeasonNames),k)
    cla

%     m_proj('satellite','lon',-40,'lat',58,'rad',1,'rot',5);
    m_proj('Miller Cylindrical','lon',[-1,1].*180)

    %plot surface image and coasts
    hold on


    %get map data, and plot
    Data = squeeze(Corrs.HadCRUT( iDir,iSeason,:,:))';
    P    = squeeze(Corrs.PHadCRUT(iDir,iSeason,:,:))'; 
    Data(P > 0.05) = NaN;

    m_pcolor(HadCRUT.lon,HadCRUT.lat,Data); shading flat
%     m_contourf(HadCRUT.lon,HadCRUT.lat,Data,size(ColourMap,1),'edgecolor','none'); shading flat

    %tidy up panel
    m_coast('color',[1,1,1].*0.3);
    m_grid('xticklabel',[],'yticklabel',[],'linestyle',':','color',[1,1,1].*0.4); 
    colormap(ColourMap)
    caxis([-1,1].*0.25)


    %contextual additions
    if iDir == 1; title(SeasonNames{iSeason},'fontsize',24); end

    if iDir == 2 && iSeason == 1
      cb = colorbar('westoutside','position',[0.93 0.1 0.015 0.3]);
      set(cb,'ytick',-0.2:0.1:0.2)
      ylabel(cb,'r')
      set(cb,'TickLength',0.11)
    end

    drawnow


  end
end

