function [] = kk_timetaken(Paths,DeSeas)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot time series of time taken over the instrument record
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load route and flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%all flight data
Flights = load([Paths.StoreDir,'/flight_data_',Paths.SourceIdentifier,'.mat']);
Flights = Flights.Flights;

%which flights are actually used?
load([Paths.StoreDir,'/routes_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find monthly-mean relative flight time for each plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time scale
[yy1,~,~]  = datevec(min(Flights.Date(Working.InSeason.All == 1)));
[yy2,~,~]  = datevec(max(Flights.Date(Working.InSeason.All == 1)));
TimeScale = datenum(yy1,-1:3:(yy2-yy1+1).*12,1);  %the -1 ensures we have seasonal chunks, by offsetting to start in DJF
clear yy1 yy2

%indices for planes
[PlaneNames,~,ic]= unique(Flights.PlaneID(Working.InSeason.All == 1));

%storage array
Store = NaN(numel(PlaneNames),numel(TimeScale),2); %2 is east or west
N     = Store; %how many flights contribute to each symbol?
for iEW=0:1
  for iPlane=1:1:numel(PlaneNames)
    ThisGroup = find(ic == iPlane & Working.Eastward(Working.InSeason.All == 1) == iEW);
    
    Store(iPlane,:,iEW+1) = bin2matN(1,Flights.Date(    ThisGroup), ...
                                       Working.tRel.All(ThisGroup), ...
                                       TimeScale,'@nanmean');
    N(    iPlane,:,iEW+1) = bin2matN(1,Flights.Date(ThisGroup), ...
                                       ones(size(ThisGroup)), ...
                                       TimeScale,'@nansum');
  end
end

%scale delays by overall median flight time
OverallMedian = 500 %temporary
Store = (Store.*OverallMedian) - OverallMedian;
% Store = Store./60; %minutes


%scale N flights to produce sensible point sizes
PointSizes = N; PointSizes(PointSizes == 0) = NaN;
PointSizes = log10(PointSizes);
PointSizes = floor(PointSizes.*10); 
PointSizes(PointSizes  == 0) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
clf
set(gcf,'color','w','position',[521 237 1277 557])

subplot = @(m,n,p) subtightplot (m, n, p, 0.04, 0.1, [0.07 0.2]);
Colours = flipud(cbrewer('qual','Set2',numel(PlaneNames)));

for iEW=1:2;
  
  
  %make panel
  subplot(2,1,iEW)
  axis([datenum([1994,2020],1,1),min(flatten(Store(:,:,iEW))).*1.4, max(flatten(Store(:,:,iEW))).*1.1])
  axis manual
  hold on
%   plot([0,1].*99e99,[0,0],'k-')  
  yvals = get(gca,'ylim');

  %get data
  Data = Store(:,:,iEW);
 
  
  
  %linear trend across all data
  a = nanmean(Data,1);
  Good = find(~isnan(a));
  p = polyfit(TimeScale(Good),a(Good),1);
  plot(TimeScale,polyval(p,TimeScale),'k-','linewi',2)





  %plot key
  %%%%%%%%%%%%%%%%%%

  %plane colours

  for iPlane=1:1:numel(PlaneNames)
    if iEW == 1;
      text(datenum(2021,15,1),max(yvals)-0.08.*range(yvals).*(iPlane-1+1),PlaneNames(iPlane),'color',Colours(iPlane,:),'fontweight','bold')
      plot(datenum(2021,05,1),max(yvals)-0.08.*range(yvals).*(iPlane-1+1),'ko','linewi',0.5,'markersize',10,'markerfacecolor',Colours(iPlane,:),'clipping','off')
    end
  end

  %symbol sizes
  for iSize=1:3:max(PointSizes(:))

    xpos = datenum(2021,05,1);
    ypos = max(yvals)-0.08.*range(yvals).*(numel(PlaneNames)-1+6+iSize./2);
    sz = 5+iSize.*6;
    scatter(xpos,ypos,sz,'k','linewi',0.5,'clipping','off','markerfacecolor',[1,1,1].*0.8)

    xpos = datenum(2021,15,1);
    thisn = round(10.^(iSize./10),1,'significant');
    if thisn > 1; plural = 's'; else; plural = ''; end
    text(xpos,ypos,texlabel([num2str(thisn),' flight',plural]),'verticalalignment','middle')
      
  end; clear iSize
  



  %plot individual flights
  %%%%%%%%%%%%%%%%%%%%%%%%%
  for iPlane=1:1:size(Store,1)
    for iTime=1:1:size(Store,2)

      %prepare to plot symbol
      Colour   = Colours(iPlane,:);
      Delay    = Data(iPlane,iTime); if isnan(Delay); continue; end
      PS       = 5+PointSizes(iPlane,iTime,iEW).*6;

      %plot symbol
      scatter(TimeScale(iTime),Delay,PS,[1,1,1].*0.2,        ...
              'markerfacecolor',Colour,'markerfacealpha',0.8,'linewi',0.5)

    end
  end


  %label plot
  xvals = get(gca,'xlim');
  if iEW == 1; Dir = 'EASTWARDS'; else Dir = 'WESTWARDS'; end
  text(max(xvals)-0.02.*range(xvals),min(yvals)+0.06.*range(yvals), ...
       Dir,'horizontalalignment','right')

  if DeSeas == 1;
    text(min(xvals)+0.02.*range(xvals),min(yvals)+0.06.*range(yvals), ...
         'deseasonalised','horizontalalignment','left','fontsize',10)
  end


  %tidy up
  datetick('keeplimits')
  if iEW == 1; set(gca,'xaxislocation','top'); end
  set(gca,'tickdir','out')
  set(gca,'xtick',datenum(1995:3:2020,1,1),'xticklabel',datestr(datenum(1995:3:2020,1,1),'yyyy'))
  ylabel('Delay [min]')
  box on
  plot(TimeScale,polyval(p,TimeScale),'k:','linewi',2)  


  drawnow
  




end
