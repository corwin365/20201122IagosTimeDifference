function [] = jj_planes_and_indices(Paths,IndicesToPlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot time series of individual aircraft points and climate indices
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/10, modified 2021/02/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, 0.01, 0.05, 0.2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%all flight data
load([Paths.StoreDir,'/flight_data_',Paths.SourceIdentifier,'.mat']);
clear Settings

%which flights are actually used?
load([Paths.StoreDir,'/routes_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat']);


%indices
Indices = load([Paths.StoreDir,'/indices_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat']);
Indices.Daily.Ranges = Indices.Ranges;
Indices.Daily.WiW = Indices.WhichIsWhich;
Indices = Indices.Daily;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subset data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Results.PlaneID = Flights.PlaneID(Working.Used == 1);
Results.Date    = Flights.Date(   Working.Used == 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find and plot the number of flights per month for each individual plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%individual planes
[Planes,ia,ic]= unique(Results.PlaneID);

%time scale
[yy1,~,~]  = datevec(min(Results.Date));
[yy2,~,~]  = datevec(max(Results.Date));
Scale = datenum(yy1,1:1:(yy2-yy1+1).*12,1);
clear yy1 yy2

%storage array
Store = NaN(numel(Planes),numel(Scale));


%grid
for iPlane=1:1:numel(Planes)
  ThisPlaneDates = Results.Date(find(ic == iPlane));
  Store(iPlane,:) = bin2matN(1,ThisPlaneDates,ones(size(ThisPlaneDates)),Scale,'@nansum');
end

%force alphabetical order
Planes = flipud(Planes);
Store = Store(end:-1:1,:);

%% plot as shaded filled cumulative contours
subplot(3+numel(IndicesToPlot),1,[1:3])
hold on

Store2 = cumsum(Store,1);
Store2 = cat(1,zeros(1,size(Store2,2)),Store2);

Colours = cbrewer('qual','Set2',numel(Planes));

for iPlane=1:1:numel(Planes)
  
  %get data
  x = [Scale,Scale(end:-1:1)];
  y = [Store2(iPlane,:),Store2(iPlane+1,end:-1:1)];
  
  %overinterpolate for plotting
  x = Scale(1):.5:Scale(end);
  y1 = interp1(Scale,Store2(  iPlane,:),x,'nearest');
  y2 = interp1(Scale,Store2(iPlane+1,:),x,'nearest');
  
  x = [x,x(end:-1:1)];
  y = [y1,y2(end:-1:1)];
  
  
  
%   y = smoothn(y,[1,5]);
  
  patch(x,y,Colours(iPlane,:),'edgecolor',[1,1,1].*0.3,'linewi',0.7)
  
  text(Scale(end)+180,10+13.*(iPlane-1),Planes(iPlane),'color',Colours(iPlane,:))
%   drawnow
  
end
ylabel('Number of flights')
axis([Scale(1) Scale(end) 0 max(nansum(Store,1)).*1.05])

box on
set(gca,'tickdir','out','xaxislocation','top')
set(gca,'xtick',datenum(1995:3:2020,1,1),'xticklabel',datestr(datenum(1995:3:2020,1,1),'yyyy'))


clear Colours Flights ia ic iPlane Planes RouteInfo 
clear y1 y2 x y ThisPlaneDates Store Store2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find and plot climate indices, for comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iIndex=1:1:numel(IndicesToPlot)

  %generate panel
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(3+numel(IndicesToPlot),1,3+iIndex)
  box on
  set(gca,'tickdir','out','color',[1,1,1].*0.7)
  set(gca,'xtick',datenum(1995:3:2020,1,1),'xticklabel',datestr(datenum(1995:3:2020,1,1),'yyyy'))
  axis([Scale(1) Scale(end) -1.1 1.1])  
  hold on

  %plot data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  t = Indices.Time;
  ToPlot = Indices.(IndicesToPlot{iIndex});
  plot(t,ToPlot,'w-','linewi',1.25)


  %label with original range
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %get values
  OR = Indices.Ranges.(IndicesToPlot{iIndex});
  RoundDP = 2;
  if strcmp(IndicesToPlot{iIndex},'Time') ~= 1;
    text(Scale(end)+180,-1,num2str(round(   OR(1),RoundDP)),'fontsize',9);
    text(Scale(end)+180, 0,num2str(round(mean(OR),RoundDP)),'fontsize',9);
    text(Scale(end)+180, 1,num2str(round(   OR(2),RoundDP)),'fontsize',9);
  else
    text(Scale(end)+180,-1,datestr(min( OR),'mm/yyyy'),'fontsize',9);
    text(Scale(end)+180, 0,datestr(mean(OR),'mm/yyyy'),'fontsize',9);
    text(Scale(end)+180, 1,datestr(max( OR),'mm/yyyy'),'fontsize',9);
  end

  %mark if raw or deseasonalised
  Which = Indices.WiW.(IndicesToPlot{iIndex});
  if strcmp(Which,'DS'); text(Scale(5),-0.8,'DS','color','w','fontsize',10); end


  %tidy up
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ylabel(IndicesToPlot{iIndex})
  if iIndex ~= numel(IndicesToPlot); set(gca,'xticklabel',[]); end
  drawnow

end

