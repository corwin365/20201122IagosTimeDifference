function [] = jj_planes_and_indices(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot time series of individual aircraft points and climate indices
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, 0.02, 0.05, 0.2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the flight metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FRAC = Settings.Frac;
load('data/flight_data.mat')


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

subplot(8,1,[1:3])
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
  
  patch(x,y,Colours(iPlane,:),'edgecolor','none')
  
  text(datenum(2020,6,1),15.*(iPlane-1),Planes(iPlane),'color',Colours(iPlane,:))
  drawnow
  
end
ylabel('Number of flights')
axis([Scale(1) Scale(end) 0 max(nansum(Store,1)).*1.05])

box on
set(gca,'tickdir','out','xaxislocation','top')
set(gca,'xtick',datenum(1995:3:2020,1,1),'xticklabel',datestr(datenum(1995:3:2020,1,1),'yyyy'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find and plot climate indices, for comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except Settings Scale Results FRAC subplot

load('data/indices.mat')


Ind = {'ENSO','HadCRUT','QBO','NAM','TSI'};

for iIndex=1:1:numel(Ind)
  subplot(8,1,3+iIndex)
  plot(Results.Date,Indices.(Ind{iIndex}),'w-','linewi',1.25)
  set(gca,'xtick',datenum(1995:3:2020,1,1),'xticklabel',datestr(datenum(1995:3:2020,1,1),'yyyy'))
   
  axis([Scale(1) Scale(end) -1.1 1.1])  
  box on
  set(gca,'tickdir','out','color',[1,1,1].*0.7)
  text(datenum(2020,6,1),0.5,Ind{iIndex})
  if iIndex ~= numel(Ind); set(gca,'xticklabel',[]); end
  ylabel(Ind{iIndex})
  
end



