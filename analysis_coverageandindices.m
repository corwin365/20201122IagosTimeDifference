function analysis_coverageandindices(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%meta information - data coverage and indices
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/05/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Data coverage and indices')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])

%load indices. Arbitrarily take the first direction, as they are unlagged
A = load([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat']);
RangeStore    = A.RangeStore.(   Settings.Choices.Directions{1});
DateIndices   = A.DateIndices.(  Settings.Choices.Directions{1});
FlightIndices = A.FlightIndices.(Settings.Choices.Directions{1});
clear A

clear FlightIndices OverallMedianTimes RouteData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot monthly data distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
clf
subplot = @(m,n,p) subtightplot (m, n, p, 0.005, [0.05,0.08], 0.08);
subplot(3+numel(Settings.Indices.List),1,[1,2,3])


%create a monthly timescale
[y,m,~] = datevec(FlightData.Date);
[~,idx] = unique(1000.*y +m);
Years = y(idx); Months = m(idx);
TimeScale = NaN(numel(Years)+1,1);
for iTime=1:1:numel(idx); TimeScale(iTime) = datenum(Years(iTime),Months(iTime),1); end
TimeScale = sort(TimeScale); TimeScale(end) = TimeScale(end-1)+28;
clear y m idx iTime

%convert IAGOS-only 'round trips' into one group
idx = find(~isnan(FlightData.OriginalW+FlightData.OriginalE));
Flag = zeros(size(FlightData,1),1);
for iRow=1:1:numel(idx)
  SourceA = FlightData.DataSource(FlightData.FlightIndex == FlightData.OriginalW(idx(iRow)));
  SourceB = FlightData.DataSource(FlightData.FlightIndex == FlightData.OriginalE(idx(iRow)));
  if strcmp(SourceA,'IAGOS') & strcmp(SourceA,SourceB); Flag(idx(iRow)) = 1; end
end
FlightData.Plane(Flag == 1) = 'Round';



%find all unique plane IDs
Planes = unique(FlightData.Plane);
Planes = Planes(~ismissing(Planes));


%hence, find all flights by each plane in each month
Results = NaN(numel(Planes),numel(TimeScale));
for iPlane=1:1:numel(Planes)
  ThisPlane = find(FlightData.Plane ==Planes(iPlane));
  for iMonth=1:1:numel(TimeScale)-1
    Results(iPlane,iMonth) = numel(inrange(FlightData.Date(ThisPlane),[TimeScale(iMonth+[0,1])]));
  end
end
clear iPlane iMonth

%find total number of flights excluding IAGOS round trips, and maximum for any single month
Sigma       = nansum(Results(find(~strcmp(Planes,'Round')),:),'all');
TheMax      = nanmax(Results(find(~strcmp(Planes,'Round')),:),[],'all');
TheFullMax  = nanmax(Results,[],'all');

%sort by total number of flights
[~,idx]= sort(nansum(Results,2),'desc');
Results = Results(idx,:); Planes = Planes(idx); clear idx


%now, plot the data
Colours = cbrew('Spectral',numel(Planes)); 
Basis = pi*20.^2; %area of a point representing the maximum value 

axis([minmax(DateIndices.Date)+[-1,1].*100,0,numel(Planes)+1])
hold on; box on
for iPlane=1:1:numel(Planes)

  plot(minmax(DateIndices.Date)+[-1,1].*100,[1,1].*iPlane,'color',[1,1,1].*0.3)

  for iMonth=1:1:numel(Months)

    %find number of points
    NPoints = Results(iPlane,iMonth);
    if NPoints == 0; continue; end


    %compute size of circle to plot
    r = ceil(sqrt(Basis.*NPoints./TheMax)./pi);
    plot(TimeScale(iMonth),1.*iPlane,'o','color','k','markersize',r,'markerfacecolor',Colours(iPlane,:))

  end

  %put total at end of row
  N = nansum(Results(iPlane,:));
  PC = nansum(Results(iPlane,:))./Sigma.*100;
  text(max(DateIndices.Date)+110,iPlane,num2str(N),'fontsize',11)
  if ~strcmp(Planes{iPlane},'Round');  text(max(DateIndices.Date)+420,iPlane,['(',num2str(round(PC,2,'significant')),'%)'],'fontsize',11); end
end
clear iPlane iMonth Colours
set(gca,'ytick',1:1:numel(Planes),'yticklabel',Planes,'ydir','reverse','xtick',datenum(1994:1:2023,1,1),'XTickLabel',{})


%overplot key
k = 0;
for iR=[1,10,20,30,40,50,TheMax,TheFullMax]
    r = ceil(sqrt(Basis.*iR./TheMax)./pi);
    plot(datenum(2015+k,1,1),-1,'o','color','k','markersize',r,'markerfacecolor',[1,1,1].*0.8,'clipping','off')
    text(datenum(2015+k,1,1),-2,num2str(iR),'clipping','off','horizontalalignment','center')
    k = k+0.65;    
end
text(datenum(2015,1,1)-150,-2,'Flights in month:','HorizontalAlignment','right')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot raw indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iIndex=1:1:numel(Settings.Indices.List)

  %get info
  ID    = Settings.Indices.List{iIndex};
  Range = RangeStore.(ID);
  Index = DateIndices.(ID);
  Time  = DateIndices.Date;


  %plot
  subplot(3+numel(Settings.Indices.List),1,3+iIndex)
  axis([minmax(Time)+[-1,1].*100,-1.1,1.1])
  hold on; box on  
  set(gca,'Color',Settings.Indices.Colours.(ID),'xtick',datenum(1994:1:2024,1,1));

  plot(Time,Index,'color','w','linewi',1)

  %label limits and index
  if strcmp(ID,'Time'); Label = datestr(Range(2)); else; Label = num2str(round(Range(2),5,'significant')); end
  text(max(Time)+110,1,Label,'fontsize',10,'clipping','off','verticalalignment','middle')
  if strcmp(ID,'Time'); Label = datestr(Range(1)); else; Label = num2str(round(Range(1),5,'significant')); end
  text(max(Time)+110,-1,Label,'fontsize',10,'clipping','off','verticalalignment','middle')
  ylabel(ID,'fontsize',18,'fontweight','bold')

  %tidy up
  if iIndex ~= numel(Settings.Indices.List); set(gca,'xticklabel',{})
  else;                                      datetick('x','keeplimits','keepticks')
  end


end; 
clear iIndex ID Range Index Time Label

