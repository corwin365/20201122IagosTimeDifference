function [] = clim_aa(Paths,Settings)

Indices    = Settings.Indices;
Seasons    = Settings.Seasons;
SubPeriods = Settings.Periods;
SigThresh  = Settings.SigThresh;

%sort periods by first year
[~,idx] = sort(SubPeriods(:,1));
SubPeriods = SubPeriods(idx,:);
clear idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%carry out regression analysis, and plot results
%this is a fork of the general version for apossible HadCRUT-only study.
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2022/01/09
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

%load indices
Index = load([Paths.StoreDir,'/indices_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat']);
FlightIndices = Index.Flight;
FlightIndices.RR = Index.Ranges.HadCRUT;
clear Index

%get season names, and add 'all' 
%list of seasons
SeasonNames = fieldnames(Settings.Seasons);
SeasonNames{end+1} = 'All';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute regressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create storage matrices
Reg = struct;
Reg.Est = NaN(3,size(SubPeriods,1),numel(SeasonNames),numel(Settings.Indices)); %3 is number of directions
Fields = {'SE','T','P','N','R2'};
for iF=1:1:numel(Fields); Reg.(Fields{iF}) = Reg.Est; end
clear Fields

%get year numbers
[yy,~,~] = datevec(Flights.Date);

%add in a whole-dataset fit
[yya,~,~] = datevec(Flights.Date(Working.InSeason.All == 1));
SubPeriods(end+1,:) = [min(yya) max(yya)];
clear yya

%do the work
for iPeriod=1:1:size(SubPeriods,1)
  for iSeason=1:1:numel(SeasonNames)

    for iDir=1:1:3 %eastward,westward,roundtrip

      %find data in season and direction
      if iDir == 1;
        ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1    ...
                        & Working.Eastward                       == true ...
                        & yy >= SubPeriods(iPeriod,1)                    ...
                        & yy <= SubPeriods(iPeriod,2)                    ...
                        );
        TotalFlightTime = Working.tRel.(SeasonNames{iSeason});
        TotalFlightTime = TotalFlightTime(ThisGroup);
        OverallMedian = OverallMedians.(SeasonNames{iSeason}).East;
      elseif iDir == 2;
        ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1     ...
                       & Working.Eastward                        == false ...
                        & yy >= SubPeriods(iPeriod,1)                     ...
                        & yy <= SubPeriods(iPeriod,2)                     ...
                        );
        TotalFlightTime = Working.tRel.(SeasonNames{iSeason});
        TotalFlightTime = TotalFlightTime(ThisGroup);
        OverallMedian = OverallMedians.(SeasonNames{iSeason}).West;
      elseif iDir == 3;
        ThisGroup = find(Working.InSeason.(SeasonNames{iSeason}) == 1 ...
                       & ~isnan(Working.Pair)                         ...
                        & yy >= SubPeriods(iPeriod,1)                 ...
                        & yy <= SubPeriods(iPeriod,2)                 ...
                        );
        TotalFlightTime = Working.tRel.(SeasonNames{iSeason});
        TotalFlightTime = TotalFlightTime(ThisGroup)+TotalFlightTime(Working.Pair(ThisGroup));
        OverallMedian = (OverallMedians.(SeasonNames{iSeason}).East+OverallMedians.(SeasonNames{iSeason}).West)/2;
      end

      %scale total flight time to "delay" metric
      TotalFlightTime = TotalFlightTime.*OverallMedian;
      if iDir == 3; TotalFlightTime = TotalFlightTime-OverallMedian.*2;
      else          TotalFlightTime = TotalFlightTime-OverallMedian;
      end
      TotalFlightTime = TotalFlightTime./60;

      %pull out indices to match the flighttime vectors, and mush into common matrix
      Beta = NaN(numel(TotalFlightTime),numel(Settings.Indices));
      for iIndex=1:1:numel(Indices)
        Index = FlightIndices.(Settings.Indices{iIndex});
        Beta(:,iIndex) = Index(ThisGroup);
      end; clear iIndex Index


      %fit a linear regression model
      mdl = fitlm(Beta,TotalFlightTime);

      %convert the coefficients to a matrix
      Coefs = table2array(mdl.Coefficients);

      %store the coefficients, the number of flights, and the R^2 value of the fit
      Reg.Est(iDir,iPeriod,iSeason,:) = Coefs(2:end,1);
      Reg.SE( iDir,iPeriod,iSeason,:) = Coefs(2:end,2);
      Reg.T(  iDir,iPeriod,iSeason,:) = Coefs(2:end,3);
      Reg.P(  iDir,iPeriod,iSeason,:) = Coefs(2:end,4);
      Reg.N(  iDir,iPeriod,iSeason,:) = mdl.NumObservations;
      Reg.R2( iDir,iPeriod,iSeason,:) = mdl.Rsquared.Adjusted;


    end; clear iDir
  end; clear iSeason
end; clear iPeriod
clear ThisGroup TotalFlightTime OverallMedian Beta Coefs mdl iPeriod yy

%remove all-data set from list of periods - will be handled manually below
SubPeriods = SubPeriods(1:end-1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scale the results so that they are **per degree of warming**
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%find rnage of index used in scaled coords
UsedRange = range(FlightIndices.HadCRUT);

%and range in raw units
RealRange = range(FlightIndices.RR);

%hence, compute scaling factor
ScaleFactor = UsedRange./RealRange; %i.e. one scaled unit is this many Kelvin of warming

%and do the scaling
Reg.Est = Reg.Est .*ScaleFactor;
Reg.SE  = Reg.SE  .*ScaleFactor;
clear ScaleFactor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare figure
% figure
clf
set(gcf,'color','w','position',[83 32 1810 964])
subplot = @(m,n,p) subtightplot (m, n, p, [0.02,0.03], 0.05, 0.1);
k = 0;


for iSeason=1:1:numel(SeasonNames);
  for iDir=[1,3,2];

   %prepare axes
   k = k+1;
   subplot(numel(SeasonNames),3,k)
   cla
   hold on; grid off


    %get data
    Est    = squeeze(Reg.Est(iDir,:,iSeason));
    P      = squeeze(Reg.P(  iDir,:,iSeason));
    SE     = squeeze(Reg.SE( iDir,:,iSeason)); 

    %plot some horizontal grid lines
    for iY=-100:5:100; plot([-1,1].*999,[1,1].*iY,'color',[1,1,1].*0.8,'linewi',0.5); end; clear iY

    %plot the all-date values as a background patch
    a = -999; b = 999;
    Colour = [255, 165, 0]./255;
    patch([a b b a a],Est(end)+[-1,-1,1,1,-1].*SE(end).*2,Colour,'edgecolor','none','facealpha',0.3)    
    patch([a b b a a],Est(end)+[-1,-1,1,1,-1].*SE(end),   Colour,'edgecolor','none','facealpha',0.8)
    %and a line
    plot([a,b],[1,1].*Est(end),'color','r')
    clear a b

    for iPeriod=1:1:size(SubPeriods,1)

% %       %get timerange
% %       TimeRange = [datenum(SubPeriods(iPeriod,1), 1, 1), ...
% %                    datenum(SubPeriods(iPeriod,2),12,31)];

      %is the data statistically significant?
      if P(iPeriod) < SigThresh; 
        Colour = 'k'; LineStyle = '-'; LineWi = 3;
      else; 
        Colour = [1,1,1].*0.7; LineStyle = '-'; LineWi = 3;
      end


      %plot data
      Top = Est(iPeriod)+SE(iPeriod); Bottom = Est(iPeriod)-SE(iPeriod);
      plot([1,1].*iPeriod,[Top,Bottom], ...
           'color',Colour,'linewi',LineWi,'linestyle',LineStyle);
      plot(iPeriod,Est(iPeriod),'o','color',Colour,'markerfacecolor',Colour,'markersize',6)


      %produce two-year identifiers of the start and end dates
      Start = sprintf('%02d',SubPeriods(iPeriod,1)-floor(SubPeriods(iPeriod,1)./100).*100);
      End   = sprintf('%02d',SubPeriods(iPeriod,2)-floor(SubPeriods(iPeriod,2)./100).*100);      


      %label data
      text(iPeriod,Bottom,Start, ...
           'fontsize',10.5, ...
           'horizontalalignment','center','verticalalignment','top');
      text(iPeriod,Top,End, ...
           'fontsize',10.5, ...
           'horizontalalignment','center','verticalalignment','bottom');
      clear Top Bottom Start End Colour LineStyle LineWi



    end; clear iPeriod

    %tidy up
%     datetick
    ylim([-1.65,1]+[min(Est-SE),max(Est+SE)]) %this is so the white line we later use to add the x-axis doens't cover the year labels
    xlim([0 size(SubPeriods,1)+1])

    %labelling
    if iDir ==1; ylabel('\Delta Delay [min/K]'); end
    if iSeason == 1;
      switch iDir; case 1; Title = 'Eastwards'; case 2; Title = 'Westwards'; case 3; Title = 'Round Trip'; end;
      title(Title,'fontsize',25)
    end; clear Title
  


    %draw my own axes
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    plot(xlims,[1,1].*min(ylims),'w','linewi',2)
    set(gca,'xtick',[]);
    plot([1,1].*max(xlims),ylims,'k'); 
    plot([1,1].*min(xlims),ylims,'k')




    drawnow



  end; clear iSeason
end; clear iDir
