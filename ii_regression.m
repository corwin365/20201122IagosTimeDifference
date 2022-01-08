function [] = ii_regression(Paths,Settings)

Indices  = Settings.Indices;
Seasons  = Settings.Seasons;
Colours  = Settings.Colours;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%carry out regression analysis, and plot resulrs
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2022/01/08
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
Reg.Est = NaN(3,numel(SeasonNames),numel(Settings.Indices)); %3 is number of directions
Fields = {'SE','T','P','N','R2'};
for iF=1:1:numel(Fields); Reg.(Fields{iF}) = Reg.Est; end
clear Fields

%do the work
for iSeason=1:1:numel(SeasonNames)

  for iDir=1:1:3 %eastward,westward,roundtrip

    %find data in season and direction
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
    Reg.Est(iDir,iSeason,:) = Coefs(2:end,1);
    Reg.SE( iDir,iSeason,:) = Coefs(2:end,2);
    Reg.T(  iDir,iSeason,:) = Coefs(2:end,3);
    Reg.P(  iDir,iSeason,:) = Coefs(2:end,4);
    Reg.N(  iDir,iSeason,:) = mdl.NumObservations;
    Reg.R2( iDir,iSeason,:) = mdl.Rsquared.Adjusted;


  end; clear iDir
end; clear iSeason
clear ThisGroup TotalFlightTime OverallMedian Beta Coefs mdl 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
clf
set(gcf,'color','w')
k = 0;


for iSeason=1:1:numel(SeasonNames)
  for iDir=[1,3,2]




    %prepare panel
    k = k+1;
    subplot(numel(SeasonNames),3,k)
    cla
    hold on; box off; grid off;
 
    plot([-999,999],[1,1].*numel(Settings.Indices)+0.5,'w-','linewi',4)        
    if iDir ~= 3; set(gca,'ytick',1:1:numel(Settings.Indices),'yticklabel',Settings.Indices); ylabel(SeasonNames{iSeason})
    else          set(gca,'ytick',1:1:numel(Settings.Indices),'yticklabel',{})
    end;
    for iY=1:1:numel(Settings.Indices); plot([-1,1].*40,[1,1].*iY,'-','color',Colours.(Settings.Indices{iY}),'linewi',0.5); end; clear iX
    for iX=-100:20:100; plot([1,1].*iX,[-999,999],'-','color',[1,1,1].*0.7,'linewi',0.25); end; clear iX
    if iSeason == numel(SeasonNames); xlabel('Delay [minutes]'); end;
    if iDir == 2; set(gca,'yaxislocation','right'); end

    axis([-40 40 0.5 numel(Settings.Indices)+0.5])
    set(gca,'ydir','reverse')

    for iIndex=1:1:numel(Settings.Indices)

      %get the value
      Value = Reg.Est(iDir,iSeason,iIndex);

      %is it statistically significant?
      Sig   = Reg.P(iDir,iSeason,iIndex);
      if Sig < 0.05; Colour = Colours.(Settings.Indices{iIndex}); else Colour = 'w'; end

      %indicate coefficient standard error
      StErr = Value + [-1,1]  .* Reg.SE(iDir,iSeason,iIndex);
      plot(StErr,[1,1].*iIndex,'k-','linewi',1); hold on
      plot([1,1].*StErr(1),[-1,1].*0.5+iIndex,'k-','linewi',1.5)
      plot([1,1].*StErr(2),[-1,1].*0.5+iIndex,'k-','linewi',1.5)    



      %plot data
      plot(Reg.Est(iDir,iSeason,iIndex),iIndex,'o', ...
           'color','k','markerfacecolor',Colour,'markersize',10)


      title(round(Reg.R2(iDir,iSeason,1),2))
      drawnow


    end; clear iIndex
  end; clear iSeason




end; clear iDir
stop