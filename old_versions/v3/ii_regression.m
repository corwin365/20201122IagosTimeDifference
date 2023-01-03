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
    Reg.Est(iDir,iSeason,:) = Coefs(2:end,1);        % Coefficient estimates for each corresponding term in the model.
    Reg.SE( iDir,iSeason,:) = Coefs(2:end,2);        % Standard error of the coefficients.
    Reg.T(  iDir,iSeason,:) = Coefs(2:end,3);        % t-statistic for each coefficient to test the null hypothesis that the corresponding coefficient is zero against the alternative that it is different from zero, given the other predictors in the model. Note that tStat = Estimate/SE.
    Reg.P(  iDir,iSeason,:) = Coefs(2:end,4);        % p-value for the t-statistic of the hypothesis test that the corresponding coefficient is equal to zero or not.
    Reg.N(  iDir,iSeason,:) = mdl.NumObservations;   % Ronseal
    Reg.R2( iDir,iSeason,:) = mdl.Rsquared.Adjusted; % adjusted coefficient of determination


  end; clear iDir
end; clear iSeason
clear ThisGroup TotalFlightTime OverallMedian Beta Coefs mdl 


%multiply values by 2 (as we defined the indices on a -1 to +1 range)
Reg.Est = Reg.Est .* 2;
Reg.SE  = Reg.SE  .* 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
clf
set(gcf,'color','w','position',[243.4 32.2 1030.4 757.6])
k = 0;
subplot = @(m,n,p) subtightplot (m, n, p, [0.03,0.03], 0.1, 0.1);


for iSeason=1:1:numel(SeasonNames)
  for iDir=[1,3,2] %eastward,roundtrip,westward

    %prepare panel
    k = k+1;
    subplot(numel(SeasonNames),3,k)
    cla
    hold on; box on; grid off;
    axis([-18 18 0.5 numel(Settings.Indices)+0.5])
    set(gca,'ydir','reverse')

    %handle x axis and vertical gridlines
    if     iSeason == 1; set(gca,'xaxislocation','top');
    else                 set(gca,'xaxislocation','bottom'); 
    end
    if iSeason ~= 1 & iSeason ~= numel(SeasonNames); set(gca,'xtick',[]); 
    else; xlabel('Delay [minutes]');
    end
    plot([1,1].*min(get(gca,'xlim')),minmax(get(gca,'ylim')),'w-','linewi',2)
    plot([1,1].*max(get(gca,'xlim')),minmax(get(gca,'ylim')),'w-','linewi',2)    
    for iX=-100:5:100; plot([1,1].*iX,[-999,999],'-','color',[1,1,1].*0.7,'linewi',0.25); end; clear iX
    plot([0,0],[-1,1].*999,'k:','linewi',2)
     
    %handle y-axes and horizontal gridlines
    set(gca,'ytick',[])
    for iY=1:1:numel(Settings.Indices); plot([-1,1].*1000,[1,1].*iY,'-','color',Colours.(Settings.Indices{iY}),'linewi',0.5); end; clear iX

    %y-axis labels
    if iDir ~= 3; set(gca,'ytick',1:1:numel(Settings.Indices),'yticklabel',Settings.Indices); ylabel(SeasonNames{iSeason})
    else          set(gca,'ytick',1:1:numel(Settings.Indices),'yticklabel',{})
    end
    if iDir == 2; set(gca,'yaxislocation','right'); end

    %label R2
    text(min(get(gca,'xlim'))+0.99.*range(get(gca,'xlim')), ...
         numel(Settings.Indices),['R^2 = ',num2str(round(Reg.R2(iDir,iSeason,1),2))], ...
         'fontsize',11,'horizontalalignment','right')    

    for iIndex=1:1:numel(Settings.Indices)

      %get the value
      Value = Reg.Est(iDir,iSeason,iIndex);

      %get the standard error
      SE = Reg.SE(iDir,iSeason,iIndex);

      %if the value falls outside the range, scale it down until it fits
      Factor = 1;
      if Value-SE < min(get(gca,'xlim')) | Value+SE > max(get(gca,'xlim'))
        Valid  = 0;
        while Valid == 0
          Factor   = Factor+1;
          NewValue = Value./Factor;
          NewSE    =    SE./Factor;
          if NewValue-NewSE < min(get(gca,'xlim')) | NewValue+NewSE > max(get(gca,'xlim')); continue; end
          Valid = 1;
        end
        Value = NewValue;
        SE    = NewSE;
      end
% %         SE    = SE./NAOFactor;
% %         Value = Value./NAOFactor
% %         disp('scaling nao')
% %       end

      %is it statistically significant?
      Sig   = Reg.P(iDir,iSeason,iIndex);
      if Sig < 0.05; Colour = Colours.(Settings.Indices{iIndex}); else Colour = 'w'; end

      %indicate coefficient standard error
      StErr = Value + [-1,1]  .* SE;
      plot(StErr,[1,1].*iIndex,'k-','linewi',1); hold on
      plot([1,1].*StErr(1),[-1,1].*0.5+iIndex,'k-','linewi',1.5)
      plot([1,1].*StErr(2),[-1,1].*0.5+iIndex,'k-','linewi',1.5)    



      %plot data
      plot(Value,iIndex,'o', ...
           'color','k','markerfacecolor',Colour,'markersize',10)

      %if we scaled the data, label it
      if Factor ~= 1;
        if Value < 0; xpos = Value+SE+1; al = 'left';
        else          xpos = Value-SE-1; al = 'right';
        end

        text(xpos,iIndex,['x',num2str(Factor)],'fontsize',10,'horizontalalignment',al)
      end


      drawnow


    end; clear iIndex
  end; clear iSeason




end; clear iDir
