function analysis_regression_lagged(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%KDF analysis of merged aircraft data
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/01/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Multilinear regression analysis, lagged')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])

%load indices
load([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat'])

%join indices onto flight data, and drop
Data = innerjoin(FlightData,FlightIndices);

clear FlightData FlightIndices RangeStore RouteData OverallMedianTimes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% allocate storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create storage matrices
Reg = struct;
Reg.Est = NaN(numel(Settings.Choices.Directions), ...
              numel(Settings.Seasons.List),       ...
              numel(Settings.Indices.List)); 
Fields = {'SE','T','P','N','R2','dt'}; %these are the output fields we'll retain from the regression
for iF=1:1:numel(Fields); Reg.(Fields{iF}) = Reg.Est; end
clear Fields iF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate a load of lag combinations
IndexNames = DateIndices.Properties.VariableNames(2:end);
LagScale = Settings.Choices.LaggedRegression.LagScale;
C = permn(LagScale,numel(IndexNames));

%combinatorics means this is big. However, we can reduce it as:
%1. it's not meaningful to lag our linear time ramp 
idx = find(strcmp(IndexNames,'Time'));
good = find(C(:,idx) == 0);
C = C(good,:);
clear idx good
%2. tests show the NAO fit is almost always best at near-zero lag, which makes sense as it's geographically local
idx = find(strcmp(IndexNames,'NAO'));
good = find(C(:,idx) == 0);
C = C(good,:);
clear idx good



%ok. Now for each season and direction...
for iSeason=1:1:numel(Settings.Seasons.List)
  for iDirection=1:1:numel(Settings.Choices.Directions)
    textprogressbar([Settings.Seasons.List{iSeason},' ',Settings.Choices.Directions{iDirection},' '])
    for iLag=1:1:size(C,1)
      if mod(iLag,100) == 1;textprogressbar(iLag./size(C,1).*100); end

      %find flights in this season and direction
      InThisSeason    = find(table2array(Data.InSeasons(:,iSeason)) == 1);
      InThisDirection = find(Data.Direction == Settings.Choices.Directions{iDirection});
      InThisSet = intersect(InThisDirection,InThisSeason);

      %hence, cut down the master table to just these entries
      FlightsInUse = Data(InThisSet,:);
      clear InThisSeason InThisDirection InThisSet

      %lag the data. This involves reinterpolating each index to the observed time
      Lags = C(iLag,:);
      Indices = DateIndices;
      for iIndex=2:1:size(Indices,2)
        if Lags(iIndex-1) ~= 0;
          Indices(:,iIndex) = circshift(Indices(:,iIndex),-Lags(iIndex-1),1);
          % % for j=1:1:abs(Lags(iIndex-1))
            % % Indices(      j,iIndex) = array2table(NaN);
            % % Indices(end-j+1,iIndex) = array2table(NaN);
          % % end
        end
      end
      for iIndex=1:1:numel(IndexNames); 
        FlightsInUse.(IndexNames{iIndex}) = interp1(Indices.Date,Indices.(IndexNames{iIndex}),FlightsInUse.Date); 
      end
      clear Lag Indices iIndex



      %select just the columns needed for the regression
      ModelInput = table(splitvars(table(FlightsInUse.Delay)).(Settings.Seasons.List{iSeason}));
      ModelInput = renamevars(ModelInput,'Var1','Delay');
      for iVar=1:1:numel(Settings.Indices.List)
        ModelInput = addvars(ModelInput,FlightsInUse.(Settings.Indices.List{iVar}), ...
                            'NewVariableNames',Settings.Indices.List{iVar});
      end; clear iVar


      %fit the linear model
      mdl = fitlm(ModelInput,'ResponseVar','Delay');

      %store what we want
      Coefs = table2array(mdl.Coefficients);

      if iLag ~=1 & mdl.Rsquared.Adjusted < Reg.R2( iDirection,iSeason,1) ; continue; end %this is a worse fit than last time, so don't store it

      Coefs = table2array(mdl.Coefficients);
      Reg.Est(iDirection,iSeason,:) = Coefs(2:end,1);        % Coefficient estimates for each corresponding term in the model.
      Reg.SE( iDirection,iSeason,:) = Coefs(2:end,2);        % Standard error of the coefficients.
      Reg.T(  iDirection,iSeason,:) = Coefs(2:end,3);        % t-statistic for each coefficient to test the null hypothesis that the corresponding coefficient is zero against the alternative that it is different from zero, given the other predictors in the model. Note that tStat = Estimate/SE.
      Reg.P(  iDirection,iSeason,:) = Coefs(2:end,4);        % p-value for the t-statistic of the hypothesis test that the corresponding coefficient is equal to zero or not.
      Reg.N(  iDirection,iSeason,:) = mdl.NumObservations;   % Ronseal
      Reg.R2( iDirection,iSeason,:) = mdl.Rsquared.Adjusted; % adjusted coefficient of determination
      Reg.dt( iDirection,iSeason,:) = C(iLag,:);



      %tidy up
      clear FlightsInUse ModelInput mdl Coefs

    end; clear iLag
    textprogressbar(100); textprogressbar('!')

  end; clear iDirection
end; clear iSeason

%multiply values by 2 (as we defined the indices on a -1 to +1 range)
Reg.Est = Reg.Est .* 2;
Reg.SE  = Reg.SE  .* 2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure
clf
k = 0;
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.03,0.03], 0.1, 0.1);

XLimit = 15; %minutes

for iSeason=1:1:numel(Settings.Seasons.List)
  for iDirection=1:1:numel(Settings.Choices.Directions)

    %prepare panel
    k = k+1;
    subplot(numel(Settings.Seasons.List),numel(Settings.Choices.Directions),k)
    cla
    hold on; box on; grid off;
    axis([-XLimit,XLimit, 0.1, numel(Settings.Indices.List)+0.5])
    set(gca,'ydir','reverse')

    %handle x axis and vertical gridlines
    if     iSeason == 1; set(gca,'xaxislocation','top');
    else                 set(gca,'xaxislocation','bottom'); 
    end
    if iSeason ~= 1 & iSeason ~= numel(Settings.Seasons.List); set(gca,'xtick',[]); 
    else; xlabel('Delay [minutes]');
    end
    plot([1,1].*min(get(gca,'xlim')),minmax(get(gca,'ylim')),'w-','linewi',2)
    plot([1,1].*max(get(gca,'xlim')),minmax(get(gca,'ylim')),'w-','linewi',2)    
    for iX=-100:5:100; plot([1,1].*iX,[-999,999],'-','color',[1,1,1].*0.7,'linewi',0.25); end; clear iX
    plot([0,0],[-1,1].*999,'k:','linewi',2)

    %handle y-axes and horizontal gridlines
    set(gca,'ytick',[])
    for iY=1:1:numel(Settings.Indices.List); 
      plot([-1,1].*1000,[1,1].*iY,'-','color',Settings.Indices.Colours.(Settings.Indices.List{iY}),'linewi',2,'linestyle',':'); 
    end; clear iY

    %y-axis labels
    if iDirection == 1 | iDirection == numel(Settings.Choices.Directions)
      set(gca,'ytick',1:1:numel(Settings.Indices.List),'yticklabel',Settings.Indices.List);
    end
    if iDirection == numel(Settings.Choices.Directions); set(gca,'yaxislocation','right'); end

    %plot indices
    for iIndex=1:1:numel(Settings.Indices.List)

      %get value, error and significance
      Value = Reg.Est(iDirection,iSeason,iIndex)./60; %minutes
      SE    = Reg.SE( iDirection,iSeason,iIndex)./60;

      %if it's outside the plot range, scale until it is
      if abs(Value) > XLimit; 
        Sign = Value./abs(Value);
          Value2 = Sign.*XLimit;        
        if Sign > 0; 
          plot(Value2+0.8,iIndex,'>','clipping','off','markerfacecolor',Settings.Indices.Colours.(Settings.Indices.List{iIndex}),'color','k')
          text(Value2+1.2,iIndex,[num2str(round(Value))],'HorizontalAlignment', 'left','VerticalAlignment','middle','FontSize',11,'color',Settings.Indices.Colours.(Settings.Indices.List{iIndex}),'fontweight','bold')
        else        
          plot(Value2-0.8,iIndex,'<','clipping','off','markerfacecolor',Settings.Indices.Colours.(Settings.Indices.List{iIndex}),'color','k')
          text(Value2-1.2,iIndex,[num2str(round(Value))],'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',11,'color',Settings.Indices.Colours.(Settings.Indices.List{iIndex}),'fontweight','bold')
        end
        Value = Value2; clear Value2
      end

      %plot
      StErr = Value + [-1,1]  .* SE;
      plot(StErr,[1,1].*iIndex,'k-','linewi',1)  
      plot([1,1].*StErr(1),[-1,1].*0.5+iIndex,'k-','linewi',1)
      plot([1,1].*StErr(2),[-1,1].*0.5+iIndex,'k-','linewi',1)     

      plot(Value,iIndex,'o',...
           'color', 'k', 'linewi',1, ...
           'markerfacecolor',Settings.Indices.Colours.(Settings.Indices.List{iIndex}), ...
           'markersize',     15)

       %print lag
       text(Value,iIndex-0.5,num2str(Reg.dt( iDirection,iSeason,iIndex)),'fontsize',10,'HorizontalAlignment','center')


    end; clear iIndex



    %label R2
    text(min(get(gca,'xlim'))+0.995.*range(get(gca,'xlim')), ...
         numel(Settings.Indices.List)+0.25,['R^2 = ',num2str(round(Reg.R2(iDirection,iSeason,1),2))], ...
         'fontsize',11,'horizontalalignment','right')    

    if iDirection == 1 
      ylabel(Settings.Seasons.List{iSeason},'fontsize',48,'fontweight','bold')
    end

  end; clear iDirection
end; clear iSeason


function combs = nmultichoosek(values, k)
%// Return number of multisubsets or actual multisubsets.
if numel(values)==1 
    n = values;
    combs = nchoosek(n+k-1,k);
else
    n = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
end
