function analysis_cruiseheight(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%is there a correlation between climate indices and cruise heights?
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/10/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Cruise height assessment')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load flight metadata
Meta = load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat']);

%load indices
if Settings.Choices.ApplyLags == 0;
      Indices = load([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat'])
else; Indices = load([Settings.Paths.DataDir,'/',Settings.ID,'_laggedindices.mat']);
end

%load flight track wind data

%merge together
Tracks = load([Settings.Paths.DataDir,'/',Settings.ID,'_flighttracks.mat']);
Tracks = Tracks.Store;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each flight, generate a histogram of fraction of time at each level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute absolute histogram
ZScale = Settings.Choices.ZScale;
Store = NaN(numel(Tracks.Z),numel(ZScale));
for iFlight=1:1:numel(Tracks.Z)

  %extract track, and retain all heights above 8km (in order to avoid ascent and descent phases)
  t = Tracks.Z{iFlight};
  t(t < 9000) = NaN;
  Store(iFlight,:) = hist(t./1000,ZScale);
end; clear iFlight Tracks

%normalise to fraction of time...
NormHeight = Store;
Sigma = sum(NormHeight,2,'omitnan');
for iFlight=1:1:numel(Sigma);
  NormHeight(iFlight,:) = NormHeight(iFlight,:)./Sigma(iFlight);
end; clear iFlight
clear Sigma Store 


%then find a modal height for each flight, including round-trips
Meta.FlightData.ModalHeight = NaN(size(Meta.FlightData.FlightIndex));

%also, check what fraction of time the modal height represents
CheckFrac = Meta.FlightData.ModalHeight;

for iFlight=1:1:numel(Meta.FlightData.ModalHeight)
  if Meta.FlightData.Direction(iFlight) == 'R'
    Norm = NormHeight(Meta.FlightData.OriginalW(iFlight),:) + NormHeight(Meta.FlightData.OriginalE(iFlight),:);
    Scaler = 2;
  else
    Norm = NormHeight(Meta.FlightData.FlightIndex(iFlight),:);
    Scaler = 1;
  end

  [val,idx] = max(Norm,[],2);
  CheckFrac(iFlight) = val./Scaler;
  Meta.FlightData.ModalHeight(iFlight) = ZScale(idx);

end
clear iFlight Norm idx val Scaler

% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %% plot figure demonstrating validity of modal approach
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % % a) histogram of fraction of time at modal height
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % cla
% % % % % xlim([0 1])
% % % % % ylim([0 800])
% % % % % hold on; box on
% % % % % title('(c) Modal Proportion')
% % % % % xlabel('Fraction of flight')
% % % % % ylabel('Number of flights')
% % % % % set(gca,'xtick',0:0.1:1)
% % % % % 
% % % % % x = 0:1/100:1;
% % % % % h = hist(CheckFrac.*1,x);
% % % % % patch([x,100,0],[h,0,0],[255,153,204]./255,'edgecolor','none')
% % % % % 
% % % % % %overplot cumulative sum
% % % % % s = cumsum(h); s = s./max(s).*800;
% % % % % plot(x,800-s,'k','linewi',3)
% % % % % 
% % % % % for iY=0:200:800;
% % % % %   text(1.01,iY,[num2str(iY./800.*100),'%'],'fontsize',11)
% % % % % end
% % % % % 
% % % % % 
% % % % % %overplot key percentiles
% % % % % for iPC=5:5:95
% % % % %   x(closest(s./800,iPC./100))
% % % % %   plot([1,1].*x(closest(s./800,iPC./100)),[0,800],'k--')
% % % % % end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, regress this modal height against the climate indices
%this follows the same logic as analysis_regression.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%split out by season
%%%%%%%%%%%%%%%%%%%%

%move tables to top level to keep as similar syntax as possible
FlightData = Meta.FlightData; 
FlightIndices = Indices.FlightIndices; 
clear Meta Indices

%and split the delays into seasons
FlightData.SeasModHeight = FlightData.Delay; %copying this because I can't remember the syntax!
for iFlight=1:1:numel(FlightData.ModalHeight);
  FlightData.SeasModHeight(iFlight,:) = FlightData.ModalHeight(iFlight).*FlightData.InSeasons(iFlight,:);
end; clear iFlight


%%allocate storage
%%%%%%%%%%%%%%%%%%%%

%create storage matrices
Reg = struct;
Reg.Est = NaN(numel(Settings.Choices.Directions), ...
              numel(Settings.Seasons.List),       ...
              numel(Settings.Indices.List)); 
Fields = {'SE','T','P','N','R2'}; %these are the output fields we'll retain from the regression
for iF=1:1:numel(Fields); Reg.(Fields{iF}) = Reg.Est; end
clear Fields iF

%%process
%%%%%%%%%%%%



%for each season and direction...
for iSeason=1:1:numel(Settings.Seasons.List)
  for iDirection=1:1:numel(Settings.Choices.Directions)

    %join indices onto flight data
    Data = innerjoin(FlightData,FlightIndices.(Settings.Choices.Directions{iDirection}));

    %find flights in this season and direction
    InThisSeason    = find(table2array(Data.InSeasons(:,iSeason)) == 1);
    InThisDirection = find(Data.Direction == Settings.Choices.Directions{iDirection});
    InThisSet       = intersect(InThisDirection,InThisSeason);

    %hence, cut down the master table to just these entries
    FlightsInUse = Data(InThisSet,:);

    %and to just the columns needed for the regression
    ModelInput = table(splitvars(table(FlightsInUse.SeasModHeight)).(Settings.Seasons.List{iSeason})); 
    ModelInput = renamevars(ModelInput,'Var1','SeasModHeight');
    for iVar=1:1:numel(Settings.Indices.List)
      ModelInput = addvars(ModelInput,FlightsInUse.(Settings.Indices.List{iVar}), ...
                           'NewVariableNames',Settings.Indices.List{iVar});
    end; clear iVar

    %fit the linear model
    mdl = fitlm(ModelInput,'ResponseVar','SeasModHeight');

    %store what we want
    Coefs = table2array(mdl.Coefficients);
    Reg.Est(iDirection,iSeason,:) = Coefs(2:end,1);        % Coefficient estimates for each corresponding term in the model.
    Reg.SE( iDirection,iSeason,:) = Coefs(2:end,2);        % Standard error of the coefficients.
    Reg.T(  iDirection,iSeason,:) = Coefs(2:end,3);        % t-statistic for each coefficient to test the null hypothesis that the corresponding coefficient is zero against the alternative that it is different from zero, given the other predictors in the model. Note that tStat = Estimate/SE.
    Reg.P(  iDirection,iSeason,:) = Coefs(2:end,4);        % p-value for the t-statistic of the hypothesis test that the corresponding coefficient is equal to zero or not.
    Reg.N(  iDirection,iSeason,:) = mdl.NumObservations;   % Ronseal
    Reg.R2( iDirection,iSeason,:) = mdl.Rsquared.Adjusted; % adjusted coefficient of determination    

    %tidy up
    clear FlightsInUse ModelInput mdl Coefs Data

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
subplot = @(m,n,p) subtightplot (m, n, p, 0.04, 0.1, 0.1);
Letters = 'abcdefghijklmnpqrstuvwxyz';

XLimit = 670; %m

for iSeason=1:1:numel(Settings.Seasons.List)
  for iDirection=1:1:numel(Settings.Choices.Directions)

    %prepare panel
    k = k+1;
    subplot(numel(Settings.Seasons.List),numel(Settings.Choices.Directions),k)
    cla
    hold on; box on; grid on;
    axis([-XLimit,XLimit, 0.5, numel(Settings.Indices.List)+0.5])
    set(gca,'ydir','reverse')

    %handle x axis and vertical gridlines
    if     iSeason == 1; set(gca,'xaxislocation','top');
    else                 set(gca,'xaxislocation','bottom'); 
    end
    if iSeason ~= 1 & iSeason ~= numel(Settings.Seasons.List); set(gca,'xtick',[]); 
    else; xlabel('Modal Height Change [m]');
    end
    plot([1,1].*min(get(gca,'xlim')),minmax(get(gca,'ylim')),'w-','linewi',2)
    plot([1,1].*max(get(gca,'xlim')),minmax(get(gca,'ylim')),'w-','linewi',2)    
    for iX=-1000:100:1000; plot([1,1].*iX,[-999,999],'-','color',[1,1,1].*0.7,'linewi',0.25); end; clear iX
    plot([0,0],[-1,1].*999,'k:','linewi',2)

    %handle y-axes and horizontal gridlines
    set(gca,'ytick',[])
    for iY=1:1:numel(Settings.Indices.List);
      x = [-1,1,1,-1].*1000;
      y = ([-1,-1,1,1].*0.15)+iY;
      % patch(x,y,Settings.Indices.Colours.(Settings.Indices.List{iY}),'edgecolor','none','facealpha',0.55)
      % plot([-1,1].*1000,[1,1].*iY,'-','color',Settings.Indices.Colours.(Settings.Indices.List{iY}),'linewi',1,'linestyle',':');
      plot([-1,1].*1000,[1,1].*iY,'-','color',[1,1,1].*0.6,'linewi',1,'linestyle',':');    end; clear iY




    %y-axis labels
    if iDirection == 1 | iDirection == numel(Settings.Choices.Directions)
      set(gca,'ytick',1:1:numel(Settings.Indices.List),'yticklabel',Settings.Indices.List);
    end
    if iDirection == numel(Settings.Choices.Directions); set(gca,'yaxislocation','right'); end

    %plot indices
    for iIndex=1:1:numel(Settings.Indices.List)

      %get value, error and significance
      Value = Reg.Est(iDirection,iSeason,iIndex).*1000; %minutes
      SE    = Reg.SE( iDirection,iSeason,iIndex).*1000;
      p     = Reg.P(  iDirection,iSeason,iIndex);

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
      plot([1,1].*StErr(1),[-1,1].*0.5+iIndex,'k-','linewi',0.5)
      plot([1,1].*StErr(2),[-1,1].*0.5+iIndex,'k-','linewi',0.5)     

      if p < 0.05; LW = 2; Colour = Settings.Indices.Colours.(Settings.Indices.List{iIndex}); 
      else         LW = 1; Colour = 'w';
      end

      % plot(Value,iIndex,Settings.Indices.Symbols.(Settings.Indices.List{iIndex}),...
      %      'color', 'k', 'linewi',LW, ...
      %      'markerfacecolor',Colour, ...
      %      'markersize',     10)
      
      plot(Value,iIndex,'o',...
           'color', 'k', 'linewi',LW, ...
           'markerfacecolor',Colour, ...
           'markersize',     10)      

    end; clear iIndex



    %label R2
    text(min(get(gca,'xlim'))+0.995.*range(get(gca,'xlim')), ...
         numel(Settings.Indices.List)+0.37,['R^2 = ',num2str(round(Reg.R2(iDirection,iSeason,1),2))], ...
         'fontsize',11,'horizontalalignment','right','verticalalignment','bottom','FontWeight','bold')    

     %label panel
     text(min(get(gca,'xlim')),0.5,['(',num2str(Letters(k)),')'],'horizontalalignment','left','VerticalAlignment','top','Clipping','off')

    if iDirection == 1 
      ylabel(Settings.Seasons.List{iSeason},'fontsize',30,'fontweight','bold')
    end

  end; clear iDirection
end; clear iSeason
