function [] = gg_regression(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%multilinearly regress relative time taken against climate indices
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/12/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load route and flight data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlightData = load('data/flight_data.mat');
RouteData  = load('data/routes.mat');
Indices    = load('data/indices.mat');

%number of airports
NPorts = numel(RouteData.Airports);

%list of seasons
Seasons = fieldnames(Settings.Seasons);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute relative flight times, and split into Eward and Wward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlightData.Results.tRel  = NaN.*FlightData.Results.t;
FlightData.Results.Eward = NaN.*FlightData.Results.t;

for iDep=1:1:NPorts
  for iArr=1:1:NPorts
   for iSeason=1:1:numel(Seasons)
    
     %find all flights meeting these criteria
     Flights = squeeze(RouteData.Flights(iSeason,iDep,iArr,:));
     Flights = Flights(~isnan(Flights));
     if numel(Flights) == 0; continue; end

     %what is the relative time taken by each flights?
     FlightData.Results.tRel(Flights) = FlightData.Results.t(Flights)./nanmedian(FlightData.Results.t(Flights));
     
     %are we heading east or west?
     if ismember(RouteData.Airports{iDep},Settings.NA); E = 1; else E = 2; end  %1 eastward, 2 westward
     FlightData.Results.Eward(Flights) = E;
     
   end
  end
end
clear iDep iArr iSeason Flights NPorts E

%also compute the mean overall flight time, so we can convert our
%coefficients from relative time to minutes
MedianFlightTime = nanmedian(FlightData.Results.t(:))./60; %MINUTES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate lag analysis data
% if lagging is switched off in the master settings
%the same code is used but with a 'lag' of zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%prepare lagging matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%

%the number of indices can vary, so we need to condense this into
%an array which can be done in a single loop
%do this by creating an array where:
  %each column is an index
  %each row is a different set of lags

if Settings.Reg.Lag == 1;
  
  %generate the list of lag-steps to be used for each index
  LagScale =Settings.Reg.Steps;
  NLags    = numel(LagScale);
  NIndices = numel(Settings.Indices);
  LagMatrix = NaN(NLags.^NIndices,NIndices);
  disp(['Testing ',num2str(NLags.^NIndices),' combinations, this may take some time'])

  %produce a lag matrix. 
  %let's take advantage of base-conversion algebra to make all the 
  %possible numeric combinations. This is *very* cheaty, and took me
  %hours to think up. It's fast though and only a few lines
  if NIndices ~= 1;
    Base10      = (0:1:NLags.^NIndices-1)';
    BaseIndices = cellstr(dec2base(Base10,NLags,NIndices));
    textprogressbar('Generating combinations ')
    for iRow=1:1:size(LagMatrix,1)
      String = BaseIndices{iRow};
      for iCol=1:1:NIndices;
        LagMatrix(iRow,iCol) = base2dec(String(iCol),NLags);
      end
      if mod(iRow,100);textprogressbar(iRow./size(LagMatrix,1).*100); end
    end
    LagMatrix= LagScale(LagMatrix+1);
    clear LagScale NLags NIndices Base10 BaseIndices iRow iCol String
    textprogressbar(100); textprogressbar('!')
  else
    LagMatrix(:) = LagScale;
  end
else
  LagMatrix = zeros(1,numel(Settings.Indices));
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regress for each season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



textprogressbar('Regressing data against indices ')
for iSeason=1:1:numel(Seasons)

  for iEast=1:2; %1 eastward, 2 westward
  
    %get all the flights in this season
    Flights = flatten(RouteData.Flights(iSeason,:,:,:));
    Flights = Flights(~isnan(Flights));
    
    %and going in the right direction
    Dir = FlightData.Results.Eward(Flights);
    Flights = Flights(Dir == iEast);    
    
    %get their relative times
    tRel = FlightData.Results.tRel(Flights);   
    

    
    %loop over possible lags, optimising for the largest R2
    R2 = 0;
    
    for iLag=1:1:size(LagMatrix,1)
      
      %get the necessary climate indices for this lag
      Beta = NaN(numel(Settings.Indices),numel(tRel));
      for iIndex=1:1:numel(Settings.Indices)
        if LagMatrix(iLag,iIndex) == 0;
          Index = Indices.Indices.(Settings.Indices{iIndex});
          Beta(iIndex,:) = Index(Flights);
          Lags(iIndex) = 0;
        else
          Index = Indices.Indices.Lagged.(Settings.Indices{iIndex});
          idx = closest(Indices.Indices.Lagged.LaggedScale,LagMatrix(iLag,iIndex));
          Beta(iIndex,:) = Index(idx,Flights);
          Lags(iIndex) = Indices.Indices.Lagged.LaggedScale(idx); %just in case the arrays don't match, this will store the actual value used
        end
      end
      Beta = Beta';
      clear iIndex Index idx
      
      %the Beta matrix is the regression indices
      %the tRel vector is the resulting times
    
    
      %fit a linear rregression model
      mdl = fitlm(Beta,tRel);
      
      %if R2 is better than any previous attempt, store
      if mdl.Rsquared.Adjusted > R2;
        
        %convert the coefficients to a matrix
        Coefs = table2array(mdl.Coefficients);
        
        
        %and store them
        Reg.Est(iEast,iSeason,:) = Coefs(2:end,1);
        Reg.SE( iEast,iSeason,:) = Coefs(2:end,2);
        Reg.T(  iEast,iSeason,:) = Coefs(2:end,3);
        Reg.P(  iEast,iSeason,:) = Coefs(2:end,4);
        
        %also store the number of flights used
        Reg.N(  iEast,iSeason,:) = mdl.NumObservations;
        
        %and the R2
        Reg.R2( iEast,iSeason,:) = mdl.Rsquared.Adjusted;
        R2 = mdl.Rsquared.Adjusted;
        
        %and the lag
        Reg.Lag(iEast,iSeason,:) = Lags;
     
      end
      
    end; clear iLag BetaLag iIndex Lags R2
    
    
  end
  textprogressbar(iSeason./numel(Seasons).*100)
end
textprogressbar('!')
clear iSeason Flights tRel mdl Coefs Beta RouteData Indices FlightData iEast Dir LagMatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.02,0.02], [0.12,0.09], 0.15);

%assign colours and symbols to the regression coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Symbols = 'sd<o^>v';
Colours = [[ 57,159,228]; ...
           [113, 69,168]; ...
           [196, 66, 79]; ...
           [ 46,148,130]; ...
           [152, 51, 91]; ...
           [255,209,107]; ...
           [ 69,174, 98]]./255;
Letters = 'acdfgijlmop';         
         
         
%plot individual panels         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



k = 0;
for iSeason=1:1:numel(Seasons)
  for iEast=1:2  %1 eastward, 2 westward
  
    %create panel
    k = k+1;
    subplot(numel(Seasons),2,k)
    
    
    %set axis lengths
    if iEast == 1; XRange = [-32,22]; TextLoc = -31;
    else           XRange = [-20,44]; TextLoc = 12;
    end
    YRange = [0.5,size(Reg.Est,3)+1.5];
    
    
    %patch backgrounds
    for iIndex=1:1:numel(Settings.Indices)
      NewColour = colorGradient(Colours(iIndex,:),[1,1,1],12);
      NewColour = NewColour(10,:);
      
      
      patch(XRange([1,2,2,1,1]),-1.5+iIndex+[1,1,2,2,1], ...
        NewColour,'edgecolor','none');%[1,1,1].*0.7)
      hold on
    end
    for x = -1000:10:1000; 
      plot([1,1].*x,[0,2+numel(Settings.Indices)],'k:','linewi',0.5)
    end
    set(gca','xtick',-1000:10:1000)
    
    %prepare axes
    if     iSeason ==              1; set(gca,'xaxislocation','top');
    elseif iSeason == numel(Seasons); set(gca,'xaxislocation','bottom'); xlabel('Delay [minutes]')
    else                              set(gca,'xticklabel',[]);
    end
    box on; grid off
    set(gca,'tickdir','out','ydir','reverse')
      
    
    %plot data
    for iIndex=1:1:numel(Settings.Indices)
            
      
      
      %get the value
      Value = Reg.Est(iEast,iSeason,iIndex).*MedianFlightTime.*2; %the 2 is due to how we normalised the indices back in routine "ff"
      
      %is it statistically significant?
      Sig   = Reg.P(iEast,iSeason,iIndex);
      if Sig < 0.05; Colour = Colours(iIndex,:); else Colour = 'w'; end      
      
      %indicate coefficient standard error
      StErr = Value + [-1,1] .* MedianFlightTime .* Reg.SE(iEast,iSeason,iIndex).*2;
      plot(StErr,[1,1].*iIndex,'k-','linewi',1); hold on
      plot([1,1].*StErr(1),[-1,1].*0.5+iIndex,'k-')
      plot([1,1].*StErr(2),[-1,1].*0.5+iIndex,'k-')      
      
      h = plot(Value,iIndex, ...
               'marker',Symbols(iIndex), ...
               'color','k','markerfacecolor',Colour, ...
               'markersize',10);     
      hold on
%       setMarkerColor(h,Colours(iIndex,:),Alpha);      
      
      
      if Settings.Reg.Lag == 1
        Shift = Reg.Lag(iEast,iSeason,iIndex);
        if Shift > 0; Shift = ['+',num2str(Shift)]; else Shift = num2str(Shift); end
        text(Value,iIndex+0.5,Shift, ...
          'horizontalalignment','center','fontsize',10)
      end
      
    end
    
    
  %tidy
  if iSeason == 1;
    if iEast == 1; title('Eastward','fontsize',24);
    else;          title('Westward','fontsize',24);
    end
  end
  if iEast ~= 1; set(gca,'yaxislocation','right'); end
  
  set(gca,'ytick',[])
  ylabel(Seasons{iSeason},'fontsize',24);
  plot([0,0],[0,size(Reg.Est,3)+2],'-','linewi',2,'color',[1,1,1].*0.5)
  xlim(XRange)
  ylim(YRange)
  text(TextLoc,size(Reg.Est,3)+0.8,['(',Letters(k),') N=',num2str(Reg.N(iEast,iSeason,1)),'; R^2=',num2str(0.01.*round(100.*Reg.R2( iEast,iSeason,1)))],'fontsize',12)
  
  drawnow    
    
  end



end


% %key - relative to last plot
% for iCoef=1:1:size(Reg.Est,3)
%   
%   x = -80+12.*iCoef;
%   plot(x,11,'clipping','off','markersize',20, ...
%        'marker',Symbols(iCoef),'color','k','markerfacecolor',Colours(iCoef,:))
%   text(x,13,Settings.Indices{iCoef},'fontsize',12,'horizontalalignment','center')  
%   
% end
% 
