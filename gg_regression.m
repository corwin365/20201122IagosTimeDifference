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
%% extract valid points, and regress for each season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    %get the corresponding climate indices
    Beta = NaN(numel(Settings.Indices),numel(tRel));
    for iIndex=1:1:numel(Settings.Indices)
      Index = Indices.Indices.(Settings.Indices{iIndex});
      Beta(iIndex,:) = Index(Flights);
    end
    clear iIndex Index
    
    %the Beta matrix is the regression indices
    %the tRel vector is the resulting times
    %fit a linear rregression model
    mdl = fitlm(Beta',tRel);
    
    %convert the coefficients to a matrix
    Coefs = table2array(mdl.Coefficients);
    
    
    %and store them
    Reg.Est(iEast,iSeason,:) = Coefs(2:end,1);
    Reg.SE( iEast,iSeason,:) = Coefs(2:end,2);
    Reg.T(  iEast,iSeason,:) = Coefs(2:end,3);
    Reg.P(  iEast,iSeason,:) = Coefs(2:end,4);
    
    %also store the number of flights used
    Reg.N(  iEast,iSeason,:) = mdl.NumObservations;
    
  end
end
clear iSeason Flights tRel mdl Coefs Beta RouteData Indices FlightData iEast Dir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create figure
figure
clf
set(gcf,'color','w')
% subplot = @(m,n,p) subtightplot (m, n, p, 0.08, 0.05, 0.05);


Symbols = 'sd<o^>v';
Colours = 'rgbycmk';

k = 0;
for iSeason=1:1:numel(Seasons)
  for iEast=1:2  %1 eastward, 2 westward

  
    %create panel
    k = k+1;
    subplot(numel(Seasons),2,k)
    
    %plot data
    for iIndex=1:1:numel(Settings.Indices)
      
      %get the value
      Value = Reg.Est(iEast,iSeason,iIndex).*MedianFlightTime;
      
      %is it statistically significant?
      Sig   = Reg.P(iEast,iSeason,iIndex);
      if Sig < 0.05; Alpha = 0.9; else Alpha = 0; end      
      
      %indicate coefficient standard error
      StErr = Value + [-1,1] .* MedianFlightTime .* Reg.SE(iEast,iSeason,iIndex);
      plot(StErr,[1,1].*iIndex,'k-','linewi',1); hold on
      plot([1,1].*StErr(1),[-1,1].*0.5+iIndex,'k-')
      plot([1,1].*StErr(2),[-1,1].*0.5+iIndex,'k-')      
      
      h = plot(Value,iIndex, ...
               'marker',Symbols(iIndex), ...
               'color','k','markerfacecolor',Colours(iIndex), ...
               'markersize',10);     
      hold on
      setMarkerColor(h,Colours(iIndex),Alpha);      
      
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
  plot([0,0],[0,size(Reg.Est,3)+2],'k--','linewi',1)
  xlim([-30,+40])
  ylim([0,size(Reg.Est,3)+1])
  text(-28,size(Reg.Est,3),['N=',num2str(Reg.N(iEast,iSeason,1))]);
  
  drawnow    
    
  end



end


%key - relative to last plot
for iCoef=1:1:size(Reg.Est,3)
  
  x = -60+12.*iCoef;
  plot(x,-3,'clipping','off','markersize',20, ...
       'marker',Symbols(iCoef),'color','k','markerfacecolor',Colours(iCoef))
  text(x,-4.6,Settings.Indices{iCoef},'fontsize',12,'horizontalalignment','center')  
  
end


