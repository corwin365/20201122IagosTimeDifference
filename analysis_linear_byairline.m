function analysis_linear_byairline(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%linear trend analysis of the data, subdividing the data by plane
%
%it turned out that there wasn't enough data from non-Lufthansa flights
%to draw any conclusions about whether different airlines behaved differently,
%which I was trying to use as a cross-check on the overall trend to see if
%it was an operational choice rather than a secular trend
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/09/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Linear trend plot')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load indices
load([Settings.Paths.DataDir,'/',Settings.ID,'_indices.mat'])

%load data
load([Settings.Paths.DataDir,'/',Settings.ID,'_flightinfo_normalised.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% split data by plane and compute linear trends
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find all unique planes
Planes = unique(FlightData.Plane(~ismissing(FlightData.Plane)));

idx = [];

% %discard any with less than 5% of the flights
% Total = numel(FlightData.Plane(~ismissing(FlightData.Plane)));
% for iPlane=1:1:numel(Planes)
% 
%   N = numel(find(FlightData.Plane == Planes(iPlane)));
%   if N > Total/20; idx(end+1) = iPlane; end
% 
% end
% Planes = Planes(idx);




clf
k = 0;
subplot = @(m,n,p) subtightplot (m, n, p, 0.01, 0.03, 0.03);

for iDirection=[1,3];%1:1:numel(Settings.Choices.Directions)
  for iPlane=1:1:numel(Planes)

    %find flights by this plane in this direction
    ByThisPlane   = find(FlightData.Plane == Planes(iPlane));
    InThisDirection = find(FlightData.Direction == Settings.Choices.Directions{iDirection});
    InThisSet = intersect(InThisDirection,ByThisPlane);
    FlightsInUse = FlightData(InThisSet,:);

    %get dates
    Dates = FlightsInUse.Date;

    %convert to minutes of delay
    t = table2array(FlightsInUse.Delay)./60;

    %for each month, show the spread of data
    Months= datenum(1994,8:1:(40*12),1); %too long, but we'll skip empty months
    
    PC = [18:1:82];
    Store = NaN(numel(Months),numel(PC));

    for iMonth=1:1:numel(Months)-1
      InThisMonth = find(Dates >= Months(iMonth) & Dates < Months(iMonth+1));
      if numel(InThisMonth) == 0; continue; end

      %remove data after august 2023, so that we start and end on the same month
      if Months(iMonth) > datenum(2023,8,1); continue; end

      Store(iMonth,:) = prctile(t(InThisMonth),PC);
    end



      %plot data
      k = k+1;
      subplot(2,numel(Planes),k)
      hold on
      title(Planes(iPlane))
      % 
      % %filled shading version
      % for iPC=1:1:numel(PC)/2
      %   x = [Months,Months(end:-1:1)];
      %   y = [Store(:,iPC);Store(end:-1:1,end-iPC+1)];
      %   Good = find(~isnan(x+y'));
      %   x = x(Good); y = y(Good);
      %   patch(x,y,[1,1,1].*0.5,'edgecolor','none','facealpha',0.02)
      % end


      %monthly line version
      x = Months;
      y = Store(:,round(numel(PC)./2));
      for iMonth=1:1:numel(x)
        plot(x(iMonth).*[1,1], ...
             Store(iMonth,[1,end]),'-','color',[1,1,1].*0.7)
      end
      plot(x,y,'ko')


      %tidy up
      datetick
      xlim([datenum(1994,8,1),datenum(2024,4,1)])
      % if iDirection == 2;  
      %   % set(gca,'yaxislocation','right')
      % end
      % ylim(minmax(y)+[-1,1].*10)
      % hold on
      % 
      % plot([datenum(1994,8,1),datenum(2024,4,1)],[0,0],'k:','linewi',0.5)
      % box on
      % grid on
      % axis on
      % 
      % if iDirection == 1;     set(gca,'xaxislocation','top');
      % elseif iDirection == 2; set(gca,'xticklabel',{});
      % end
      % 
      % ylabel('Delay [minutes]')


    % % % %compute and overplot SEASONAL trend
    % % % subplot(numel(Settings.Choices.Directions),1,iDirection)
    % % % 
    % % % switch Settings.Seasons.List{iSeason}
    % % %   case 'All';
        LineColour = 'k';
        LineStyle  = '-';
        LineWidth  = 3;
    % % %   case 'DJF'; 
    % % %     LineColour = [255,102,178]./255;
    % % %     LineStyle  = '-';
    % % %     LineWidth  = 2;   
    % % %   case 'MAM'; 
    % % %     LineColour = [255,128,0]./255;
    % % %     LineStyle  = '-';
    % % %     LineWidth  = 2;             
    % % %   case 'JJA'; 
    % % %     LineColour = [152, 51, 91]./255;
    % % %     LineStyle  = '-';
    % % %     LineWidth  = 2;           
    % % %   case 'SON'; 
    % % %     LineColour = [ 69,174, 98]./255;
    % % %     LineStyle  = '-';
    % % %     LineWidth  = 2;             
    % % %   otherwise
    % % %     LineColour = [1,1,1].*0.5;
    % % %     LineStyle  = '-';
    % % %     LineWidth  = 2; 
    % % % end
    % % % 
    % % % 
    MedianSeries = Store(:,round(numel(PC)./2));
    Good = find(~isnan(MedianSeries));

    %linear fit for season, and label
    [p,S] = polyfit(Months(Good),MedianSeries(Good),1);
    y = polyval(p,Months(Good));
    plot(Months(Good),y,'linewi',LineWidth,'color',LineColour,'LineStyle',LineStyle)
    % text(datenum(1989+6.*iSeason,6,1),min(gca().YLim)+0.02.*range(gca().YLim),[num2str(round(p(1)*365*10,2,'significant')),' min/dec'], ...
    %      'VerticalAlignment','bottom','color',LineColour,'FontSize',14)
    % % % 
    % % % 
    % % % %monthnames
    % % % if iDirection == numel(Settings.Choices.Directions)
    % % %   text(datenum(1993+5.*iSeason,1,1),min(gca().YLim)-0.33.*range(gca().YLim),Settings.Seasons.List{iSeason}, ...
    % % %      'VerticalAlignment','bottom','color',LineColour,'FontSize',20)
    % % % end

    drawnow
  end
end
