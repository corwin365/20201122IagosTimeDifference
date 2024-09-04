function analysis_routestats(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%meta information - statistics for route flight numbers
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/06/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('+++++++++++++++++++++++++++')
disp('Table of flights on each route')
disp('+++++++++++++++++++++++++++')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load([Settings.Paths.DataDir,'/',Settings.ID,'_routes.mat'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% identify all routes and the number of flights on them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = unique(RouteData.Dep);
b = unique(RouteData.Arr);

List = unique([a;b]);
clear a b

Grid = NaN(numel(List),numel(List));
for iX=1:1:numel(List)
  for iY=1:1:numel(List)
    % if iY >= iX; continue; end

    Pair = find(RouteData.Dep == List(iX) & RouteData.Arr == List(iY));
    if numel(Pair) == 0; continue; end
    if RouteData.NFlights(Pair) < Settings.Choices.MinFlights; continue; end
    Grid(iX,iY) = RouteData.NFlights(Pair);

  end
end
clear iX iY Pair


%drop empty rows/cols after above filters
a = find(nansum(Grid,1) > 0);
b = find(nansum(Grid,2) > 0);
idx = union(a,b);
Grid = Grid(:,idx);
Grid = Grid(idx,:);
List = List(idx);
clear a b idx 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% hence, make a table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix2latex(Grid,'test.tex',  ...
             'rowLabels',List, ...
             'columnLabels',List);

disp('Outputted city-pair flight numbers to test.tex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% also spit out a list of the most popular city pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('========================')
disp('Most popular city pairs:')
disp('========================')

[a,idx] = sort(Grid(:),'desc');
idx = idx(~isnan(Grid(idx)));
[ix,iy] = ind2sub(size(Grid),idx);

for iPair=1:1:10;
  pc = round((Grid(idx(iPair))./nansum(Grid(:)).*1000))./10;

  disp([List(ix(iPair))+" -> "+List(iy(iPair))+"  - "+num2str(Grid(idx(iPair)))+" flights ("+num2str(pc)+"%)"])
end

disp('========================')
disp('========================')

clear a idx ix iy iPair pc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot map of the European and North American airports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get location of the airports we're considering
GeoLoc = NaN(numel(List),2);
InfoArray = readtable([Settings.Paths.DataDir,'/airports.csv']);
InfoArray = InfoArray(:,[5,6,7,8,9,14,]); %lat,lon,elev,continent,country,iata code

for iAirport=1:1:numel(List)
  idx = find(contains(InfoArray.iata_code,List(iAirport)));
  GeoLoc(iAirport,1) = InfoArray.longitude_deg(idx);
  GeoLoc(iAirport,2) = InfoArray.latitude_deg(idx);
end
clear InfoArray iAirport idx



Lon = -120:0.01:30;
Lat =  0:0.01:90;
[xi,yi] = meshgrid(Lon,Lat);
Context = get_context(xi,yi,'SurfaceImage',true);

%%

%Europe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
clf

subplot(2,2,2)
m_proj('robinson','lat',[42,55],'lon',[-5,20]);
title('Europe')
m_image(Lon,Lat,Context.SurfaceImage);
hold on

for iAirport=1:1:numel(List)
  m_text(GeoLoc(iAirport,1), ...
         GeoLoc(iAirport,2), ...
         List(iAirport), ...
         'color','r','horizontalalignment','center', ...
         'fontsize',10,'fontweight','bold', ...
         'clipping','on')

end

 m_grid('xtick',[-185:5:180],'ytick',[-90:5:90])

%North America
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1)
m_proj('robinson','lat',[25,48],'lon',[-100,-68]);
title('North America')
m_image(Lon,Lat,Context.SurfaceImage);
hold on

for iAirport=1:1:numel(List)
  m_text(GeoLoc(iAirport,1), ...
         GeoLoc(iAirport,2), ...
         List(iAirport), ...
         'color','b','horizontalalignment','center', ...
         'fontsize',10,'fontweight','bold', ...
         'clipping','on')

end

 m_grid('xtick',[-185:10:180],'ytick',[-90:5:90])
