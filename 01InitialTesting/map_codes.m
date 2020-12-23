clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get a list of aircraft codes present in the IAGOS dataset
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/11/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot flight departure and arrival airports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('airport_codes.mat')

Lats = [];
Lons = [];
for iCoord = 1:1:numel(Coords)
  
  C = strsplit(Coords{iCoord});
  Lons(end+1) = str2num(C{1});
  Lats(end+1) = str2num(C{2});

end
clear iCoord C Coords

clf
m_proj('robinson');

for iA=1:1:numel(Lons)
  
  m_text(Lons(iA),Lats(iA),Codes{iA},'fontsize',6);
  hold on
  
end
m_coast('color','k');
m_grid;

