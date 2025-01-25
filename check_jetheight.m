clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute mean jet centre height and latitude
%
%Corwin Wright, c.wright@bath.ac.uk
%2024/11/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/joc.8095

Settings.LatRange  = [30,70];
Settings.LonRange  = [-20,0];
Settings.PrsRange  = [150,400];
Settings.TimeScale = datenum(2012,1,273:302);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create storage arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%at heights below 300 hPa, the terrain following b-coefficient of ERA5 data
%is less than 0.1 out of 1.0. Thus, we will neglect the lnsp part of the 
%pressure calculation

p = ecmwf_prs_v3(137);
pidx = inrange(p,Settings.PrsRange);

Store = NaN(numel(Settings.TimeScale),numel(Settings.LatRange),numel(pidx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.TimeScale)

  %get data for this day
  Data = rCDF(era5_path(Settings.TimeScale(iDay)));
  
  %daily average and select lat/lon range
  latidx = inrange(Data.latitude, minmax(Settings.LatRange));
  lonidx = inrange(Data.longitude,minmax(Settings.LonRange));
  u = nanmean(Data.u,1);
  u = u(:,pidx,:,:);
  u = u(:,:,latidx,:);
  u = squeeze(u(:,:,:,lonidx));
  
  %take a zonal mean
  u = squeeze(nanmean(u,3));

  %plot as test
  pc(Data.latitude(latidx),p(pidx),u); set(gca,'yscale','log','ydir','reverse')
  drawnow
end




