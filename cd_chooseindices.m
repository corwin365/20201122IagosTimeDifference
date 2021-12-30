function [] = cd_chooseindices(Paths,Settings)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose which indices we want to keep raw and which deseasonalised, and 
%write to a new file to make them easier to call in this form
%Corwin Wright, c.wright@bath.ac.uk, 2021/12/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load all indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([Paths.StoreDir,'/indices_',Paths.SourceIdentifier,'.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% copy over the chosen fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reate structs to keep chosen indices in
Flight = struct();
Daily  = struct();
WhichIsWhich = struct();
Ranges = struct();


%first, copy all indices we intend to use over as raw
%we'll then replace some with deseasonalised values
for iIndex=1:1:numel(Settings.Indices)
  Flight.(Settings.Indices{iIndex}) = FlightIndices.Raw.(Settings.Indices{iIndex});
  Daily.( Settings.Indices{iIndex}) =   DateIndices.Raw.(Settings.Indices{iIndex});
  WhichIsWhich.(Settings.Indices{iIndex}) = 'Raw';
  Ranges.(Settings.Indices{iIndex}) = DateIndices.Raw.Ranges.(Settings.Indices{iIndex});
end

%now, replace with deasonalised if requested
for iIndex=1:1:numel(Settings.ChosenDS)
  Flight.(Settings.ChosenDS{iIndex}) = FlightIndices.DS.(Settings.ChosenDS{iIndex});
  Daily.( Settings.ChosenDS{iIndex}) =   DateIndices.DS.(Settings.ChosenDS{iIndex});
  WhichIsWhich.(Settings.ChosenDS{iIndex}) = 'DS';
  Ranges.(Settings.ChosenDS{iIndex}) = DateIndices.DS.Ranges.(Settings.ChosenDS{iIndex});  
end



%finally, copy over the times
Flight.Time = FlightIndices.Time;
Daily.Time   = DateIndices.Time;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tidy up, rename, and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([Paths.StoreDir,'/indices_',Paths.SourceIdentifier,'_',Paths.PPIdentifier,'.mat'],'Flight','Daily','Ranges','WhichIsWhich')
