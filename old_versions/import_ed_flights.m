function FlightData = import_ed_flights(filename, dataLines)
%IMPORTFILE Import data from a text file
%  TFMSOCEANIC20170601 = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  TFMSOCEANIC20170601 = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  TFMSoceanic20170601 = importfile("C:\Data\corwin\ed_flights\TFMS_oceanic_20170601.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 27-Jul-2022 13:06:37

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
  dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 20);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["AIRCRAFT_ID", "TFMS_FLIGHT_INDEX", "UNIQUE_FLIGHT_ID", "DATE_UTC", "DEPT_APRT", "DEPT_DATE_UTC", "DEPT_TIME_UTC", "ARR_APRT", "ARR_DATE_UTC", "ARR_TIME_UTC", "USER_CLASS", "PHYSICAL_CLASS", "ACFT_TYPE", "TRACK_POINT_DATE_UTC", "TRACK_POINT_TIME_UTC", "LATITUDE", "LONGITUDE", "ALTITUDE_x100_FT", "GROUNDSPEED_kts", "CENTER"];
opts.VariableTypes = ["double", "double", "double", "datetime", "categorical", "double", "datetime", "categorical", "double", "datetime", "categorical", "categorical", "double", "double", "datetime", "double", "double", "double", "double", "categorical"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["DEPT_APRT", "ARR_APRT", "USER_CLASS", "PHYSICAL_CLASS", "CENTER"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "DATE_UTC", "InputFormat", "dd-MMM-yy");
opts = setvaropts(opts, "DEPT_TIME_UTC", "InputFormat", "HH:mm:ss");
opts = setvaropts(opts, "ARR_TIME_UTC", "InputFormat", "HH:mm:ss");
opts = setvaropts(opts, "TRACK_POINT_TIME_UTC", "InputFormat", "HH:mm:ss");
opts = setvaropts(opts, ["AIRCRAFT_ID", "ACFT_TYPE"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["AIRCRAFT_ID", "ACFT_TYPE"], "ThousandsSeparator", ",");

% Import the data
TFMSoceanic20170601 = readtable(filename, opts);

