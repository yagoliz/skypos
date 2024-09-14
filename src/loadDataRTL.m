function [aircraft, sensors] = loadDataRTL(mDir, mCsv, sCsv, f, fraction)
% Load aircraft and sensor data, and compute the TDOA for all the pairs of
% sensors and aircraft

% We will always filter data
if nargin < 4
    f = true;
end

% Start by reading all the dataset
aircraft = readtable(sprintf('%s/%s', mDir, mCsv));
indices = randsample(1:height(aircraft),round(fraction*height(aircraft)));
aircraft = aircraft(indices,:);
sensors = readtable(sprintf('%s/%s', mDir, sCsv));

% This is where filtering starts
if f
    % Filter by sensor type (Radarcape has GPS, thus they are well
    % localized and synchronized
    
    % Filter aircrafts
    % First by NaN
    [row, ~] = find(~isnan([aircraft.latitude, aircraft.longitude, aircraft.geoAltitude]));
    uniquerow = unique(row);
    aircraft = aircraft(uniquerow,:);
    
    % Now remove rows that have not been seen by our sensor set
    serials = sensors.serial;
    keep = zeros(size(aircraft,1),1) == 1;
    measurements_sensors = cell(size(aircraft,1),1);
    measurements_ts = cell(size(aircraft,1),1);
    measurements = aircraft.measurements;
    N = size(aircraft,1);
    for ii = 1:N
        measurement = str2num(measurements{ii});
        measurement = reshape(measurement,3,[])';
        measurement_serials = measurement(:,1);
        
        % If all serials are in the radarcape set, we will keep this
        % measurement
        radarcape_serials = intersect(measurement_serials, serials);
        if numel(radarcape_serials) == numel(measurement_serials)
            keep(ii) = true;
            measurements_sensors{ii} = measurement_serials';
            measurements_ts{ii} = measurement(:,2)';
        end
    end
    
    aircraft.sensors = measurements_sensors;
    aircraft.ts = measurements_ts;
    aircraft.measurements = [];
    aircraft = aircraft(keep,:);
end
