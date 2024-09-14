function [point_x, point_y, mse_doa ] = htmp(doa_meters, rx, xrange, yrange, combinations, resolution)

%    returns: 
%    heat_long: longitudes of heatmap points
%    heat_lat: latitudes of heatmap points
%    mse_doa: heatmap magnitudes

    if nargin == 7
       resolution = 100; 
    end

    disp('creating heatmap... ');

    num_points = resolution; % points in one dimension (creates squared area)

    % defines the area, where the heatmap is displayed (around the geodetic
    % reference point)
    xmin = xrange(1);
    xmax = xrange(2);
   
    ymin = yrange(1);
    ymax = yrange(2);
    
    % create heatmap
    heat_x  = linspace(xmin,  xmax,  num_points);
    heat_y = linspace(ymin, ymax, num_points);
    point_x = zeros(num_points);
    point_y = zeros(num_points);
    mse_doa = zeros(num_points);
    
    for idx = 1:num_points
        for idy = 1:num_points
            % calculate mean squared error of current point in terms of tdoa

            % distance current point to receivers
            dist_to_rxs = sqrt(sum(([heat_x(idx), heat_y(idy)] - rx).^2, 2));
            
            % current doa in meters
            current_doa = zeros(length(combinations),1);
            for ii = 1:length(combinations)
                s1 = combinations(ii,1);
                s2 = combinations(ii,2);
                current_doa(ii) = dist_to_rxs(s1) - dist_to_rxs(s2);
            end
            
            % error doa
            doa_error = sum((current_doa - doa_meters).^2);
            point_x(idx, idy) = heat_x(idx);
            point_y(idx, idy) = heat_y(idy);
            mse_doa(idx, idy) = doa_error;
        end
    end
    
    mse_doa = 10.*log10(1./mse_doa);
    m = min(min(mse_doa));
    M = max(max(mse_doa));
    mse_doa = (mse_doa -m)./(M-m);
   
    disp('creating heatmap done! ');
    point_x = point_x(:);
    point_y = point_y(:);
    mse_doa = mse_doa(:);
end
