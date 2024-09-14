function [err, grad] = minimizer(x, y, z, serial, aircraft, ddoa)

% Prepare data in tables 
sensors = array2table([x, y, z, serial], 'VariableNames', {'x', 'y', 'z', 'serial'});

[~,acId] = ismember(ddoa.id, aircraft.id, 'rows');
[~,s1Id] = ismember(ddoa.s1, sensors.serial, 'rows');
[~,s2Id] = ismember(ddoa.s2, sensors.serial, 'rows');

acd = table2array(aircraft(acId,{'x', 'y', 'z'}));
s1d = table2array(sensors(s1Id,{'x', 'y', 'z'}));
s2d = table2array(sensors(s2Id,{'x', 'y', 'z'}));

% Computation of TDOAs
d1 = ecef_distance(acd, s1d);
d2 = ecef_distance(acd, s2d);

exp = d1 - d2;
res = exp - ddoa.ddoam;

ddoa.d1 = d1;
ddoa.d2 = d2;
ddoa.res = res;

% Removing outliers (not sure how correct this is)
% lowerlim = quantile(ddoa.res,0.01);
% upperlim = quantile(ddoa.res,0.99);
% mask = ddoa.res > lowerlim & ddoa.res < upperlim;
% 
% data = ddoa(mask,:);
% acd = acd(mask,:);
% s1d = s1d(mask,:);
% s2d = s2d(mask,:);
% res = res(mask);
% d1 = d1(mask);
% d2 = d2(mask);
data = ddoa;

% Computing the actual error
err = sum(res.^2);

% Analytical computation of the gradient (if requested)
if nargout > 1
    % The derivatives are (x_i - x_p)/d1 (for first term)
    % And (x_p - y_i)/d2 (for second)
    dx1 = (s1d(:,1) - acd(:,1)) ./ d1;
    dy1 = (s1d(:,2) - acd(:,2)) ./ d1;
    dz1 = (s1d(:,3) - acd(:,3)) ./ d1;
    
    dx2 = -(s2d(:,1) - acd(:,1)) ./ d2;
    dy2 = -(s2d(:,2) - acd(:,2)) ./ d2;
    dz2 = -(s2d(:,3) - acd(:,3)) ./ d2;
    
    dx = zeros(height(sensors),1);
    dy = zeros(height(sensors),1);
    dz = zeros(height(sensors),1);
    
    for ii = 1:height(sensors)
        s = sensors.serial(ii);
        m1 = data.s1 == s;
        m2 = data.s2 == s;
        
        dx(ii) = sum(2 .* res(m1) .* dx1(m1)) + ...
                 sum(2 .* res(m2) .* dx2(m2));
             
        dy(ii) = sum(2 .* res(m1) .* dy1(m1)) + ...
                 sum(2 .* res(m2) .* dy2(m2));
             
        dz(ii) = sum(2 .* res(m1) .* dz1(m1)) + ...
                 sum(2 .* res(m2) .* dz2(m2));
    end
    
    grad = [dx; dy; dz];
end
end