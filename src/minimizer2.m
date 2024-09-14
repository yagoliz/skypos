function [err, grad] = minimizer2(x, y, z, serial, aircraft, ddoa, ground)

% Prepare data in tables 
sensor = array2table([x, y, z, serial], 'VariableNames', {'x', 'y', 'z', 'serial'});

mask_s = ismember(ddoa.s1, sensor.serial, 'rows') | ismember(ddoa.s2, sensor.serial, 'rows');
mask_g = ismember(ddoa.s1, ground.serial, 'rows') | ismember(ddoa.s2, ground.serial, 'rows');
mask = mask_s & mask_g;

comb = [sensor;ground];
data = ddoa(mask,:);

if height(data) < 4
    warning('Not enought data for sensor');
    err = 0; grad = [0;0;0];
    return
end

[~,acId] = ismember(data.id, aircraft.id, 'rows');
[~,s1Id] = ismember(data.s1, comb.serial, 'rows');
[~,s2Id] = ismember(data.s2, comb.serial, 'rows');

acd = table2array(aircraft(acId,{'x', 'y', 'z'}));
s1d = table2array(comb(s1Id,{'x', 'y', 'z'}));
s2d = table2array(comb(s2Id,{'x', 'y', 'z'}));

% Computation of TDOAs
d1 = ecef_distance(acd, s1d);
d2 = ecef_distance(acd, s2d);

exp = d1 - d2;
res = exp - data.ddoam;

data.d1 = d1;
data.d2 = d2;
data.res = res;

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
    
    
    s = sensor.serial(1);
    m1 = data.s1 == s;
    m2 = data.s2 == s;

    dx = sum(2 .* res(m1) .* dx1(m1)) + ...
         sum(2 .* res(m2) .* dx2(m2));

    dy = sum(2 .* res(m1) .* dy1(m1)) + ...
         sum(2 .* res(m2) .* dy2(m2));

    dz = sum(2 .* res(m1) .* dz1(m1)) + ...
         sum(2 .* res(m2) .* dz2(m2));
    
    grad = [dx; dy; dz];
end
end