function [err, grad] = minimizeroffsetdrift(x, y, z, o, drift, serial, dtime, aircraft, data, ground)

% Prepare data in tables 
sensor = array2table([x, y, z, o, drift, serial], 'VariableNames', {'xc', 'yc', 'zc', 'offsets', 'drift', 'serial'});
comb = [sensor;ground];

if height(data) < 5
%     warning('Not enough data for sensor');
    err = 0; grad = [0;0;0;0;0];
    return
end

[~,acId] = ismember(data.id, aircraft.id, 'rows');
[~,s1Id] = ismember(data.s1, comb.serial, 'rows');
[~,s2Id] = ismember(data.s2, comb.serial, 'rows');

acd = table2array(aircraft(acId,{'xc', 'yc', 'zc'}));
s1d = table2array(comb(s1Id,{'xc', 'yc', 'zc'}));
s2d = table2array(comb(s2Id,{'xc', 'yc', 'zc'}));
s1o = comb.offsets(s1Id);
s2o = comb.offsets(s2Id);
d1o = comb.drift(s1Id);
d2o = comb.drift(s2Id);

% Computation of TDOAs
d1 = ecef_distance(acd, s1d);
d2 = ecef_distance(acd, s2d);

exp = d1 - d2 + (s1o - s2o) + d1o.*dtime - d2o.*dtime;
res = exp - data.ddoam;

data.d1 = d1;
data.d2 = d2;
data.o1 = s1o;
data.o2 = s2o;
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
     
    do = sum(2 .* res(m1)) - sum(2 .* res(m2));
    
    dt = sum(2 .* res(m1) .* dtime(m1)) - sum(2 .* res(m2) .*dtime(m2));
    
    grad = [dx; dy; dz; do; dt]';
end
end
