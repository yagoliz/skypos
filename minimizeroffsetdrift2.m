function [err, grad] = minimizeroffsetdrift2(x, y, z, o, drift, serial, dtime, aircraft, data, ground)

% Prepare data in tables 
sensor = array2table([x, y, z, o, drift, serial], 'VariableNames', {'x', 'y', 'z', 'offsets', 'drift', 'serial'});
comb = [sensor;ground];

if height(data) < 4
%     warning('Not enough data for sensor');
    err = 0; grad = [0;0;0;0];
    return
end

[~,acId] = ismember(data.id, aircraft.id, 'rows');
[~,s1Id] = ismember(data.s1, comb.serial, 'rows');
[~,s2Id] = ismember(data.s2, comb.serial, 'rows');

acd = table2array(aircraft(acId,{'x', 'y', 'z'}));
s1d = table2array(comb(s1Id,{'x', 'y', 'z'}));
s2d = table2array(comb(s2Id,{'x', 'y', 'z'}));
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
    epsilon = 1e-10;
    
    xp = x*(1+epsilon); xm = x*(1-epsilon);
    dx = (minimizeroffsetdrift2(xp,y,z,o,drift,serial,dtime,aircraft,data,ground) ...
            - minimizeroffsetdrift2(xm,y,z,o,drift,serial,dtime,aircraft,data,ground)) / (2 * x * epsilon);
        
    yp = y*(1+epsilon); ym = y*(1-epsilon);
    dy = (minimizeroffsetdrift2(x,yp,z,o,drift,serial,dtime,aircraft,data,ground) ...
            - minimizeroffsetdrift2(x,ym,z,o,drift,serial,dtime,aircraft,data,ground)) / (2 * y * epsilon);
        
    zp = z*(1+epsilon); zm = z*(1-epsilon);
    dz = (minimizeroffsetdrift2(x,y,zp,o,drift,serial,dtime,aircraft,data,ground) ...
            - minimizeroffsetdrift2(x,y,zm,o,drift,serial,dtime,aircraft,data,ground)) / (2 * z * epsilon);
        
    op = o*(1+epsilon); om = o*(1-epsilon);
    do = (minimizeroffsetdrift2(x,y,z,op,drift,serial,dtime,aircraft,data,ground) ...
            - minimizeroffsetdrift2(x,y,z,om,drift,serial,dtime,aircraft,data,ground)) / (2 * o * epsilon);
        
    dtp = drift*(1+epsilon); dtm = drift*(1-epsilon);
    dt = (minimizeroffsetdrift2(x,y,z,o,dtp,serial,dtime,aircraft,data,ground) ...
            - minimizeroffsetdrift2(x,y,z,o,dtm,serial,dtime,aircraft,data,ground)) / (2 * drift * epsilon);
    
    
    grad = [dx; dy; dz; do; dt]';
end
end