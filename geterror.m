unknownllh = ecef2llh(position(1:3)'+center);
realposllh = table2array(sensor(:,{'latitude','longitude','height'}));
error = havdist(unknownllh, realposllh)