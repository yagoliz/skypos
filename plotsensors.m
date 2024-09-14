% prepare the map
figure
ax = worldmap([20 70],[-30 40]);
hold on
setm(ax, 'MapProjection',"mercator", 'Frame',"off", 'Grid',"off")
%% zoom in
ax.XLim = [-1.7846e+06 1.5681e+06];
ax.YLim = [4.3808e+06 7.1748e+06];
% plot borders / lands
geoshow(ax, 'shapefiles/CNTR_RG_10M_2016_4326.shp', 'FaceColor',[0.85 0.85 0.85]);
% plot good sensor locations (green dots)
scatterm(sensors.latitude, sensors.longitude, [], 'green', 'filled','MarkerEdgeColor',[0 0 0]);
% plot bad sensor locations (yellow dots)
scatterm(ground.latitude, ground.longitude, [], 'yellow', 'filled','MarkerEdgeColor',[0 0 0]);

scatterm(sensor.latitude, sensor.longitude, [], 'red', 'filled','MarkerEdgeColor',[0 0 0]);
