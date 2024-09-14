ms1 = ismember(data.s1,ground.serial(4)) | ismember(data.s2,ground.serial(4));
m1 = data(ms1,:);
m1 = sortrows(m1,'id');
ac1 = sortrows(aircraft,'id');
msa1 = ismember(ac1.id,m1.id);
t1 = aircraft.timeAtServer(msa1);
hold on; grid on; scatter(t1,m1.err/3e8 * 1e9);

% ms2 = ismember(data.s1,ground.serial(2)) | ismember(data.s2,ground.serial(2));
% m2 = data(ms2,:);
% m2 = sortrows(m2,'id');
% ac2 = sortrows(aircraft,'id');
% msa2 = ismember(ac2.id,m2.id);
% t2 = aircraft.timeAtServer(msa2);
% scatter(t2,m2.err/3e8 * 1e9);

% axis square; legend('Reference sensor 1', 'Reference sensor 2');
% xlabel('Time at server (s)');xlim([min([t1;t2]),max([t1;t2])]);
% axis square;
xlabel('Time at server (s)');xlim([min([t1]),max([t1])]);
% ylim([0,10e-9])
ylabel('Error in expected TDOA and measured TDOA (ns)');


%% When sorted
ms1 = ismember(dt.s1,ground.serial(1)) | ismember(dt.s2,ground.serial(1));
m1 = dt(ms1,:);
m1 = sortrows(m1,'id');
ac1 = sortrows(ac,'id');
msa1 = ismember(ac1.id,m1.id);
t1 = ac1.timeAtServer(msa1);
hold on; grid on; scatter(t1,m1.err);

ms2 = ismember(dt.s1,ground.serial(2)) | ismember(dt.s2,ground.serial(2));
m2 = dt(ms2,:);
m2 = sortrows(m2,'id');
ac2 = sortrows(ac,'id');
msa2 = ismember(ac2.id,m2.id);
t2 = ac2.timeAtServer(msa2);
scatter(t2,m2.err);

axis square; legend('Reference sensor 1', 'Reference sensor 2');
xlabel('Time at server (s)');xlim([min([t1;t2]),max([t1;t2])]);
ylabel('Error in expected TDOA and measured TDOA (m)');

%% Filtered
ms1 = ismember(dtm.s1,ground.serial(1)) | ismember(dtm.s2,ground.serial(1));
m1 = dtm(ms1,:);
m1 = sortrows(m1,'id');
ac1 = sortrows(acm,'id');
msa1 = ismember(ac1.id,m1.id);
t1 = ac1.timeAtServer(msa1);
hold on; grid on; scatter(t1,m1.err/3e8-m1.err(1)/3e8);

% ms2 = ismember(dtm.s1,ground.serial(2)) | ismember(dtm.s2,ground.serial(2));
% m2 = dtm(ms2,:);
% m2 = sortrows(m2,'id');
% ac2 = sortrows(acm,'id');
% msa2 = ismember(ac2.id,m2.id);
% t2 = ac2.timeAtServer(msa2);
% scatter(t2,m2.err/3e8);

% axis square; 
% legend('Reference sensor 1', 'Reference sensor 2');
xlabel('Time at server (s)');
% xlim([min([t1;t2]),max([t1;t2])]); ylim([quantile(m2.err,0.05),quantile(m2.err,0.95)]);
xlim([min([t1]),max([t1])]); 
% ylim([quantile(m1.err/3e8,0.05),quantile(m1.err/3e8,0.95)]);
ylabel('Error in expected TDOA and measured TDOA (s)');

%% Histogram
histogram(sensors.offsets/3e8 * 1e9,'Normalization','probability'); grid on; xlabel('Computed offset (ns)'); ylabel('Probability');ylim([0.,1.]);
