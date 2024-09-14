function ddoa = getddoaRTL(aircraft)
% From this function we can compute the table that contains:
% m | s1 | s2 | tdoa

N = size(aircraft,1);
ddoa_tmp = cell(N,1);

for ii = 1:N
    aircraft_row = aircraft(ii,:);
    combinations = nchoosek(1:aircraft_row.numMeasurements, 2);
    
    tmp_mat = zeros(size(combinations,1),4); % 4 columns of the table
    for jj = 1:size(combinations,1)
        s = aircraft_row.sensors{1};
        t = aircraft_row.ts{1};
        
        s1 = s(combinations(jj,1));
        s2 = s(combinations(jj,2));
        
        if s1 < s2
            ddoa = t(combinations(jj,1)) - t(combinations(jj,2));
        else
            ddoa = t(combinations(jj,2)) - t(combinations(jj,1));
            [s1, s2] = swap(s1,s2);
        end
        tmp_mat(jj,:) = [aircraft_row.id, s1, s2, ddoa];
    end
    ddoa_tmp{ii} = tmp_mat;
end

ddoa = cell2mat(ddoa_tmp);
ddoa = array2table(ddoa,'VariableNames', {'id','s1','s2', 'ddoa'});

end

% Small trick
function [b,a] = swap(a,b)
end