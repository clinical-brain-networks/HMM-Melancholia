
function [mat, T] = prepare_hmm(search_pattern)

% d=dir('../../timeseries/*movie*.txt');
d=dir(search_pattern);


T = [];
mat = [];
% all_K = [12 14 16];

% for i_this_K = 1:numel(all_K)
% this_K = all_K(i_this_K);

for i_d=1:numel(d)
    m = load([d(i_d).folder filesep d(i_d).name]);
    
    % remove first 10 datapoints
    m(1:10,:) = [];
    
    % how many elements?
    T = [T size(m, 1)];
    
    % detrend
    m = detrend(m, 'constant');
    
    % rescape with the SD (each roi individually)
    
    m = m ./ (ones(size(m, 1), 1) * std(m, [], 1));
    
    mat = [mat; m];
    
end






