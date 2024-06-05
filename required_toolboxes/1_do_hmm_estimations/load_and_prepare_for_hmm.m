
function [mat, T] = load_and_prepare_for_hmm(search_pattern)

d=dir(search_pattern);


T = [];
mat = [];

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
    
    % figure;plot(m);
    
end
    
    


