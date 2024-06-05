
function [mat, T] = load_and_prepare_for_hmm(search_pattern, varargin)


max_size = Inf;
cut_first = 10;

if numel(varargin) > 0
    cut_first = varargin{1};
    max_size = varargin{2};
end

% d=dir('../../timeseries/*movie*.txt');
d=dir(search_pattern);


T = [];
mat = [];
% all_K = [12 14 16];

% for i_this_K = 1:numel(all_K)
% this_K = all_K(i_this_K);

for i_d=1:numel(d)
    disp(d(i_d).name);
    m = load([d(i_d).folder filesep d(i_d).name]);
    
    
    % to take into account our lower-than-expected-person
    if size(m,1) > max_size
        m = m(1:max_size,:);
    end
    
    
    % remove first 10 datapoints
    m(1:cut_first,:) = [];
    
    
    
    
    % how many elements?
    T = [T size(m, 1)];
    
    % detrend
    m = detrend(m, 'constant');
    
    % rescape with the SD (each roi individually)
    
    m = m ./ (ones(size(m, 1), 1) * std(m, [], 1));
    
    mat = [mat; m];
    
end






