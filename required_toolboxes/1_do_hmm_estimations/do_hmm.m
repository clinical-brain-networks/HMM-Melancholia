n_sub = numel(subs_to_use); % no. of subjects(used for concatenation)
K = this_NSTATES; % no. states
reps = 5; % to run it multiple times (saves all the results
% as a seperate mat file)

TR = 0.81;  % TR of timeseries
K=10;

% options for model estimation; I have set the usual choices here
% for explanation see the HMM-MAR Wiki page at:
% https://github.com/OHBA-analysis/HMM-MAR/wiki
options = struct();
options.K = K; % number of states
options.order = 0; % no autoregressive components
options.zeromean = 0; % model the mean
options.covtype = 'full'; % full covariance matrix
options.Fs = 1/TR;
options.verbose = 1;
options.standardise = 1;
options.inittype = 'HMM-MAR';
options.cyc = 500;
options.initcyc = 10;
options.initrep = 3;

% run the HMM multiple times - and save the results to .mat
% files



for K=[2:2:20]
    for r = 1:reps
        disp(['RUN ' num2str(r)]);
        [hmm, Gamma, ~, vpath, ~, ~, fe] = hmmmar(dat,T,options);
        save([save_dir '/HMMrun_rep_' num2str(r) '.mat'],'Gamma','vpath',...
            'hmm','T','J','K','n_sub','fe', 'options');

    end
end