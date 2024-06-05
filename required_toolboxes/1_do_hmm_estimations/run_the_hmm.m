function run_the_hmm(mat, T, to_write_folder, Ks, reps)

if ~exist(to_write_folder,'dir')
    mkdir(to_write_folder);
end


for K = Ks
    
    
    % K = 10; % no. states
    % reps = 5; % to run it multiple times (saves all the results
    % as a seperate mat file)
    
    TR = 0.81;  % TR of rest data
    
    % options for model estimation; I have set the usual choices here;
    % no need to change anything in geenral
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
    
    % run the HMM multiple times
    
    
    for r = 1:reps
        disp('doing: ');
        disp([K, r]);
        
        
        try
            
            disp(['RUN ' num2str(r)]);
            [hmm, Gamma, ~, vpath, ~, ~, fe] = hmmmar(mat,T,options);
            save([to_write_folder '/HMMrun_K' num2str(K) '_rep_' num2str(r) '.mat'],'Gamma','vpath',...
                'hmm','T','K','fe', 'options');
            
        catch
            disp('failed!')
            disp([K, r]);
            % keyboard;
            disp(lasterr);
            
            
        end
        
    end
    
end



