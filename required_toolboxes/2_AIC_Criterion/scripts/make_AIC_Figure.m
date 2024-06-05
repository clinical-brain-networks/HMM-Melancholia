% Load in for each inference all the HMM's, so we can work with those


basedir = '/mnt/data/johanv/Phil/analysis/2_AIC_Criterion/';
ds=struct();
% hmmdata=load(sprintf('%s/results_%d/aroma/all/HMMrun_rep_%d_data.mat',i_inference,i_run));
whichrun = 'AIChmms_for_movie';

[mat, T] = prepare_hmm('../../timeseries/*movie*.txt');

for i_inference=[2:2:38]

    all_existing_files = dir(sprintf('%s/%s/HMMrun_K%d_rep_*.mat',basedir, whichrun, i_inference));
    all_existing_runs = cell2mat(cellfun(@(x) str2num(regexprep(x, '.*_([0-9]+).mat$','$1')), {all_existing_files.name}, 'uniformoutput', false));
    
    for i_run=all_existing_runs
        
        hmmvars=load(sprintf('%s/%s/HMMrun_K%d_rep_%d.mat',basedir, whichrun, i_inference,i_run));
        % summeasures=load(sprintf('%s/results_%d/%s/Summary_measures_rep_%d.mat',basedir, i_inference,i_run));
        
        ds(i_inference,i_run).hmmvars=hmmvars;
        % ds(i_inference,i_run).summeasures=summeasures;
        
    end
end



% addpath(genpath('/home/johan/mnt/hpcworking/projects/hmm-movie/hmm-mar-matlab/revisions/R12/HMM-MAR'));


output=struct();

for i_inference=[2:2:38]
    
    % take the first inference of the HMM where the maximal number of
    % states are occupied (which are the NSTATES.. usually) -- note down
    % the number of states too.
    all_existing_files = dir(sprintf('%s/%s/HMMrun_K%d_rep_*.mat',basedir, whichrun, i_inference));
    all_existing_runs = cell2mat(cellfun(@(x) str2num(regexprep(x, '.*_([0-9]+).mat$','$1')), {all_existing_files.name}, 'uniformoutput', false));
   
    
    uvals=[];
    for i_run=all_existing_runs
        vpath=ds(i_inference,i_run).hmmvars.vpath;
        
        unique_vals = numel(unique(vpath));
        uvals(end+1) = unique_vals;
    end
    
    
    % another fix:
    which_inference_to_choose = all_existing_runs(find(uvals==max(uvals),1,'first'));
    NSTATES = uvals(find(uvals==max(uvals),1,'first'));
    

    % now we chose (out of 15) the inference - calculate its fe:
    T = ds(i_inference,which_inference_to_choose).hmmvars.T;

    data = mat;
    hmm=ds(i_inference,which_inference_to_choose).hmmvars.hmm;
    
    
    
    [fe ll]=hmmfe(data,T,hmm);
    
    
    output(i_inference).which_inference_to_choose = which_inference_to_choose;
    output(i_inference).NSTATES = NSTATES;
    output(i_inference).fe = fe; % free energy;
    output(i_inference).AIC_option_FreeEnergy = 2*i_inference - 2*log(fe);
    output(i_inference).AIC_option_MeanLLPerTimePoint = 2*i_inference - 2*mean(ll);
    output(i_inference).meanll = mean(ll);
    
    nstates_option_2 = NSTATES*(NSTATES-1) + (NSTATES-1) + (NSTATES)*13;
    nstates_option_3 = 2*NSTATES*14;
    
    
    output(i_inference).AIC_option_2_MeanLLPerTimePoint = nstates_option_2 - 2*mean(ll);
    output(i_inference).AIC_option_3_MeanLLPerTimePoint = nstates_option_3 - 2*mean(ll);
    
    output(i_inference).AIC_option_2_SumLLPerTimePoint = nstates_option_2 - 2*sum(ll);
    output(i_inference).AIC_option_3_SumLLPerTimePoint = nstates_option_3 - 2*sum(ll);
    
    output(i_inference).AIC_option_NSTATEDIV14 = 2*NSTATES/14 - 2*mean(ll);
    
    % thematrix(end+1,:) = [i_inference, which_inference_to_choose, NSTATES, fe, 2*i_inference - 2*log(fe), 2*i_inference - 2*mean(ll), mean(ll)];
    
    
end

figure;
subplot(2,2,1);
plot([2:2:38], [output(:).AIC_option_2_SumLLPerTimePoint])
title('sum of the LL per time point');
subplot(2,2,3);
plot([2:2:38], -1*100*[NaN diff([output(:).AIC_option_2_SumLLPerTimePoint])]./[output(:).AIC_option_2_SumLLPerTimePoint])
title('% change of the sum of the LL per increase step in K');
subplot(2,2,2);
plot([2:2:38], [output(:).fe])
title('fe');
subplot(2,2,4);
plot([2:2:38], -1*100*[NaN diff([output(:).fe])] ./ [output(:).fe]);
title('% change fe per increase step in K');



    
    
    
    