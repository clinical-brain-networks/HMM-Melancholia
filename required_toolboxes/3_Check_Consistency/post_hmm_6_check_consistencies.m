% we need to define the iterations where all 10 states were expressed

% I need the directory of all the outputs, and K, the # of states to
% consider
function post_hmm_6_check_consistencies(in_dir, K)

% in_dir = 'movie';
% K = 10;

d=dir(sprintf('%s/HMMrun_K%d_rep_*.mat', in_dir, K));


% iteration's states.
clear mapping
mapping=struct();
mapping(:)=[];

% extract the numbers of the iterations:
iters = cell2mat(cellfun(@(x) str2num(regexprep(x, '.*_rep_([0-9]+).*','$1')), {d(:).name}', 'uniformoutput', false))';

for i=1:numel(iters)
    
    fprintf('Doing: %s\n', [d(i).folder filesep d(i).name]);
    
    % load the states of i
    hmmvars_i=load([d(i).folder filesep d(i).name]);
    allw_i=[];
    for ii=1:K
        allw_i=[allw_i hmmvars_i.hmm.state(ii).W.Mu_W'];
    end
    
    mapping(i).i = i;
    mapping(i).iter_val = iters(i);
    
    mapping(i).assignmatrix=zeros(K, numel(iters));
    
    mapping(i).corrvals =zeros(K, numel(iters));
    mapping(i).assign_iter_vals = zeros(1, numel(iters));
    
    vpath_i = hmmvars_i.vpath;
    
    for j=setxor(1:numel(iters), i)
        
        % load the states of i
        hmmvars_j=load([d(j).folder filesep d(j).name]);
        
        vpath_j = hmmvars_j.vpath;
        allw_j=[];
        for jj=1:K
            allw_j=[allw_j hmmvars_j.hmm.state(jj).W.Mu_W'];
        end
        
        % bookkeeping
        mapping(i).assign_iter_vals(j) = iters(j);
        
        % do the Hungarian Algorithm to find closest match in the other
        % set; cost = 1-corr.
        
        % calculate the cost(s) associated iwth the 2 vpaths:
        JO=[];
        for iState=1:K
            for jState=1:K
                binned_i = vpath_i==iState;
                binned_j = vpath_j==jState;
                JO(iState, jState) = numel(intersect(find(binned_i),find(binned_j))) / numel(union(find(binned_i),find(binned_j)));
            end
        end
        
        % this is the way to go for RS analyses. OR (!) For analyses if high-order brain regions
        % that do computations that are irreconcilable with outside perturbations, but rather with internal
        % processes such as 'what do I make of this', 'increase vigilance, htis might get interesting',etc...
        % We cannot assume that everyone is in the same state at the same time. Therefore, the
        % Vpath kind of loses its significance.
        
        % [reassignments, costs] = munkres((1-corr(allw_i, allw_j)));
        % [reassignments, costs] = munkres((1-corr(allw_i, allw_j)) + (1-JO));
        [reassignments, costs] = munkres((1-JO));
        
        % just check if I didn't get it wrong...
        % tmp3=sortrows([1:10; reassignments]',2);
        % alternative = tmp3(:,1);
        % reassignments = alternative;
        
        % bookeeping
        mapping(i).assignmatrix(:, j) = reassignments;
        
        % reshuffle the allw_j's
        reshuffled_all_j = allw_j(:, reassignments);
        
        % calculate correlation between allw_i and this allw_j:
        mapping(i).corrvals(:, j) = diag(corr(allw_i, reshuffled_all_j));
        % keyboard;
        
        
        mapping(i).vpaths(:, i) = vpath_i;  % we do not reshuffle the i
        
        % reshuffle the vpath according to reassignments:
        reshuffled_vpath_j = [];
        tmp3=sortrows([1:K; reassignments]',2);
        alternative = tmp3(:,1);
        for ir=1:numel(reassignments)
            reshuffled_vpath_j(vpath_j==ir) = alternative(ir);
        end
        
        mapping(i).vpaths(:, j) = reshuffled_vpath_j;
        mapping(i).unshuffled_vpaths(:, j) = vpath_j;
        
        
    end
    
    % average correlation with matched state in the other hmm inference:
    tmp=mapping(i).corrvals;
    
    
    mapping(i).avg_corr_with_others = mean(tmp);
    tmp2=mean(tmp);
    tmp2=mean(tmp2(tmp2~=0));
    mapping(i).overall_corr = tmp2;
    
end


y_values = [mapping(:).overall_corr];


max_y_value = y_values(y_values == max(y_values));
i_max_y_value = find(y_values == max(y_values));
i_file = mapping(i_max_y_value).iter_val;


fh=figure;lh=plot(y_values);
hold on;
set(lh,'linewidth',2, 'marker','+','color','k');
xlim([0 numel(iters)+1]);
xticks(1:numel(iters));
xticklabels([mapping.iter_val]');
xtickangle(90);
plot(i_max_y_value, max_y_value, 'rO', 'linewidth',2, 'markersize', 20);
lh = line(get(gca,'xlim'), [1 1] * max_y_value);
set(lh,'linestyle','--','color',[0.4 0.4 0.4]);

box off

xl=xlabel(sprintf('iteration file number: X in HMMrun_K%d_rep_X.mat', K));
set(xl,'interpreter','none');

old_pos = get(gcf,'position');
set(gcf, 'position', [old_pos(1:2) 1106         476]);


ylabel('Average correlation with other inferences');
title(sprintf('Consistency of BOLD loadings of HMM inferences across Inferences\nMax corr %2.2f at iteration: %d\nStates are matched between Inferences using Jaccard overlap between State Paths\nCorrelations are calculated using weightings\nNSTATES = %d', max_y_value, i_file, K));
set(gcf,'color','w');
set(gcf,'paperunits','centimeters');
set(gcf,'papersize', [25 20]);
set(gcf,'paperposition',[0 0 25 20]);

saveas(gcf,sprintf('%s/figures/hmm_inferences_consistency_K%d.pdf', in_dir, K)); % sprintf('hmm_inferences_consistency_K%d.pdf', K));


save(sprintf('HungarianMapping_K%d.mat', K),'mapping');
