% clear all;close all;clc;
addpath(genpath('../../scripts/NBSdirected'));


file_in = 'movie/HMMrun_K10_rep_50.mat';

ismel = find(v==1);
isnonmel = find(v==0);

% COMP_A = 'mov1a';
% COMP_B = 'resta';

% load valid_inferences_all.mat
% RUN=valid_inferences(1); % you might change this depending on which HMM inference you wish to check out.


% ANALYSIS='all';

% run=sprintf('run%d',RUN);
% subs_to_use_real = [2     3     7     8     9    10    11    12    14    16    17    18    19    20];


out_mel = post_hmm_4_investigate_state_stats(file_in, ismel);
out_nonmel = post_hmm_4_investigate_state_stats(file_in, isnonmel);



out_mat_a = [];
for i=1:numel(ismel)
    out_mat_a = cat(3, out_mat_a, out_mel.trans_matrices{i});
    diag_elelents = out_mel.trans_diags{i};
    for j=1:size(out_mat_a, 1)
        out_mat_a(j,j,end) = 0; %diag_elelents(j);
    end
end


out_mat_b = [];
for i=1:numel(isnonmel)
    out_mat_b = cat(3, out_mat_b, out_nonmel.trans_matrices{i});
    diag_elelents = out_nonmel.trans_diags{i};
    for j=1:size(out_mat_b, 1)
        out_mat_b(j,j,end) = 0; %diag_elelents(j);
    end
end






C=cat(3,out_mat_a,out_mat_b);

nsubs=size(C,3);

% design matrix:
X=[ones(numel(ismel),1) zeros(numel(ismel),1); zeros(numel(isnonmel),1) ones(numel(isnonmel),1)];  % NOTE: Modeling the subject variance yields similar results; see also Exchange Blocks (below)




%Any significant results are stored in variable called con_mat

%Connectivity matrices (regions x regions x subjects)
% C=randn(10,10,6);
% C(1,2,1:3)=C(1,2,1:3)+10;
% C(1,3,1:3)=C(1,3,1:3)+10;
% C(1,4,1:3)=C(1,4,1:3)+10;
% C(4,3,1:3)=C(4,3,1:3)+10;
% C(4,1,1:3)=C(4,3,1:3)+10;

%Total number of permutations to generate
GLM.perms=5000;
%Design matrix
GLM.X=X;
 %Contrast
 GLM.contrast=[1 -1 ]; 
 %Type of test
 GLM.test='ttest'; % 'ttest' or 'ftest'
 %Exchange block for repeated measures. NOTE: Constraining the variance to within-subject effects produces highly similar results.
 %GLM.exchange=[]; STAT
 STATS.size='Extent'; %'Intensity' or 'Extent'
 %Threshold
 STATS.thresh=1.0; 
 %Significance (usually 0.05)
 STATS.alpha=0.8; 

 %NO NEED TO CHANGE ANYTHING BEYOND THIS POINT
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 bgl=0; 
 nDisp=100; 
 
 %Number of nodes
 N=size(C,1);
 ind_uplo=union(find(triu(ones(N,N),1)),find(tril(ones(N,N),-1)));
 
 GLM.y=zeros(size(GLM.X,1),length(ind_uplo));
 for i=1:size(C,3)
     tmp=C(:,:,i);
     GLM.y(i,:)=tmp(ind_uplo);
 end

%Precompute test stat
STATS.test_stat=NBSglm(GLM); 

%Number of edges
J=length(ind_uplo); 

%Determine whether test statistics have been precomputed and determine
%index of edges exceeding the primary threshold
if ~isempty(STATS.test_stat)
    %Precomputed test statistics
    ind=ind_uplo(STATS.test_stat(1,:)>STATS.thresh); 
    %Number of permutations
    K=size(STATS.test_stat,1)-1; 
else
    %Never get to this case
end

%Size of a component measured using extent or intensity? 
Intensity=0;
if strcmp(STATS.size,'Intensity')    
    %If size measure using intensity, create an N x N matrix cotaining the 
    %test statistic for each edge minus the test statistic threshold
    %(primary threshold)
    Intensity=1; 
    %Compute a test statistic matrix
    test_stat_mat=zeros(N,N); 
    if ~isempty(STATS.test_stat)
        %Precomputed
        test_stat_mat(ind_uplo)=STATS.test_stat(1,:)-STATS.thresh;
%       test_stat_mat=(test_stat_mat+test_stat_mat');
    else
        %Never reach this case. 
    end
end
    
adj=spalloc(N,N,length(ind));
adj(ind)=1; 
%Only consider components comprising more than one node, equivalent to at
%least one edge
if bgl==1
    [a,sz]=components((adj+adj')/2); 
else
    [a,sz]=get_components((adj+adj')/2); 
end
ind_sz=find(sz>1);
sz_links=zeros(1,length(ind_sz));
max_sz=0; 
for i=1:length(ind_sz)
    nodes=find(ind_sz(i)==a);
    if Intensity
        %Measure size as intensity
        sz_links(i)=sum(sum(adj(nodes,nodes).*test_stat_mat(nodes,nodes))); %/2;
    else
        %Measure size as extent
        sz_links(i)=sum(sum(adj(nodes,nodes))); %/2;
    end
    adj(nodes,nodes)=adj(nodes,nodes)*(i+1);
    if max_sz<sz_links(i)
        max_sz=sz_links(i);
    end
end

%Subtract one to remove edges not part of a component
%Although one is also subtracted from edges comprising a component, this is 
%compensated by the (i+1) above
adj(~~adj)=adj(~~adj)-1;

%Repeat above for each permutation
%Empirical null distribution of maximum component size
null_dist=zeros(K,1); 
str1='| Permutation | Max Size | Max Size | Lowest  |';
str2='|             |  Random  |  Actual  | p-value |';
try tmp=get(H,'string'); set(H,'string',[{str1};{str2};tmp]); drawnow;
catch;  fprintf([str1,'\n',str2,'\n']); end 
p_approx=0;
%Store what is already displayed in the listbox
try pre_str=get(H,'string'); catch; end
new_str={};
%First row of test_stat is the observed test statistics, so start at the
%second row
for i=2:K+1
    if ~isempty(STATS.test_stat)
        %Precomputed test statistics 
        ind=ind_uplo(STATS.test_stat(i,:)>STATS.thresh); 
    else

    end
    if Intensity 
        %Compute a test statistic matrix
        test_stat_mat=zeros(N,N); 
        if ~isempty(STATS.test_stat)
            test_stat_mat(ind_uplo)=STATS.test_stat(i,:)-STATS.thresh;
        else

        end    
    end
    adj_perm=spalloc(N,N,length(ind));
    adj_perm(ind)=1;
    if bgl==1
        [a,sz]=components((adj_perm+adj_perm')/2); 
    else
        [a,sz]=get_components((adj_perm+adj_perm')/2); 
    end
    ind_sz=find(sz>1);
    max_sz_perm=0; 
    for j=1:length(ind_sz)
        nodes=find(ind_sz(j)==a);
        if Intensity
            tmp=sum(sum(adj_perm(nodes,nodes).*test_stat_mat(nodes,nodes))); %/2;
        else
            tmp=sum(sum(adj_perm(nodes,nodes))); %/2;
        end
        if tmp>max_sz_perm
            max_sz_perm=full(tmp);
        end   
    end
    null_dist(i-1)=max_sz_perm; 
    if max_sz_perm>=max_sz
        p_approx=p_approx+1;
    end
   % str=sprintf('|   %5d/%5d |     %4d |     %4d |   %0.3f |',...
   %v1.1.2 Changed to %6.0f to %6.1f to allow fractional component sizes
   %that arise when component size is measured with intensity. 
   str=sprintf('| %5d/%5d |   %6.1f |   %6.1f |  %0.4f |',...
            i-1,K,max_sz_perm,max_sz,p_approx/(i-1));
        %Display no mare than nDisp most recent permutations
        new_str=[str,{new_str{1:min(nDisp,length(new_str))}}]';
        try set(H,'string',[new_str;pre_str]); drawnow; 
        catch;  fprintf([str,'\n']); end 
end
str1='| Permutation | Max Size | Max Size | Lowest  |';
str2='|             |  Random  |  Actual  | p-value |';
try tmp=get(H,'string'); set(H,'string',[{str1};{str2};tmp]); drawnow;
catch;  fprintf([str1,'\n',str2,'\n']); end 


%Determine components satisfying alpha significance threshold
n_cnt=0; 
for i=1:length(sz_links)
    tmp=sum(null_dist>=sz_links(i))/K;
    if tmp<=STATS.alpha
        n_cnt=n_cnt+1;
        ind=find(adj==i);
        con_mat{n_cnt}=spalloc(N,N,length(ind)*2);
        con_mat{n_cnt}(ind)=1; 
        con_mat{n_cnt}=con_mat{n_cnt};
        pval(n_cnt)=tmp;
    end
end
if n_cnt==0
    pval=[]; con_mat=[]; 
end

m=full(con_mat{1});
if numel(con_mat) > 1
    for i=2:numel(con_mat)
        m = m + con_mat{i};
    end
end

save Step6a_calculate_NBS_mov_rest_matrix_sesa.mat m
