function pval = permtest_aux_NPC(X,D,Nperm,confounds,pairs,index_subjects)
% permutation testing routine, integrating tests through the NPC algorithm
% X: data
% D: design matrix
% Nperm: no. of permutations
% confounds
% grouping: the first grouping(1) rows belong to one group, 
%    the next grouping(2) rows belong to the second group, 
%    and so on and so on - DEPRECATED
%
% See Vidaurre et al. (2018) Human Brain Mapping
%
% Diego Vidaurre, University of Oxford (2018)


[N,p] = size(X); q = size(D,2);  
if (nargin>3) && ~isempty(confounds)
    confounds = confounds - repmat(mean(confounds),N,1);
    X = bsxfun(@minus,X,mean(X));   
    X = X - confounds * pinv(confounds) * X;
    %D = bsxfun(@minus,D,mean(D));   
    %D = D - confounds * pinv(confounds) * D;    
end
paired = (nargin>4) && ~isempty(pairs);
permute_per_subject = (nargin>5) && ~isempty(index_subjects);

%D = bsxfun(@minus,D,mean(D));   
X = zscore(X);
%proj = (D' * D + 0.001 * eye(size(D,2))) \ D';  
%rss0 = sum(X.^2);

T = zeros(Nperm,1);
for perm=1:Nperm
    if perm==1
        Xin = X;
    elseif permute_per_subject
        Xin = zeros(size(X));
        Nsubj = max(index_subjects);
        for j = 1:Nsubj
           jj = find(index_subjects == j); 
           Ntrials = length(jj);
           r = randperm(Ntrials);
           Xin(jj,:) = X(jj(r),:);
        end
    elseif paired
        r = 1:N;
        for j = 1:max(pairs)
            jj = find(pairs==j);
            if length(jj) < 2, continue; end
            if rand<0.5
                r(jj) = r(jj([2 1]));
            end
        end
        Xin = X(r,:);
    else
        Xin = X(randperm(N),:);
    end
    pv = zeros(1,p);
    %beta = proj * Xin;
    for i = 1:p
        [~,~,~,~,pv_i] = regress(Xin(:,i),D);
        pv(i) = pv_i(3);
        %rss = sum((D * beta(:,i) - Xin(:,i)).^2);
        %F = (rss0(i)-rss)/rss * ((N-q)/q) ; 
        %pv(i) = 1 - fcdf(F,q,N-q);
    end
    T(perm) = -2 * sum(log(pv));
end

pval = sum( T(1) <= T ) ./ (1 + Nperm) ;

end


