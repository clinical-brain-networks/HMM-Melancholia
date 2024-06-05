function Xi = approximateXi(Gamma,T,par)
if nargin==3
    if ~isstruct(par)
        order = par;
    else
        [~,order] = formorders(par.train.order,par.train.orderoffset,...
            par.train.timelag,par.train.exptimelag);
    end
else
    order = 0;
end
K = size(Gamma,2);
Xi = zeros(sum(T-1-order),K,K);
for j = 1:length(T)
    indG = (1:(T(j)-order)) + sum(T(1:j-1)) - (j-1)*order;
    indXi =  (1:(T(j)-order-1)) + sum(T(1:j-1)) - (j-1)*(order+1);
    for t = 1:length(indXi)
        xi = Gamma(indG(t),:)' * Gamma(indG(t+1),:);
        xi = xi / sum(xi(:));
        Xi(indXi(t),:,:) = xi;
    end
end

end