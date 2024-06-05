
function [lower_lims, upper_lims, this_average, count_matrix_all, observed_count] = make_null_for_vpath_average(m, sort_index, time_lag)





% all_numbers = 1:size(,m 1);
this_sort = sort_index;

nperm = 1000;
fprintf('permuting memberships, making summary vpath null distribution:\n');
X=5;

count_matrix_all = [];
for i_perm = 1:(nperm+1)
    
    % pick random
    perm_this_sort = randperm(size(m, 1), numel(this_sort));

    if i_perm == nperm+1
        perm_this_sort = this_sort;
    end
    

    new_m = m(perm_this_sort, :);
    fprintf('.');
    if i_perm > nperm/100*X
        fprintf('%d', X);
        fprintf('\n');
        X = X + 5;
    end
    
    
    
    outs=[];
    outc=[];
    for im=1:size(new_m,2)
        
        b=im-floor(time_lag/2);
        e=im+floor(time_lag/2);
        
        
        if b < 1 || e > size(new_m,2)
            out(im)=NaN;
            
        else
            
            
            all_vals=new_m(:,b:e);
            
            uv=unique(all_vals);
            ct=[];
            for iuv=1:numel(uv)
                ct(end+1)=sum(sum(all_vals==uv(iuv),2)>0);
            end
            
            um=sortrows([uv ct'],2,'descend');
            
            outs(end+1,:)=um(1:2,1);
            outc(end+1,:)=um(1:2,2);
            
        end
    end
    count_matrix_all(:, end+1) = outc(:, 1);
end

observed_count = count_matrix_all(:, end);
count_matrix_all(:, end) = [];


% keyboard;


% upper_lims = mean(count_matrix_all') + 1.7*std(count_matrix_all');
% lower_lims = mean(count_matrix_all') - 1.7*std(count_matrix_all');
this_average = mean(count_matrix_all');

lower_lims = prctile(count_matrix_all',5);
upper_lims = prctile(count_matrix_all',95);

% 
% figure;
% plot(this_average,'k');
% hold on;
% plot(upper_lims,'r');
% plot(lower_lims,'r');
% plot(observed_count,'k','linewidth',2)
% 
% plot(prctile(count_matrix_all',5),'m');
% plot(prctile(count_matrix_all',95),'m');



% figure;
% axes;
% patch([1:numel(lower_lims) numel(lower_lims):-1:1], [lower_lims, upper_lims(end:-1:1)],'b')
% hold on;
% plot(this_average,'k','linewidth',2)
% plot(observed_count,'r','linewidth',2)
% set(gcf,'position',[39         758        1867          80]);
% title(sprintf('%d',numel(uv)));
% 
% 
% figure;plot(observed_count - this_average');
% hold on;
% plot( 1.7*std(count_matrix_all'));
% plot( -1.7*std(count_matrix_all'));
% set(gcf,'position', [39         618        1867          80]);
% set(gca,'position',[0 0 1 1]);
% title(sprintf('%d',numel(uv)));
% 



    