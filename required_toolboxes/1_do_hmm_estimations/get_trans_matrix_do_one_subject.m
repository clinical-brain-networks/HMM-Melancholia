function [matrices, diags] =  get_trans_matrix_do_one_subject(vpath, T, K, tmask)


% keyboard;


thematrix = zeros(K,K);
diags = zeros(1,K);

% assuming that each T == 1 subject:
e=0;
for i_T = 1:numel(T)
    b=e+1;
    e = b + T(i_T) - 1;
    this_string = sprintf('calculating transitions: %d - %d; %d samples.', b, e, e-b+1);
    %if i_T > 1
        % fprintf('\b', 1, prev_string)
    %else
        fprintf('%s\n', this_string)
        prev_string = this_string;
    %end
    
        
    try
       this_vpath = vpath(b:e);
    catch
        keyboard;
    end
  
    
    % thematrix = zeros(K, K);
    
    for i_K = 1:numel(K)
        for j_K = 1:numel(K)
            for i_vpath = 2:numel(this_vpath)
                this_state = this_vpath(i_vpath);
                prev_state = this_vpath(i_vpath-1);
                
                thematrix(prev_state, this_state) = thematrix(prev_state, this_state) + 1;
            end
        end
    end
end


all_transitions_from = sum(thematrix, 2);
for i_row = 1:K
    thematrix(i_row, :) = thematrix(i_row, :) / all_transitions_from(i_row);
end

my_diagonal = diag(thematrix);
my_pruned_matrix = thematrix;
for i_K = 1:K
    my_pruned_matrix(i_K, i_K) = nan;
end
matrices = my_pruned_matrix;
diags = my_diagonal;


