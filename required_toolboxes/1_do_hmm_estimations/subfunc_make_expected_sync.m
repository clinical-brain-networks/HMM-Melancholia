function [lower_lim, upper_lim, thr_mean] = subfunc_make_expected_sync(NSUB, transitions, NTIMEPOINTS, NPERMS)

NSTATES = size(transitions, 1);
% keyboard;



target_intervals = cumsum(transitions,2);


thresholds_perms = zeros(NPERMS, NTIMEPOINTS);
for i_perm = 1:NPERMS

    
    
    mat_timeseries_subs = [];
    for i_sub = 1:NSUB
        
        
        % randomize the initial state:
        initial_state = randi(NSTATES);
        timeseries = zeros(1,NTIMEPOINTS);
        timeseries(1) = initial_state;
        
        for i_timepoints = 2:NTIMEPOINTS
            % this, truly is, the blackest of magicks.
            % to attain this infernal wisdom, several bunnies were
            % sacrificed and their entrails studiously investigated
            timeseries(i_timepoints) = find(rand > [0 target_intervals(timeseries(i_timepoints-1),:)],1, 'last');
        end
        mat_timeseries_subs(end+1,:) = timeseries;
        
    end
    
    
    
    time_lag=7;
    
    % keyboard;
    % what's the threshold?
    % STHR=floor(nsub/100*20);
    % OFFSET=floor(nsub/100*15);
    
    
    m=mat_timeseries_subs;
    % check out, then, the threshold, pls:
    outs=[];
    outc=[];
    for im=1:size(m,2)
        
        b=im-floor(time_lag/2);
        e=im+floor(time_lag/2);
        
        
        if b < 1 || e > size(m,2)
            out(im)=NaN;
            
        else
            
            
            all_vals=m(:,b:e);
            
            uv=unique(all_vals);
            ct=[];
            for iuv=1:numel(uv)
                ct(end+1)=sum(sum(all_vals==uv(iuv),2)>0);
            end
            
            um=sortrows([uv ct'],2,'descend');
            try
            outs(end+1,:)=um(1:2,1);
            outc(end+1,:)=um(1:2,2);
            catch
                outs(end+1,:)=um(1);
                outc(end+1,:)=um(2);
            end
            
        end
    end
    THR_NSUB = [NaN; NaN; nan; outc(:, 1) ;NaN; NaN; NaN];
    THR_NSUB_RAND = [NaN; NaN; nan; randi(NSTATES, numel(outc(:, 1)), 1)  ;NaN; NaN; NaN];
    
    % disp(i_perm);
    
    thresholds_perms(i_perm, :) = THR_NSUB;
end

% keyboard;
thr_mean = nanmedian(nanmean(thresholds_perms));
lower_lim = prctile(thresholds_perms(:), 5);
upper_lim = prctile(thresholds_perms(:), 95);
 % thr_std = nanmedian(nanstd(thresholds_perms,[],1));
    