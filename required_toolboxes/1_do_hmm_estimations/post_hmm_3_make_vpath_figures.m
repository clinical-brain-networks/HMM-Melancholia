% clear all; close all; clc;

function post_hmm_3_make_vpath_figures(file_in, varargin)

if numel(varargin) > 0
    sort_index = varargin{1};
    sort_it = 1; % TRUE!
else
    sort_it = 0;
end

% RATIO_UPPER_LOWER = [0.6 0.4];
MY_RESIZE_FACTOR = 0.08;

append_this = '';
tmask=[];

make_nullmodel=0;
if numel(varargin) > 1
    make_nullmodel = varargin{2};
end


if numel(varargin) > 2
    append_this = ['_' varargin{3}];
end

if numel(varargin) > 3
    append_this = [append_this '_masked'];
    tmask = varargin{4};
  
end

added_ext='';
if numel(varargin) > 4
    added_ext = varargin{5};
    % keyboard;
    append_this = [append_this '_' added_ext];
  
end

% keyboard






% which_file = 'rest_run/HMMrun_K10_rep_1.mat';
[my_p, my_f, my_e] = fileparts(file_in);

output_dir = [my_p filesep 'figures'];
disp(output_dir);
if ~exist(output_dir,'dir')
    % keyboard;
    mkdir(output_dir);
end

output_filename = [my_p filesep 'figures' filesep 'Fig2_vpath_' my_f append_this '.jpg'];



try
a=load(file_in);
TR = 1/a.options.Fs;
NSTATES = a.K;
nvols = a.T;
nsub = numel(a.T);
if sort_it
    nsub = numel(sort_index);
end
vpath = a.vpath;
catch
    keyboard;
end
CLOSE_FIGS = 0;

% apply the nice colormap for maximum contrast:
cmat=[[230, 25, 75];...
    [60, 180, 75];...
    [255, 225, 25];...
    [0, 130, 200];...
    [245, 130, 48];...
    [70, 240, 240];...
    [240, 50, 230];...
    [250, 190, 190];...
    [0, 128, 128];...
    [230, 190, 255];...
    [170, 110, 40];...
    [255, 250, 200];...
    [128, 0, 0];...
    [170, 255, 195];...
    [0, 0, 128];...
    [128, 128, 128];...
    [255, 255, 255];...
    [0, 0, 0]];

% the parula colormap is used for the binary masks (i.e. story annotations)
parulamap=colormap('parula');

cmat=cmat(1:NSTATES,:);
cmat=cmat./255;


% container for all data -- PER SUBJECT
sd=struct(); % subject data

my_cond_names = {'rest'};
new_statemat = struct(); new_statemat(:)=[]; % this is for counting/keeping track of the states.



all_anno_axes=[];
fh=figure('color','w','visible','off');
set(fh,'position',[50         120        1700         800]);
if ~CLOSE_FIGS
    set(fh,'visible','on');
end



% the state path is one big vector of everything - here we
% subdivide by subject/session again. The b and e are running
% indices that change number upon each iteration, so that b:e
% 'slices' the correct partition of vpath.
mat={};
e_vpath=0;
b_vpath=1;
for invols=1:numel(nvols)
    nvol=nvols(invols);
    e_vpath = e_vpath+  nvol;
    tot_el=numel(b_vpath:e_vpath);
    mat{end+1} = reshape(vpath(b_vpath:e_vpath),tot_el/1,1);
    b_vpath=e_vpath+1;
end




% now, let's apply the mask:

if numel(tmask) > 0
    new_mat = mat;
    make_lines_here = {};
    for i_mat = 1:numel(new_mat)
        my_keep_these = intersect(1:numel(new_mat{i_mat}), find(tmask));
        my_logical_array = zeros(1,numel(new_mat{i_mat}));
        my_logical_array(my_keep_these) = 1;
        some_data = my_find_logical_clusters(my_logical_array);
        % make_lines_here{end+1} = [some_data(1:end-1).e]-0.5;
        make_lines_here{end+1} = find(diff(my_keep_these) ~= 1)+0.5;
        new_mat{i_mat} = new_mat{i_mat}(my_keep_these);
    end

    mat=new_mat;
end
% keyboard;




%
%
%
% We can also (!!) try to figure out, from the transition matrix (which we
% have), what the av. expected threshold is...
%
%

% we have the mat, and we have the transitions...
% keyboard;

% keyboard;
if sort_it
    % nsubs to check; transition matrix; how many timepoints we expect (not
    % all have the same # of timepoints... so median); Number of
    % permutations to do.
    
    % keyboard;
    my_new_out = post_hmm_4_investigate_state_stats(file_in, sort_index);
    
    [sync_lower_lim, sync_upper_lim, sync_thr_mean] = subfunc_make_expected_sync(numel(sort_index), my_new_out.av_trans_matrix_full, median(cellfun(@(x) numel(x), mat)), 1000);
else
    [sync_lower_lim, sync_upper_lim, sync_thr_mean] = subfunc_make_expected_sync(numel(mat), a.hmm.P, median(cellfun(@(x) numel(x), mat)), 1000);
end



    


%
%
% make the figure.
%
%

nsub1=nsub;
nsub2=numel(a.T) - nsub;
use_this_max_nsub = max([nsub1, nsub2]);


fcols=1; %numel(mat);
frows=use_this_max_nsub;

spacingx=0.012;
spacingy=0.002;
bigspacingxl = 0.10;
bigspacingyl = 0.05;

bigspacingxu = 0.10;
bigspacingyu = 0.20 + MY_RESIZE_FACTOR;


factorx = (1-bigspacingxl-bigspacingxu-(fcols-1)*spacingx)/fcols;
factory = (1-bigspacingyl-bigspacingyu-(frows-1)*spacingy)/frows;
ah=[];
im=[];


% so now -- some magic to scale the factors, so as to have it
% even nicer.

extra_factorx=[1]; %1 1 1 1];

total_scans=[];
for i_gr=1:numel(mat)
    total_scans(end+1) = size(mat{i_gr}, 1);
end
for i_gr=1% :numel(mat)
    extra_factorx(i_gr) = 1/(1/numel(mat) / (total_scans(i_gr)/sum(total_scans)));
end



%
%         % in which sequence to present figures:
%         if numel(scan_names) == 4
%             new_order=[3 1 4 2];
%         else
%             new_order=[1 2];
%         end
new_order = [1];

for i_gr2=1:numel(new_order)
    i_gr = new_order(i_gr2);
    
    x_offset = 0;
    if i_gr2 > 1
        for ii_gr=1:(i_gr2-1)
            x_offset = x_offset + factorx*extra_factorx(new_order(ii_gr)) + spacingx;
        end
    end
    
    for i=1:use_this_max_nsub % *numel(nvols)
        
        
        [fj, fi] = ind2sub([frows fcols],i+numel(a.T)*(i_gr-1));
        
        
        %
        old_xp = 0; % old_position(1);
        old_yp = 0; % old_position(2);
        old_xs = 1; % old_position(3);
        old_ys = 1; % old_position(4);
        
        new_xp = old_xp * factorx*extra_factorx(i_gr) + x_offset + bigspacingxl;
        new_yp = 1-(old_yp * factory + (fj) * (factory + spacingy) + bigspacingyl);
        new_xs = old_xs * factorx*extra_factorx(i_gr);
        new_ys = old_ys * factory;
        
        try
        new_position = [new_xp new_yp new_xs new_ys];
        
        
       %  keyboard;
        ah(end+1) = axes('parent',fh,'position',new_position);
        catch
            keyboard;
        end
        %
      
            

        if sort_it
            % this_sub = sort_index(i);

            % this_sub = i; % chosen_subs(nsub-i+1);
            this_sub_i = i - (use_this_max_nsub - nsub); % chosen_subs(nsub-i+1);
            if this_sub_i < 1
                set(ah(end), 'visible','off');
                draw_axis = false;
            else
                draw_axis = true;
                try
                this_sub = sort_index(this_sub_i);
                catch
                    keyboard;
                end
            end
        else
            this_sub_i = i;
            draw_axis=true;
        end

               
        if draw_axis
                

            % disp(this_sub);


            % v=mat(:,nsub-(i-1))';
            try
            v=mat{this_sub};% (:,this_sub)';
            catch
                keyboard;
            end


            % save v into big data struct:
            % sd.(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(subs_to_use(i))]).(regexprep(analyses{i_analyses,1}{i_gr},'-','_')).hmmdata=v;

            NSTATES = a.K;
            % convert state indices into 0-1 matrix
            % save into big data struct:
            statemat=zeros(numel(v),NSTATES);
            % statemat(:)=NaN;
            for itstate=1:NSTATES
                statemat(v==itstate,itstate) = 1;
                %if sum(v==itstate)==0whos v
                %    statemat(:,itstate)=NaN;
                %end
            end
            % sd.(regexprep(preprocessing_dir,'-','_')).(regexprep(analysis,'-','_')).(['sm' num2str(SUMMARY_MEASURES_FILE)]).(['s' num2str(subs_to_use(i))]).(regexprep(analyses{i_analyses,1}{i_gr},'-','_')).hmmmat=statemat;



            if numel(v) < 1460
               %  keyboard;
            end

            N_ELEMENTS_TO_PLOT = max(cellfun(@(x) numel(x), mat));

            im(this_sub_i) = imagesc(1:numel(v), 1, v');

            % make some lines...
            % keyboard;
            if numel(tmask)>0
                these_make_lines_here = make_lines_here{this_sub_i};
                pa = get(im(this_sub_i),'parent');
                for i_these_make_lines_here = 1:numel(these_make_lines_here)
                    this_line_x = these_make_lines_here(i_these_make_lines_here);
                     % keyboard;
                    my_ylims = get(pa,'ylim');
                    % annotations are pretty useless - they require [0-1] interval,
                    % bleh
                    % annotation('line', this_line_x*[1 1], my_ylims + [-1.5 1.5]*diff(my_ylims),'color','k','linewidth',3);
                    %  annotation('line', this_line_x*[1 1], my_ylims + [-1.5 1.5]*diff(my_ylims),'color','k','linewidth',3);
                     line(this_line_x*[1 1], get(pa,'ylim'),'color','k','linewidth', 3);
                    line(this_line_x*[1 1], get(pa,'ylim'),'color','r','linewidth',0.5);

                end
            end




            old_xlims = get(get(im(this_sub_i),'parent'),'xlim');
            set(get(im(this_sub_i),'parent'),'xlim', [old_xlims(1), N_ELEMENTS_TO_PLOT]);

            set(ah(end),'visible','off');
            set(ah(end),'colormap',cmat);
            set(ah(end),'clim',[0.5 NSTATES+0.5]);
            % set(ah(end),'colormap',cmat);

            if this_sub_i==1
                % keyboard;
                yl=get(ah(end),'ylim');
                xl=get(ah(end),'xlim');
                tha=text(mean(xl),yl(1),sprintf('vpath figure %s', append_this));
                set(tha,'fontweight','bold');
                set(tha,'horizontalalignment','center','verticalalignment','bottom');
                set(tha,'parent',ah(end));
                set(tha,'interpreter','none');
            end

            if i_gr2==1
                th=text(ah(end),-10,0,sprintf('P%d - %d',this_sub,''));
                set(th,'horizontalalignment','right','units','normalized');
                set(th,'position',[-0.001 0.5]);
            end

        end
        
        
    end
    
    
    % keyboard;
    ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-factory-factory+0.005 new_xs 0.010]);
    %set(ah(end),'ytick',[],'yticklabel',[],'ycolor','w');
    set(ah(end),'xlim',[0 numel(v)*TR]);
    time_v = TR*(1:numel(v));
    ph=plot(time_v, ones(length(time_v)));
    set(ph,'visible','on','color','w');
    set(ah(end),'yticklabel',{},'ytick',[]);
    set(ah(end),'xlim',[0 numel(v)*TR]);
    xlh=xlabel('time (s)');
    % set(xlh,'horizontalalignment','right');
    box off
    
    
    
    %
    %
    %      EXTRA FIG -- STATES (the most common ones)
    %
    %
    %
    
    
    
    time_lag=7;
    
    % keyboard;
    % what's the threshold?
    STHR=floor(nsub/100*20);
    OFFSET=floor(nsub/100*15);
    
    
    % pre-check:
    
    mat_lengths = cellfun(@(x) numel(x), mat);
    fix_this = find(mat_lengths ~= median(mat_lengths));
    for i_fix_this = 1:numel(fix_this)
            %     keyboard;

        this_fix_this = fix_this(i_fix_this);
        
        old_v = mat{i_fix_this};
        new_v = mat{1} * 0;
        mat{this_fix_this} = new_v;
        mat{this_fix_this}(1:numel(old_v)) = old_v;
    end
    
    
    try
    % keyboard;
        
    m = [mat{:}]';
    if sort_it
        m = [mat{sort_index}]';
    end
    
    catch
        this_mat = mat{end};
        prev_mat = mat{end-1};
    end
    
    % calculate the lower & upper lims, for this one:
    % this is a tricky call - you want everything to put in, not just the
    % subset of m!!
    if sort_it
        if make_nullmodel
            [lower_lims, upper_lims, this_average, count_matrix_all, observed_count] = make_null_for_vpath_average([mat{:}]', sort_index, time_lag);


            to_save.lower_lims = lower_lims;
            to_save.upper_lims = upper_lims;
            to_save.count_matrix_all = count_matrix_all;
            to_save.observed_count = observed_count;


            output_filename_vpath_examination = [my_p filesep 'figures' filesep 'Fig2_vpath_data_' my_f append_this '.mat'];
            save(output_filename_vpath_examination,'to_save');
        end
    end
    
    
    outs=[];
    outc=[];
    for im=1:size(m,2)
        
        b=im-floor(time_lag/2);
        e=im+floor(time_lag/2);
        
        
        if b < 1 || e > size(m,2)
            try
                out(im)=NaN;
            catch
                keyboard;
            end
            
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
    
    selection=outc(:,1)>STHR;
    
    vout=outs(:,1);
    
    
    vout(selection==0)=NaN;
    
    % keyboard;
    ah(end+1) = axes('parent',fh,'position',[new_xp bigspacingyu-factory-0.05-0.05-0.03+0.02-MY_RESIZE_FACTOR new_xs 0.05+MY_RESIZE_FACTOR-0.01]);
    
    
    
    inpatch=0;
    vout = [nan(floor(time_lag/2),1);vout];
    for i=2:numel(vout)
        
        
        if ~inpatch
            line_px=[(i-0.5)*TR];
            line_py=[OFFSET];
        end
        
        if ~isnan(vout(i))
            
            inpatch=1;
            
            if vout(i) ~= vout(i-1) && ~isnan(vout(i-1))
                line_px(end+1) = (i-0.5)*TR;
                line_py(end+1) = OFFSET;
                line_px(end+1) = line_px(1);
                line_py(end+1) = OFFSET;
                draw_it;
                % keyboard;
                
                line_px=[(i+[-0.5])*TR];
                
                line_py=[OFFSET];
            end
            
            
            line_px(end+1)=(i+[-0.5])*TR;
            line_px(end+1)=(i+[ 0.5])*TR;
            
            line_py(end+1)=outc(i-floor(time_lag/2),1)*[1 ];
            line_py(end+1)=outc(i-floor(time_lag/2),1)*[1 ];
            
            % keyboard;
            
            val=vout(i);
            confidence=outc(i-floor(time_lag/2),1);
            
            % ph_back=patch(TR*(i+[-0.5 0.5 0.5 -0.5]), [0 0 1 1]*confidence,[1 1 1]);
            ph=patch(TR*(i+[-0.5 0.5 0.5 -0.5]), [0 0 1 1]*confidence,cmat(val,:));
            set(ph,'linestyle','none');
            set(gca,'xlim',[0 size(m,2)*TR]);
            set(gca,'ylim',[OFFSET-0.2 size(m,1)]);
            set(gca,'xtick',[],'xticklabel',[]);
            if i_gr2 > 1
                set(gca,'ytick',[],'yticklabel',[]);
            end
            
            if i_gr2 == 1
                set(gca,'ytick',[OFFSET nsub],'yticklabel',{[sprintf('%d',round(OFFSET/nsub*100)) '%'],'100%'});
            end
            
        end
        
        
        if isnan(vout(i)) && numel(line_px)>1 && inpatch && i ~=numel(vout)
            line_px(end+1) = (i-0.5)*TR;
            line_py(end+1) = OFFSET;
            line_px(end+1) = line_px(1);
            line_py(end+1) = OFFSET;
            draw_it;
            inpatch=0;
            % keyboard;
        end
        
        if i==numel(vout) && inpatch
            line_px(end+1) = (i+0.5)*TR;
            line_py(end+1) = OFFSET;
            line_px(end+1) = line_px(1);
            line_py(end+1) = OFFSET;
            draw_it;
            inpatch=0;
        end
        
        
        time_which_consistency = [ (1:(numel(vout) + floor(time_lag/2)))'*TR [[zeros(floor(time_lag/2), 1); outc(:,1); zeros(floor(time_lag/2),1)] [vout; zeros(floor(time_lag/2), 1)]]];
        new_statemat(1).(my_cond_names{i_gr}) = time_which_consistency;
        
    end
    
    if i_gr2 == 1
        % keyboard;
    end
    
    
end

set(fh,'visible','on');

set(gca,'nextplot','add');
if sort_it
    if make_nullmodel
        conf_x_values = (1:size(m,2))' * TR;
        to_get_rid_of = floor(time_lag/2);

        conf_x_values(1:to_get_rid_of) = [];
        conf_x_values((end-to_get_rid_of+1):end) = [];
        conf_y_values_l = lower_lims;
        conf_y_values_u = upper_lims;

        make_path_grey = patch([conf_x_values; conf_x_values(end:-1:1)], [conf_y_values_l conf_y_values_u(end:-1:1)],[0.2 0.2 0.2]);
        set(make_path_grey,'facecolor', [0.5 0.5 0.5], 'linestyle','-','edgecolor','k', 'facealpha', 0.5);
        uistack(make_path_grey,'bottom');
    end
end

% keyboard;

line(get(gca,'xlim'), [1 1] * sync_upper_lim, 'color', 'k');
% my_thr_mean + 1.96*my_thr_std


% keyboard

if sort_it
    if make_nullmodel
        underwhelming_here = observed_count' < conf_y_values_l;
        % exceeding_here = observed_count' > conf_y_values_u;
        % keyboard;
        exceeding_here = (observed_count' > conf_y_values_u)  & (observed_count' > sync_upper_lim);
        for p = find(underwhelming_here)
            % let us only do the exceeding model(s):
            % text(conf_x_values(p)-floor(time_lag/2)*TR, get(gca,'ylim')*[1 0]', '↓','verticalalignment','top', 'horizontalalignment','center');
        end

        for p = find(exceeding_here)
            % text(conf_x_values(p)-floor(time_lag/2)*TR, get(gca,'ylim')*[0 1]', '↑','verticalalignment','bottom');
            % text(conf_x_values(p)-floor(time_lag/2)*TR, max(observed_count) + 0.1 * (max(observed_count) - min(observed_count)), '↑','verticalalignment','bottom');
            text(conf_x_values(p), max(observed_count) + 0.1 * (max(observed_count) - min(observed_count)), '↑','horizontalalignment','center','verticalalignment','bottom');
        end
    end
end


ah(end+1) = axes('parent',fh,'position',[new_xp + new_xs+0.02, bigspacingyu, bigspacingxu, 1-bigspacingyl-bigspacingyu]);

set(ah(end),'visible','off');
cb=colorbar('West');
set(cb,'Limits',[0.5 0.5+NSTATES]);
set(ah(end),'clim',[0.5 0.5+NSTATES]);
cb.Position=cb.Position + [0 0 0.0018 0];
cb.Ticks=[1:NSTATES];
cb.Label.String='HMM States';
%set(cb,'colormap',cmat);



colormap(cmat);
set(all_anno_axes,'colormap',parulamap);

calc_size=0;
for imat=1:1 %numel(mat)
    calc_size = calc_size+size(mat{imat},1);
end
calc_size= calc_size/50;
calc_size= calc_size + (numel(mat)-1) *0.8;


set(fh,'paperunits','centimeters');
set(fh,'papersize',[10+ calc_size, 20]*1.2);
set(fh,'paperposition',[0 0 10+ calc_size 20]*1.2);

% 
% % fig_fname=['Fig2_vpath_' preprocessing_dir '-' analysis '-rep-' num2str(SUMMARY_MEASURES_FILE) '.jpg'];
% 
% output_filename = [scriptpath filesep '..' filesep 'figures' filesep 'Fig2_' preprocessing_dir '-' analysis_dir '-' num2str(SUMMARY_MEASURES_FILE) '.jpg'];
% 
% 

% 
% 

% 
% 
% % Matlab tables are definitely not pandas dataframes:
% Participant = {};
% for i=1:numel(subs_to_use)
%     Participant{end+1} = sprintf('Participant_%d',i);
% end
% Participant=Participant';
% Participant{end+1} = 'Consistency';
% Participant{end+1} = 'MostConsistentState';


% now we have made several variables; so we need to call Table with
% that; i need to say myTable, otherwise it might interfere with
% another Table already in the namespace. Furthermore, we need to
% use eval, since Table names columns according to names in the
% workspace.

% Making Tables in Matlab. It's so simple and elegant, a true joy
% to behold, and an example to strive towards to. So much better
% than Pandas dataframes!


%         Circumvent the following error:
%         Error using writetable (line 124)
% The data block starting at cell 'A1' exceeds the sheet boundaries by 0 row(s) and 1235 column(s).
%
% Error in Make_Fig2_state_path_figure (line 648)
%         writetable(myTable, output_filename_values)
% 
% 
% 
% for i=1:numel(my_cond_names)
%     s='';
%     s=[s sprintf('%s = mat{%d}\''\n',my_cond_names{i}, i)];
%     s=[s sprintf('%s=[%s; round(new_statemat.%s(:,2)/numel(subs_to_use)*100)\'']\n',my_cond_names{i},my_cond_names{i},my_cond_names{i})];
%     s=[s sprintf('%s=[%s; new_statemat.%s(:,3)\'']\n',my_cond_names{i},my_cond_names{i},my_cond_names{i})];
%     
%     s=[s sprintf('%s = num2str(%s)', my_cond_names{i}, my_cond_names{i})];
%     
%     eval(s);
%     
% end


% 
% 
% 
% s='';
% s=[s sprintf('myTable = table(')];
% s=[s sprintf('Participant,')];
% my_conditions = fieldnames(new_statemat);
% for i=1:numel(my_conditions)
%     s=[s sprintf('%s,',my_conditions{i})];
%     
% end
% s(end)=[]; % remove comma
% s=[s sprintf(');')];
% eval(s); % this will make the table.
% % fsource=fig_fname;

% 
% output_filename_values = [scriptpath filesep '..' filesep 'figures' filesep 'Fig2_' preprocessing_dir '-' analysis_dir '-' num2str(SUMMARY_MEASURES_FILE) '_values.xls'];
% 
% myTable;
% writetable(myTable, output_filename_values)
% 
saveas(fh, output_filename);
print('-djpeg','-r900', regexprep(output_filename, '.jpg', '_highres.jpg'));


if CLOSE_FIGS
    close(fh);
end







