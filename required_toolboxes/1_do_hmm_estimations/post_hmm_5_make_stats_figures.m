
function out = post_hmm_5_make_stats_figures(file_in, do_this, v, varargin)

% 'FO',  is_ie_melancholic, pos_mask
label_text_left = 'NON-MEL';
label_text_right = 'MEL';

if numel(varargin)>0
    label_text_left = varargin{1}{1};
    label_text_right = varargin{1}{2};
end
do_the_distributions = 1;
if numel(varargin)>1
   do_the_distributions = varargin{2}; 
end

tmask=[];
my_extension = '';
if numel(varargin)>2
   tmask = varargin{3}; 
   my_extension = '_masked';
end

added_ext = '';
if numel(varargin)>3
   added_ext = varargin{4}; 
   my_extension = [my_extension '_' added_ext];
end

% keyboard;

[f_p, f_n, f_e] = fileparts(file_in);

% keyboard;

out_mel = post_hmm_4_investigate_state_stats(file_in, find(v==1), tmask);
out_nonmel = post_hmm_4_investigate_state_stats(file_in, find(v==0), tmask);


% do_this = 'FO';
disable_all_but_first = 0;
switch do_this
    
    case 'FO'
        
        plot_values_left = out_nonmel.FO(:,:);
        plot_values_right = out_mel.FO(:,:);
        my_y_label_text = 'Fractional Occupancy [0-1]';

        
    case 'lifetimes'
        
        plot_values_left = cellfun(@(x) mean(x), out_nonmel.lifetimes);
        plot_values_right = cellfun(@(x) mean(x), out_mel.lifetimes);
        my_y_label_text = 'State Dwell Times [sec.]';

        
    case 'intervals'
        
        plot_values_left = cellfun(@(x) mean(x), out_nonmel.intervals);
        plot_values_right = cellfun(@(x) mean(x), out_mel.intervals);
        my_y_label_text = 'State Interval Times [sec.]';
        
    case 'switchingRate'
        
        plot_values_left = 60*repmat(out_nonmel.switchingRate(:,:),1, out_mel.K);
        plot_values_right = 60*repmat(out_mel.switchingRate(:,:), 1, out_mel.K);
        my_y_label_text = 'Switching Rate [switches per minute]';
        disable_all_but_first = 1;
        
        
        
        
        
        
        
end

if ~strcmp(do_this, 'switchingRate')


    % keyboard;

    diagnosis = [repmat({'nonmel'},size(plot_values_left,1),1);repmat({'mel'},size(plot_values_right,1),1)];

    m = [plot_values_left; plot_values_right];


    my_labels = {}; 
    for i_my_labels = 1:size(plot_values_left, 2)
        my_labels{end+1} = sprintf('state%d', i_my_labels);
    end

    my_tmp_variablenames = {};
    for i_my_labels = 1:size(plot_values_left, 2)
        my_tmp_variablenames{end+1} = sprintf('state%d', i_my_labels);
    end
    my_tmp_variablenames = ['diagnosis', my_tmp_variablenames];
    my_s = '';
    my_s = [my_s sprintf('t=table(diagnosis, ')];
    for i_my_labels=1:size(plot_values_left, 2)
        my_s = [my_s sprintf('m(:, %d), ',  i_my_labels)];
    end
    my_s = [my_s sprintf('\''VariableNames\'', my_tmp_variablenames)')];

    eval(my_s);


    % t=table(diagnosis, m(:,1), 'VariableNames', my_tmp_variablenames);
    Meas = table(my_labels', 'VariableNames',{'States'});

    try
    rm = fitrm(t, sprintf('state1-state%d~diagnosis', size(plot_values_left, 2)) ,'WithinDesign', Meas);

        [manovatbl, A, C, D] = manova(rm);
        ranovatbl = ranova(rm);


        % disp(this_measure);
        save(sprintf('ranovatbl_%s%s.mat', do_this, my_extension), 'ranovatbl');
        save(sprintf('manovatbl_%s%s.mat', do_this, my_extension), 'manovatbl');
    catch
        warning('The Stats Failed!');
    end
    
end


disp(file_in);
disp(plot_values_left);

disp(file_in);
disp(plot_values_right);



this_big_title = sprintf('%s between Depressed and Melancholic Depressed', do_this);
if numel(varargin)>0
    this_big_title = sprintf('%s between %s and %s', do_this, varargin{1}{1}, varargin{1}{2});
end


color_lines_left = [0 0 0];
color_lines_right = [0.35 0.35 0.35];
% label_text_left = 'NON-MEL';
% label_text_right = 'MEL';
write_figsize_centimeter = [35 25];
% keyboard;
if strcmp(f_p, '')
    this_path = pwd;
else
    this_path = f_p;
end

% UGLY
if ~exist(sprintf('%s/figures', this_path),'dir')
    mkdir(sprintf('%s/figures', this_path))
end
to_save_filename = sprintf('%s/figures/compare_%s_%s%s', this_path, do_this, f_n, my_extension);
% my_y_label_text = 'Fractional Occupancy [0-1]';




cb=[[230, 25, 75];...
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
    [0, 0, 0]]/255;


% addpath('../../scripts/RainCloudPlots/tutorial_matlab/')


nstates =out_mel.K;

nrows = 2;
ncols = nstates/2;
xspacing = 0.001;
yspacing = 0.05;
x_pad_low = 0.075;
y_pad_low = 0.1;
x_pad_high = 0.025;
y_pad_high = 0.025;



% make a matrix of positions...
axes_positions = {};
do_axes_labels = [];
do_labels = [];

for i_nrows = 1:nrows
    for i_ncols = 1:ncols
        x_size = (1 - x_pad_low - x_pad_high - xspacing*(ncols-1)) / ncols;
        y_size = (1 - y_pad_low - y_pad_high - yspacing*(nrows-1)) / nrows;
        
        x_pos = x_pad_low + (i_ncols-1)*(xspacing + x_size);
        y_pos = 1 - y_pad_high - i_nrows* (y_size + yspacing) + yspacing;
        
        axes_positions{end+1} = [x_pos y_pos x_size y_size];
        
        if i_ncols == 1
            do_axes_labels(end+1) = 1;
        else
            do_axes_labels(end+1) = 0;
        end
        
        if i_nrows == nrows
            do_labels(end+1) = 1;
        else
            do_labels(end+1) = 0;
        end
        
        
        
        
    end
end




fh=figure;

xlims = [prctile(sort([plot_values_left(:); plot_values_right(:)]), 1), prctile(sort([plot_values_left(:); plot_values_right(:)]), 99)];
x_step = diff(xlims)/4;
add_extra_percent_to_xlims = 5;
xlims = [(xlims(1) - diff(xlims)/100*add_extra_percent_to_xlims) (xlims(2) + diff(xlims)/100*add_extra_percent_to_xlims*2)];
if disable_all_but_first
    xlims = [xlims(1) xlims(2)*2];
end

if disable_all_but_first
    nstates = 1;
end

stats_th_fdr = [];
pvalues_fdr = [];
for i_state = 1:nstates
    %this_plot = subplot(nstates/3, 3, i_state);
    % keyboard;
    
    axes_left = axes('parent',fh,'position', axes_positions{i_state});
    h1 = raincloud_plot(plot_values_left(:, i_state), 'box_on', 1, 'color', cb(i_state,:), 'alpha', 0.4,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 0);
    set(h1{1}, 'Ydata', get(h1{1},'YData') / 4);
    set(get(h1{1},'Baseline'),'linestyle','none')
    
    if ~do_the_distributions
        set(h1{1},'visible','off');
    end
    
    % move it:
    move_it_by = get(h1{3},'position')*[0 1 0 0]' * -1 + diff(get(get(h1{1},'parent'),'ylim'))/100*5; % 1 percent?;
    old_basevalue = get(get(h1{1},'Baseline'), 'BaseValue');
    old_ydata = get(h1{1},'YData');
    
    % keyboard;
    set(get(h1{1},'Baseline'),'BaseValue', old_basevalue + move_it_by);
    set(h1{1},'YData', old_ydata + move_it_by);
    set(h1{2}, 'YData', get(h1{2},'YData') + move_it_by);
    set(h1{3}, 'position', get(h1{3},'position') + [0 move_it_by 0 0]);
    set(h1{4}, 'YData', get(h1{4},'YData') + move_it_by);
    set(h1{5}, 'YData', get(h1{5},'YData') + move_it_by);
    set(h1{6}, 'YData', get(h1{6},'YData') + move_it_by);
    
    set(h1{1},'EdgeColor',color_lines_left);
    set(h1{2},'MarkerFaceColor',color_lines_left);
    set(h1{3},'EdgeColor', color_lines_left);
    set(h1{4},'Color', color_lines_left);
    set(h1{5},'Color', color_lines_left);
    set(h1{6},'Color', color_lines_left);
    
    if ~do_the_distributions
        % set(h1{2},'MarkerFaceColor','r');
        set(h1{3},'visible','off');
        set(h1{4},'visible','off');
        set(h1{5},'visible','off');
        set(h1{6},'visible','off');
        
        
        the_y_data = get(h1{2},'ydata');
        the_x_data = get(h1{2},'xdata');
        for i_the_data = 1:numel(the_x_data)
            text(the_x_data(i_the_data), the_y_data(i_the_data), sprintf('%d', i_the_data));
        end

        
        
    end
    
    
    
    view([-90, 90]);
    set(get(h1{1},'parent'),'visible','off');
    set(get(h1{1},'parent'), 'xlim', xlims);
    
    axes_right = axes('parent',fh,'position', axes_positions{i_state});
    
    h2 = raincloud_plot(plot_values_right(:, i_state), 'box_on', 1, 'color', cb(i_state,:), 'alpha', 0.4,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 0);
    set(h2{1}, 'Ydata', get(h2{1},'YData') / 4);
    set(get(h2{1},'Baseline'),'linestyle','none')
    
    if ~do_the_distributions
        set(h2{1},'visible','off');
    end
    
    % move it:
    move_it_by = get(h2{3},'position')*[0 1 0 0]' * -1 + diff(get(get(h2{1},'parent'),'ylim'))/100*5; % 1 percent?;
    old_basevalue = get(get(h2{1},'Baseline'), 'BaseValue');
    old_ydata = get(h2{1},'YData');
    
    set(get(h2{1},'Baseline'),'BaseValue', old_basevalue + move_it_by);
    set(h2{1},'YData', old_ydata + move_it_by);
    set(h2{2}, 'YData', get(h2{2},'YData') + move_it_by);
    set(h2{3}, 'position', get(h2{3},'position') + [0 move_it_by 0 0]);
    set(h2{4}, 'YData', get(h2{4},'YData') + move_it_by);
    set(h2{5}, 'YData', get(h2{5},'YData') + move_it_by);
    set(h2{6}, 'YData', get(h2{6},'YData') + move_it_by);
    
    
    set(h2{1},'EdgeColor',color_lines_right);
    set(h2{2},'MarkerFaceColor',color_lines_right);
    set(h2{3},'EdgeColor', color_lines_right);
    set(h2{4},'Color', color_lines_right);
    set(h2{5},'Color', color_lines_right);
    set(h2{6},'Color', color_lines_right);

    if ~do_the_distributions
        % set(h1{2},'MarkerFaceColor','r');
        set(h2{3},'visible','off');
        set(h2{4},'visible','off');
        set(h2{5},'visible','off');
        set(h2{6},'visible','off');
        the_y_data = get(h2{2},'ydata');
        the_x_data = get(h2{2},'xdata');
        for i_the_data = 1:numel(the_x_data)
            text(the_x_data(i_the_data), the_y_data(i_the_data), sprintf('%d', i_the_data));
        end

        
        
    end
   %  keyboard;

    
    view([90, -90]);
    set(get(h2{1},'parent'),'visible','off');
    set(get(h2{1},'parent'), 'xlim', xlims);
    % set(axes_right,'xlim',get(axes_left,'xlim'));
    % set(axes_right,'ylim',get(axes_left,'ylim'));
    
    
    anno_ax = axes('parent', fh, 'position', axes_positions{i_state});
    my_xlims = get(axes_left,'xlim');
    my_ylims = get(axes_left,'ylim');
    set(anno_ax,'xlim',my_xlims);
    set(anno_ax,'ylim',my_ylims);
    view([90, -90]);
    % if disable_all_but_first
        th=text(my_xlims(2) - diff(my_xlims)/100*7.5, mean(my_ylims(:)), sprintf('S %d', i_state));
    % else
    %     th=text(my_xlims(2) - diff(my_xlims)/100*7.5, mean(my_ylims(:)), 'all');
    % end
    
    set(th,'horizontalalignment','center','verticalalignment','top');
    
    
    [my_h, my_p, my_ci, my_stats] = ttest2(plot_values_left(:, i_state), plot_values_right(:, i_state));
    pvalues_fdr(end+1) = my_p;
    stats_th = text(my_xlims(2) - diff(my_xlims)/100*15, mean(my_ylims(:)), sprintf('p = %2.2f', my_p));
    
    set(stats_th,'horizontalalignment','center','verticalalignment','top');
    set(stats_th, 'fontsize', get(stats_th,'fontsize') - 2);
    stats_th_fdr(end+1) = text(my_xlims(2) - diff(my_xlims)/100*20, mean(my_ylims(:)), '');
    set(stats_th_fdr(end),'horizontalalignment','center','verticalalignment','top');
    set(stats_th_fdr(end), 'fontsize', get(stats_th,'fontsize') - 2);
    
    
    
    lh = line(anno_ax, 'XData', [my_xlims(1) my_xlims(1) + diff(my_xlims)*0.95], 'YData', my_ylims(1)*[1 1]);
    tick_values = [xlims(1):x_step:xlims(2)];
    tick_values(end) = [];
    for i_tick_values = 1:numel(tick_values)
        tlh = line(anno_ax, 'XData', tick_values(i_tick_values)*[1 1], 'YData', my_ylims(1)+ [0 diff(my_ylims)/100*5]);
        tth = text(tick_values(i_tick_values), my_ylims(1) - diff(my_ylims)/100*0.5, sprintf('%.2f', tick_values(i_tick_values)));
        set(tth,'horizontalalignment','right','verticalalignment', 'middle');
        set(tth,'fontsize', get(tth,'fontsize') - 0);
        
        if ~do_axes_labels(i_state)
            set(tth,'visible','off');
        end
        
    end
    
    if do_labels(i_state)
        lthl=text(my_xlims(1) - diff(my_xlims)/100*2, my_ylims(1) + diff(my_ylims)/100*25, label_text_left);
        set(lthl,'horizontalalignment','center','verticalalignment','top');
        set(lthl,'color', color_lines_left);
        
        lthr=text(my_xlims(1) - diff(my_xlims)/100*2, my_ylims(1) + diff(my_ylims)/100*75, label_text_right);
        set(lthr,'horizontalalignment','center','verticalalignment','top');
        set(lthr,'color', color_lines_right);
    end
    
    
    if do_axes_labels(i_state)
        lthlyl=text(mean(my_xlims),  my_ylims(1) - diff(my_ylims)/100*35, my_y_label_text);
        set(lthlyl,'horizontalalignment','center','verticalalignment','top');
        set(lthlyl,'rotation', 90);
    end
    
    
    
    
    set(anno_ax,'visible','off');
    
    
end


if nstates > 1
    % keyboard;
    fdr_corr_pvalues = mafdr(pvalues_fdr', 'BHFDR', true);
    
    % [h, crit_p, adj_ci_cvrg, fdr_corr_pvalues]=fdr_bh(pvalues_fdr,0.05);
    
    for i=1:numel(fdr_corr_pvalues)
        try
            set(stats_th_fdr(i),'string',sprintf('FDR p: %2.2f', fdr_corr_pvalues(i)));
        catch
            keyboard;
        end
    end
end


set(fh,'color','w');


anno_ax = axes('parent',fh,'position',[0 0 1 1],'visible','off');
th=text(0.5, 0.999, this_big_title);
set(th,'horizontalalignment','center','verticalalignment','top');
set(th,'fontsize', 16);

set(fh,'paperunits','centimeters');
set(fh,'papersize',write_figsize_centimeter);
set(fh,'paperposition',[0 0 write_figsize_centimeter]);
print('-djpeg','-r250',to_save_filename);
print('-dpdf',to_save_filename);

out=fh;
