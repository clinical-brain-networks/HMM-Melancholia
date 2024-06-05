

function post_hmm_2_make_state_figures(which_file, varargin)


if numel(varargin) > 0
    network_descriptions = varargin{1};
else
    network_descriptions = {'','','','','','','','','','','',''};
    
end




fbgcolor = [1 1 1]; % [0.95 0.95 0.95 ];
BOUTLINE = 1; % brain-outline
VISIBILITY = 'on';
CLOSE_FIGS = false;


% which_file = 'rest_run/HMMrun_K10_rep_1.mat';
[my_p, my_f, my_e] = fileparts(which_file);

output_dir = [my_p filesep 'figures'];
if ~exist(output_dir,'dir')
    try
    mkdir(output_dir);
    catch
        keyboard;
    end
end
output_filename = [my_p filesep 'figures' filesep 'Fig1_v2_states_' my_f '.jpg'];






% phil's SPM is probably fine.
addpath /mnt/data/johanv/Phil/scripts/spm12/
% these are stil fine

roi_base = '/mnt/data/johanv/Phil/scripts/craft_brainfig_hmm/Positive_Network_III_smoothed.nii';

% make the stupid figure - 1
% for printing the network -- in COLORS...
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



% first let's get all of the ROIs pls - in T1 space (2009c version)
% cmap_activities=fireice;

% cmap_activities = [cmap_activities(33:end,:); cmap_activities(1:32,:)];
% trim it a little bit:
% cmap_activities(32-7:33+7,:)=[];
% cmap_activities(1:9,:)=[];
% cmap_activities(end-8:end,:) = [];
% trim white:
% cmap_activities = [cmap_activities(33:end,:); cmap_activities(1:32,:)];
% cmap_activities=colormap(jet);
load('bw.mat');
cmap_activities = bw;


roi_vols = spm_vol(roi_base);
roi_files={};
roi_vol={};
for i=1:numel(roi_vols)
    %roi_files{i} = [roi_base 's' sprintf('%.2d',i) '-to-2009c-in-2009c-dims.nii'];

    roi_data{i} = spm_read_vols(roi_vols(i));
    
    
    % write 0's as NaN's:
    % roi_data{i}(roi_data{i}==0)=NaN;
    % roi_vol{i} = spm_write_vol(roi_vol{i},roi_data{i});
    % V = spm_write_vol(V,Y)
    
    
end
t1_file =  '/mnt/data/johanv/Phil/scripts/craft_brainfig_hmm/MNI152_T1_1mm_brain.nii';
t1_vol = spm_vol(t1_file);
t1_data = spm_read_vols(t1_vol); %-- it's a 193x229x193


t1_outline_file = '/mnt/data/johanv/Phil/scripts/craft_brainfig_hmm/MNI152_T1_1mm_brain-dil1.nii';
t1_outline_vol = spm_vol(t1_outline_file);
t1_outline_data = spm_read_vols(t1_outline_vol); %-- it's a 193x229x193


% check Summary_measures_1
a=load(which_file);

% get the coeffs of all the states
% this yields a roi-by-state matrix. each col is a representation of which
% roi's/networks are represented in that state.
% we go through states later and make nice T1 figs.
allw=[];
try
for i=1:a.K
    allw=[allw a.hmm.state(i).W.Mu_W']; % beware of shoddy concatenation!!
end
NSTATES = a.K;
% keyboard;
end




% spm fmri
% nake 12
%coronal_slice=100;
sl=100;
thr_fig=0.25;
clims=[-0.5 0.5];





fh=[];  % to store figure handles. We basically make 10
% separate figures (one for each state). In each figure
% there are 2 + 14 (T1 + Outline + 14 brain network)
% images drawn over each other with transparency set
% properly. Then one big figure is made, and for each
% of the 10 figures, the contents are copied into the
% big figure.

% B2POS = [0 0 0.2 1]+ [0.3 0 0 0]; % this becomes... [0 0.1 0.9*0.2 0.9];
% B1POS = [0.2 0 0.23 1]+ [0.3 0 0 0];
% BARPOS = [0.45 0 0.53 0.1]
TEXTCOLOR = 'w';
% TEXTCOLOR = [0.95 0.95 0.95];
TEXTPOS = [0 0.2 0.1 0.7] + [0.0375 0.5 0 -0.5];
% TEXTPOS = [0.45 0.2 0.15 0.7];

% B2POS = [0 0.1 0.9*0.2 0.9] + [0.15 0 0 0];
% B1POS = [0.2 0.1 0.9*0.2 0.9] + [0.15 0 0 0];
BARPOS = [0 0 0.53 0.095]  + [0.0375 0 0 0];
% TEXTCOLOR = 'w';
% TEXTPOS =


sl_and_pos = {...
    49, [0 0.1 0.9*0.2 0.9] + [0.15 0 0 0];...
    72, [0.2 0.1 0.9*0.2 0.9] + [0.15 0 0 0];...
    };

slice_x = [1:size(t1_outline_data, 1)];
slice_y = [(1+9):(size(t1_outline_data, 2) - 9)];
% keyboard;

for state=1:NSTATES
    
    fh(state)=figure('position',[400,400,1000,200],'color',fbgcolor,'visible',VISIBILITY);
    set(fh(state),'units','normalized')
    
    for i_sl_and_pos = 1:size(sl_and_pos,1)
        this_sl = sl_and_pos{i_sl_and_pos,1};
        this_pos = sl_and_pos{i_sl_and_pos,2};
        
        % brain outline?
        if BOUTLINE
            % brain outline - to make the networks more pronounced if
            % they're white and on the edge of the brain.
            ah1axoutl=axes('parent',fh(state),'position',this_pos);
            set(ah1axoutl,  'view', [180 -90]);
            t1data=0.5 * squeeze(t1_outline_data(slice_x,slice_y,this_sl))';
            im1=imagesc(ah1axoutl,1-t1data);
            set(ah1axoutl,'clim',[0 1]);
            % set(im1,'alphadata',t1data>0);
            axis xy
            colormap(ah1axoutl,'bone')
            set(ah1axoutl,'visible','off');
        end

        % make an axial section?
        ah1ax=axes('parent',fh(state),'position',this_pos);
        set(ah1ax, 'view', [180 -90]);
        t1data=squeeze(t1_data(slice_x,slice_y,this_sl))';
        im1=imagesc(ah1ax,t1data);
        set(im1,'alphadata',t1data>0);
        axis xy
        colormap(ah1ax,'bone')
        set(ah1ax,'visible','off');

        % do for each ROI the following:
        % translate_mat = 1:14; %[4 9 8 1 11 14 6 5 2 10 7 12 3 13];
        phs=[];
        for i=1:numel(roi_vols)

            this_w = allw(i, state);
            ah2ax=axes('parent',fh(state),'position',this_pos);

            roidata=squeeze(roi_data{i}(slice_x,slice_y,this_sl))';
            im2ax=imagesc(ah2ax,(roidata>thr_fig)*this_w);
            axis xy
            colormap(ah2ax,cmap_activities);
            set(ah2ax,'clim',clims);
            set(im2ax,'alphadata',roidata>thr_fig);
            set(ah2ax,'visible','off');
            colormap(ah2ax);
            % colormap(im2,cmap_activities);
            set(ah2ax, 'view', [180 -90]);


        end
    
    end
    
     
    ah_state = axes('parent', fh(state), 'position', BARPOS);
    % set(ah_state, 'view', [180 -90]);
    set(ah_state,'xlim', [0 numel(roi_vols)], 'ylim', [0, 1]);
    colormap(ah_state, cmap_activities);
    set(ah_state,'visible','off');
    set(ah_state,'clim', clims);
    
    
    phs=[];
    translate_mat = 1:numel(roi_vols);
    for i=1:numel(roi_vols)

        this_w = allw(i, state);
        patch_coords_x = [0 1 1 0] + (translate_mat(i)-1);
        patch_coords_y = [0 0 1 1];
        phs(translate_mat(i)) = patch(ah_state,patch_coords_x,patch_coords_y,this_w);

    end

    
%     keyboard;
    
    ah_text = axes('parent',fh(state),'position',TEXTPOS);
    set(ah_text,'visible','off');
    ph=patch(ah_text,[0 1 1 0],[0 0 1 1],TEXTCOLOR);
    set(ph,'edgecolor',TEXTCOLOR);
    thxl=get(gca,'xlim');
    thyl=get(gca,'ylim');
    th3=text(thxl(1)-diff(get(gca,'xlim'))*0.0005,thyl(2),sprintf(' State %d',state),'fontweight','bold','verticalalignment','top');
    th3=text(thxl(1)+diff(get(gca,'xlim'))*0.0005,thyl(2),sprintf(' State %d',state),'fontweight','bold','verticalalignment','top');
    th3=text(thxl(1),thyl(2)-diff(get(gca,'ylim'))*0.0005,sprintf(' State %d',state),'fontweight','bold','verticalalignment','top');
    th3=text(thxl(1),thyl(2)+diff(get(gca,'xlim'))*0.0005,sprintf(' State %d',state),'fontweight','bold','verticalalignment','top');
    th3=text(thxl(1),thyl(2)+diff(get(gca,'xlim'))*0.0005,sprintf(' State %d',state),'fontweight','bold','verticalalignment','top','color',cmat(state,:)/255);
    %th=text(thxl(1),thyl(2),sprintf(' Network %d',state),'fontweight','bold','Color',cmat(state,:)/255);
    
    
    
    
    % keyboard;
    if numel(varargin) > 0 %  == 12
        th2=text(5.5, 0.5, sprintf('\n\n%s',varargin{1}{state}),'verticalalignment','top');
        set(th2,'interpreter','none');
    end
    
    
end




% keyboard;




% the magic -- make the Big Plot - subdivide into 5-by-2 figs
% but we need an fh array of fig handles.
bigf=figure('color',fbgcolor,'visible',VISIBILITY);
set(bigf,'position',[50 50 800 800])
fcols=2;
frows=6;
offx=0; %% how much space at bottom of figure?
offy=0.1;

spacing=0.001;
factorx = (1-offx-(fcols-1)*spacing)/fcols;
factory = (1-offy-(frows-1)*spacing)/frows;

% keyboard;


for i=1:NSTATES %:-1:1
    [fj, fi] = ind2sub([frows fcols],i);
    fj=frows+1-fj;
    
    objs = findobj(fh(i),'type','axes');
    
    for i_o = numel(objs):-1:1
        o=objs(i_o);
        old_position = get(o, 'position');
        old_xp = old_position(1);
        old_yp = old_position(2);
        old_xs = old_position(3);
        old_ys = old_position(4);
        
        % calculate new position according to my notes:
        new_xp = old_xp * factorx + (fi-1) * (factorx + spacing) + offx;
        new_yp = old_yp * factory + (fj-1) * (factory + spacing) + offy;
        new_xs = old_xs * factorx;
        new_ys = old_ys * factory;
        
        new_position = [new_xp new_yp new_xs new_ys];
        
        % then copy that obj into the new figure:
        newo = copyobj(o, bigf);
        set(newo,'position',new_position);
    end
end


% keyboard;

cb_ax = axes('parent',bigf);
set(cb_ax,'position',[0.2 0.02 0.6 0.06]);
cb=colorbar('North');
set(cb,'Limits',clims);
set(cb_ax,'clim',clims);
colormap(cb_ax,cmap_activities);
set(cb_ax,'visible','off');

annot_ax = axes('parent',bigf);
set(annot_ax,'position',[0.15 0.03 0.7 0.05],'visible','off');
set(annot_ax,'xlim',[0 14]);
set(annot_ax,'ylim',[0 2]);
% 
% DESCRIPTORS = {'dDMN','PRE','vDMN','ASN','PSN','LECN','RECN','BGN','AUD','pVIS','hVIS','SMN','VSN','LAN'};
% phs=[];
% for i=1:14
%     patch_coords_x = [0 1 1 0] + (i-1);
%     patch_coords_y = [0 0 1 1]+1;
%     ph=patch(annot_ax,patch_coords_x,patch_coords_y,0);
%     phs(end+1) = ph;
%     set(ph,'facealpha',0);
%     pdiff=0.06;
%     switch i
%         case {3, 5, 8, 13}
%             th=text(i-0.5-pdiff/2,1.5,DESCRIPTORS{i});
%             % disp('doingit');
%         case {4, 6, 9, 14}
%             th=text(i-0.5+pdiff/2,1.5,DESCRIPTORS{i});
%             % disp('doingit2');
%         otherwise
%             th=text(i-0.5-pdiff/2,1.5,DESCRIPTORS{i});
%             % disp('doingit3');
%     end
%     
%     
%     
%     set(th,'fontsize',7,'horizontalalignment','center');
% end

% 
% % DMN line
% annots_labels={'DMN','SAL','EXEC','SENS'};
% annots_begins = [0, 3, 5, 8];
% annts_ends = [3, 5, 8, 13];
% 
% for iannot=1:4
%     b=annots_begins(iannot);
%     e=annts_ends(iannot);
%     t=annots_labels{iannot};
%     
%     line([b, b, e, e] ,[0.75, 0.5, 0.5, 0.75],'color','k');
%     
%     th=text(mean([b, e]), 0.4,t);
%     set(th,'horizontalalignment','center','verticalalignment','top','fontsize',7,'color','k');
%     % line([b b],[1 2],'linewidth',2,'color','k');
%     % line([e e],[1 2],'linewidth',2.5,'color','k');
%     % line([e e],[1 2.05],'linewidth',0.5,'color','w');
%     % line([b e e b b],[1 1 2 2 1],'linewidth',2,'color','k');
% end
% 
% set(phs(3),'XData', get(phs(3),'XData') - [0 pdiff pdiff 0]')
% set(phs(4),'XData', get(phs(4),'XData') + [pdiff 0 0 pdiff]')
% set(phs(5),'XData', get(phs(5),'XData') - [0 pdiff pdiff 0]')
% set(phs(6),'XData', get(phs(6),'XData') + [pdiff 0 0 pdiff]')
% set(phs(8),'XData', get(phs(8),'XData') - [0 pdiff pdiff 0]')
% set(phs(9),'XData', get(phs(9),'XData') + [pdiff 0 0 pdiff]')
% set(phs(13),'XData', get(phs(13),'XData') - [0 pdiff pdiff 0]')
% set(phs(14),'XData', get(phs(14),'XData') + [pdiff 0 0 pdiff]')
% 



set(cb_ax,'visible','off');

if CLOSE_FIGS
    close(fh);
end
cb.Position=cb.Position + [0 0 0 0.009];

cb.Position=cb.Position - [0 0.01 0 0];


a_spacingx = 0.1;
a_spacingy = 0.9;

a1_x = [cb.Position(1) + a_spacingx * cb.Position(3), cb.Position(1)];
a1_y = cb.Position(2) + (1+a_spacingy) * cb.Position(4) * [1 1];

a2_x = [cb.Position(1) + (1-a_spacingx) * cb.Position(3), cb.Position(1) + cb.Position(3)];
a2_y = a1_y;

a1_h = annotation('textarrow',a1_x, a1_y,'String','  Low');
a2_h = annotation('textarrow',a2_x, a2_y,'String','High  ');


anno_ax = axes('parent',bigf,'position',[0 0 1 1],'visible','off');
a3_h = text(cb.Position(1) + 0.5*cb.Position(3),a2_y(1),'Mean','HorizontalAlignment','Center');
set(a3_h,'Parent',anno_ax);


set(cb,'Ticks', [-0.5 0 0.5]);
set(cb,'TickLabelInterpreter','None');
set(cb,'TickLabels',{'-50%','0','+50%'});

set(cb,'FontSize',10);
set(cb,'FontWeight','bold');


fh_overview_states=bigf;


set(fh_overview_states,'paperunits','centimeters');
set(fh_overview_states,'papersize',1.25 * [22 16]);
set(fh_overview_states,'paperposition',1.25 * [0 0 22 16]);


% saveas(fh_overview_states,['~/' preprocessing_dir '-' analysis_dir '-' num2str(SUMMARY_MEASURES_FILE) '.jpg']);
% saveas(fh_overview_states,[preprocessing_dir '-' analysis_dir '-' num2str(SUMMARY_MEASURES_FILE) '.jpg']);



% fprintf('Done making Fig1 for %s, %s, file: %d\n',preprocessing_dir, analysis_dir, SUMMARY_MEASURES_FILE);

% the following can be uncommented toautomatize that one of the
% figures is copied automatically to a 'figures' directory.
% fsource=[preprocessing_dir '-' analysis_dir '-' num2str(SUMMARY_MEASURES_FILE) '.jpg'];

% ftarget_file = 'f1_state_10.jpg';

% ftarget1=[base_cwd '/../figures/' ftarget_file];
% copyfile(fsource,ftarget1);

% print('-djpeg','-r300', output_filename);

% if CLOSE_FIGS
%     close(fh_overview_states);
% end
% 
% % we will ALSO ... make text files with the values for that figure!
% Regions={'roi-1_ASN'; 'roi-2_AUD'; 'roi-3_BGN'; 'roi-4_dDMN'; 'roi-5_hVIS'; 'roi-6_LAN'; 'roi-7_LECN'; 'roi-8_PSN'; 'roi-9_PRE'; 'roi-10_pVIS'; 'roi-11_RECN'; 'roi-12_SMN'; 'roi-13_vDMN'; 'roi-14_VSN'};
% 
% % warning: Extremely Ugly code needed to make the Matlab Tables
% % play nice.
% for i=1:NSTATES
%     s=sprintf('State_%d_Weights = allw(:,%d);\n',i,i);
%     eval(s);
% end



% keyboard;
saveas(fh_overview_states, output_filename);

% keyboard;

if CLOSE_FIGS
    close(fh_overview_states);
end

fprintf('Done making Fig1: %s\n', output_filename);

