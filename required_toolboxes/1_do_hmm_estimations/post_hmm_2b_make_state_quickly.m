function post_hmm_2b_make_state_quickly(which_file, atlaslabelsfile)



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
output_filename = [my_p filesep 'figures' filesep 'Fig1_nobrain_states_' my_f '.jpg'];


if ~exist([my_p filesep 'figures'],'dir')
    mkdir([my_p filesep 'figures']);
end

CLOSE_FIGS=0;


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


% keyboard;
labels = {};
fid=fopen(atlaslabelsfile);
while ~feof(fid)
    labels{end+1} = fgetl(fid);
end

labels = cellfun(@(x) regexprep(x,'_',' '), labels, 'uniformoutput', false);


state_legends = {};
for i_state = 1:a.K
    state_legends{end+1} = sprintf('S%3d', i_state);
end


fh=figure;

plot(allw, 'linewidth', 1.5);

xlim([0 numel(labels)+1]);
set(gca,'xtick', 1:numel(labels));
set(gca,'xticklabels', labels);
set(gca,'xticklabelrotation', 60);

legend(state_legends,'location', 'northeastoutside');





set(fh,'paperunits','centimeters');
set(fh,'papersize',1 * [22 16]);
set(fh,'paperposition',1 * [0 0 22 16]);


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
saveas(fh, output_filename);

% keyboard;

if CLOSE_FIGS
    close(fh);
end

fprintf('Done making Fig1: %s\n', output_filename);

