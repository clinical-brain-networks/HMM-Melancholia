function out = post_hmm_7_prepare_NBS_stats(file_in, v, varargin)

% file_in = 'movie/HMMrun_K10_rep_50.mat';


include_diag = false;
if numel(varargin) > 0
    include_diag = varargin{1};
end

my_extension = '';
tmask = [];
if numel(varargin)>1
    tmask = varargin{2};
    my_extension = '_masked';
end

normalize_mats = false;
if numel(varargin)>2
    normalize_mats= varargin{3};
end

added_ext = '';
if numel(varargin)>3
    added_ext= varargin{4};
    my_extension = [my_extension '_' added_ext];
end





ismel = find(v==1);
isnonmel = find(v==0);

% COMP_A = 'mov1a';
% COMP_B = 'resta';

% load valid_inferences_all.mat
% RUN=valid_inferences(1); % you might change this depending on which HMM inference you wish to check out.


% ANALYSIS='all';

% run=sprintf('run%d',RUN);
% subs_to_use_real = [2     3     7     8     9    10    11    12    14    16    17    18    19    20];


out_mel = post_hmm_4_investigate_state_stats(file_in, ismel, tmask);
out_nonmel = post_hmm_4_investigate_state_stats(file_in, isnonmel, tmask);

K = out_mel.K;


out_mat_a = [];
for i=1:numel(ismel)
    out_mat_a = cat(3, out_mat_a, out_mel.trans_matrices{i});
    diag_elelents = out_mel.trans_diags{i};
    
    if include_diag
        for j=1:size(out_mat_a, 1)
            out_mat_a(j,j,end) = diag_elelents(j); % we can zero-fill? Or?
        end
    end
end


out_mat_b = [];
for i=1:numel(isnonmel)
    out_mat_b = cat(3, out_mat_b, out_nonmel.trans_matrices{i});
    diag_elelents = out_nonmel.trans_diags{i};
    if include_diag
        for j=1:size(out_mat_b, 1)
            out_mat_b(j,j,end) = diag_elelents(j);
        end
    end
end





this_nbs_dir = sprintf('for_nbs_K%d%s',K, my_extension);
if isfolder(this_nbs_dir)
    % keyboard;
    rmdir(this_nbs_dir,'s');    
end
mkdir(this_nbs_dir);

this_nbs_matdir = sprintf('%s/matrices',this_nbs_dir);
if isfolder(this_nbs_matdir)
    % keyboard;
    rmdir(this_nbs_matdir,'s');
end
mkdir(this_nbs_matdir);

% keyboard;
X = [[ones(numel(isnonmel), 1); zeros(numel(ismel), 1)] [zeros(numel(isnonmel), 1); ones(numel(ismel), 1)]];
% we could add in their age, but I am not going to do thatnow.

% labels:
labels = '';
for i=1:K
    labels = [labels sprintf('State_%d\n', i)];
end

% we don't really have ROI locations, eh

% save the labels:
to_save_labels_filename = sprintf('%s/labels.txt', this_nbs_dir);
fid=fopen(to_save_labels_filename,'w+');
fprintf(fid, labels);
fclose(fid);

% save the text(s):
save(sprintf('%s/X.txt',this_nbs_dir), 'X', '-ascii');

% save the matrices; first the NONMEL, then the MEL.
% keyboard;
all_mat = cat(3, out_mat_a, out_mat_b);

% keyboard;

new_all_mat= [];
for i=1:size(all_mat,3)
    new_all_mat(:,:,end+1) = all_mat(:,:,i) ./ nansum(all_mat(:,:,i),2);
end


for i=1:size(all_mat,3)
    if normalize_mats
        m = new_all_mat(:,:,i);
    else
        m = all_mat(:,:,i);
    end
    save(sprintf('%s/mat_%03d.txt', this_nbs_matdir, i), 'm','-ascii');
end


this_nbs_matdir_bigmatdir = sprintf('%s/bigmat', this_nbs_matdir);
if isfolder(this_nbs_matdir_bigmatdir)
    rmdir(this_nbs_matdir_bigmatdir,'s');
end
mkdir(this_nbs_matdir_bigmatdir);

m=all_mat;
save(sprintf('%s/all_matmat.mat', this_nbs_matdir_bigmatdir), 'm');

% keyboard;
% make the NODE COORDS (needed for the GUI):
delta_angle = 2*pi/K;
all_angles = 0:delta_angle:2*pi;all_angles(end) = [];
% keyboard;
coords = [];
begin_z_coord = -20;
% end_z_coord = 20;
for i_all_angles = 1:numel(all_angles)
    coords(i_all_angles, :) = 30*[sin(all_angles(i_all_angles)), cos(all_angles(i_all_angles)) ,(begin_z_coord + (i_all_angles-1)/(K-1)*40)/30];
end


save(sprintf('%s/coords.txt', this_nbs_dir), 'coords', '-ascii');


out = this_nbs_dir;





