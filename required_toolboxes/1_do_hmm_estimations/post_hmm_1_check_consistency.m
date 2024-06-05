% # check occupancy of states script


function post_hmm_1_check_consistency(indir)


% get all of the hmm files in that dir:
d=dir([indir filesep 'HMMrun_*.mat']);



s = sprintf('%s\n', indir);



for i_d = 1:numel(d)
    
    
    fname=sprintf([d(i_d).folder filesep d(i_d).name]);
    a=load(fname);
    
    
    totscans = sum(a.T);
    nsub = numel(a.T);
    
    if i_d == 1
        s = sprintf('file\t');
        for i_state = 1:a.K
            s = sprintf('%socc_st_%d\t', s, i_state);
        end
        
    s = sprintf('%stotscans\t', s);
    s = sprintf('%snsub\n', s);

    
        
    end


    s = sprintf('%s%s\t', s,  d(i_d).name);
    
    
    
    si=[];for i=1:a.K;si(i) = sum(a.vpath==i);end
    s=[s sprintf('%d\t', si)];
    s=sprintf('%s%d\t', s, totscans);
    s=sprintf('%s%d\t', s, nsub);
    
    s(end) = [];
    s = sprintf('%s\n', s);
    
end

fprintf(s)

% fid=fopen('vpath-diagostic.txt','w+');
% fprintf(fid,s);
% fclose(fid);

% disp('file written to vpath-diagnostic.txt -- check THAT!')





