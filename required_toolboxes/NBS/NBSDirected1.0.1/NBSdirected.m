function [n_cnt,con_mat,pval]=NBSDirected(varargin)
%NBSDirected Computes directed network components among edges that survive a primary 
%test statistic threshold. Assigns a corrected p-value to each component 
%using permuted data that has been supplied as part of the structure STATS. 
%
%   [N_CNT,CON_MAT,PVAL]=NBSDirected(STATS) operates on the network and
%   associated statistical data in the structure STATS. 
%
%   [...]=NBSDirected(STATS,H) writes out progress to a listbox with handle H. 
%   If writing progress to H fails, progress is written to screen instead. 
%
%   [...]=NBSDirected(STATS,H,GLM) GLM is mandatory if the permuted data has 
%   not been precomputed, which is the case when STATS.test_stat is empty. 
%   In this situation, dNBSglm is repeatedly called to compute a test
%   statistic for each permutation. Slower than precomputation, but saves 
%   memory. 
%
%   A STATS structure contains the following fields:
%       STATS.thresh:     Primary test statistic threshold    
%       STATS.alpha:      Corrected significance (user specified), network 
%                         components not satisfying alpha signficance are 
%                         not reported
%       STATS.N:          Number of nodes in network
%       STATS.test_stat:  K+1 x J array of test statistics. The first 
%                         row is the oberved test statistic. The remaining 
%                         rows are samples of the test statistic under the 
%                         null hypothesis. Each column corresponds to a 
%                         seperate edge. K is number of permutations. J is
%                         the number of edges. The test statistics can be
%                         computed with dNBSglm. Columns are mapped to
%                         edges such that column i=1:J corresponds to the
%                         edge with index ind_upper(i), where ind_upper are
%                         the indexes of the upper trianguler elements. 
%                         ind_upper = find(triu(ones(N,N),1)); 
%       STATS.size        'extent' | 'intensity' 
%                         Measure used to assess size of a network 
%                         component  
%                          
%   Outputs:
%       N_CNT:            Number of network components satisfying alpha 
%                         significance
%       CON_MAT:          1 x N_CNT cell array of adjacency matrices. Each
%                         cell holds a N x N upper-triangular adjacency 
%                         matrix specifying a network component satisfying 
%                         alpha signifcance
%       PVAL:             1 x N_CNT array of corrected p-values 
%   
%   Remarks:
%       If no network components satisfy alpha significance, CON_MAT and
%       PVAL are returned empty and N_CNT = 0.  
%
%       STATS.test_stat is empty if the number of permutations is too large
%       to precompute. See Limit parameter in dNBSrun for details. In this 
%       situation, dNBSglm is repeatedly called (J times) to compute test 
%       statistics for each permutation. 
%
%   azalesky@unimelb.edu.au

%Number of most recent permutations to display in listbox

%Number of most recent permutations to display in listbox
nDisp=5;

STATS=varargin{1}; 
if nargin==3
    %Handle to listbox
    H=varargin{2}; 
elseif nargin==4
    %Handle to GLM
    H=varargin{2};
    GLM=varargin{3};
    NBS = varargin{4};
end

%Is BGL available?
bgl=0;
if exist('components','file')==2
    %Use components.m provided by MatlabBGL, otherwise use get_components_dNBS.m
    bgl=1;
end

[currentFolder,name,ext] = fileparts(which('dNBS.m'));
%currentFolder = pwd;
GroupFolder= strrep(datestr(datetime), ' ', '_');
GroupFolder= strrep(GroupFolder,':','_');



%Connectivity matrices (regions x regions x subjects)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for n = 1:numel(STATS.Folder)%pathFolders)
cd(NBS.FolderConnMat);
searchStr = {'*.txt'};


pathparts = strsplit(NBS.FolderConnMat,'\');
strFolder = char(fullfile(GroupFolder, pathparts(end)));%%%%%%%5
%Test


%Test END
[filepath2,name2,ext2] = fileparts(NBS.FolderConnMatFile);

%Test
if ext2 == '.mat'
    C = readUI(NBS.FolderConnMatFile);
     DIMS.observations=size(C,3);
else 
    for f=1:size(searchStr,2)    
        fileNames = dir(searchStr{f});

        for files=1:size(fileNames,1)
            disp(fileNames(files,1).name);
            data = importdata(fileNames(files,1).name);
           %Summary.files = data';     % wenn files nicht alle exakt lang sind (e.g. resting state)
           %X=data.data([50:150],:)'
            C(:,:,files) = data(:,:);  %wenn files alle gleich lang sind (e.g. exportierte Segmente aus Analyzer)
        end
        DIMS.observations=size(C,3);;
        %save(['Summary_',outputNames{f}],'Summary');
    end
end
cd (currentFolder);
%GLM.X = load(pathGroup);


%Number of matrices

DIMS.preditcors=size(GLM.X,2);%Number of Groups
DIMS.nodes=size(C,2);


cd(currentFolder);

 


 %Contrast

 for k = STATS.thresh(1):STATS.thresh(2):STATS.thresh(3)
     %NO NEED TO CHANGE ANYTHING BEYOND THIS POINT
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Threshold = k;
     bgl=0; 
     nDisp=5; 
     nbs.NBS.n =0;
     nbs.NBS.test_stat=[];
     nbs.NBS.pval = [];
     nbs.NBS.con_mat = [];
     %Number of nodes
     

     N=size(C,1);
     ind_uplo=union(find(triu(ones(N,N),1)),find(tril(ones(N,N),-1)));

     GLM.y=zeros(size(GLM.X,1),length(ind_uplo));
     for i=1:size(C,3)
         tmp=C(:,:,i);
         GLM.y(i,:)=tmp(ind_uplo);
     end

    %Precompute test stat
    STATS.test_stat=dNBSglm(GLM); 

    %Number of edges
    J=length(ind_uplo); 

    %Determine whether test statistics have been precomputed and determine
    %index of edges exceeding the primary threshold
    if ~isempty(STATS.test_stat)
        %Precomputed test statistics
        ind=ind_uplo(STATS.test_stat(1,:)>Threshold); 
        %Number of permutations
        K=size(STATS.test_stat,1)-1; 
    else
        %Never get to this case
    end

    %Size of a component measured using extent or intensity? 
    Intensity=0;
    if strcmp(STATS.size,'Intensity')    
        %If size measure using intensity, create an N x N matrix cotaining the 
        %test statistic for each edge minus the test statistic threshold
        %(primary threshold)
        Intensity=1; 
        %Compute a test statistic matrix
        test_stat_mat=zeros(N,N); 
        if ~isempty(STATS.test_stat)
            %Precomputed
            test_stat_mat(ind_uplo)=STATS.test_stat(1,:)-Threshold;
    %       test_stat_mat=(test_stat_mat+test_stat_mat');
        else
            %Never reach this case. 
        end
    end

    adj=spalloc(N,N,length(ind));
    adj(ind)=1; 
    %Only consider components comprising more than one node, equivalent to at
    %least one edge
    if bgl==1
        [a,sz]=components((adj+adj')/2); 
    else
        [a,sz]=get_components_dNBS((adj+adj')/2); 
    end
    ind_sz=find(sz>1);
    sz_links=zeros(1,length(ind_sz));
    max_sz=0; 
    for i=1:length(ind_sz)
        nodes=find(ind_sz(i)==a);
        if Intensity
            %Measure size as intensity
            sz_links(i)=sum(sum(adj(nodes,nodes).*test_stat_mat(nodes,nodes))); %/2;
        else
            %Measure size as extent
            sz_links(i)=sum(sum(adj(nodes,nodes))); %/2;
        end
        adj(nodes,nodes)=adj(nodes,nodes)*(i+1);
        if max_sz<sz_links(i)
            max_sz=sz_links(i);
        end
    end

    %Subtract one to remove edges not part of a component
    %Although one is also subtracted from edges comprising a component, this is 
    %compensated by the (i+1) above
    adj(~~adj)=adj(~~adj)-1;

    %Repeat above for each permutation
    %Empirical null distribution of maximum component size
    null_dist=zeros(K,1); 
    str1='| Permutation | Max Size | Max Size | Lowest  |';
    str2='|             |  Random  |  Actual  | p-value |';
    try tmp=get(H,'string'); set(H,'string',[{str1};{str2};tmp]); drawnow;
    catch;  fprintf([str1,'\n',str2,'\n']); end 
    p_approx=0;
    %Store what is al ready displayed in the listbox
    try pre_str=get(H,'string'); catch; end
    new_str={};
    %First row of test_stat is the observed test statistics, so start at the
    %second row
    for i=2:K+1
        if ~isempty(STATS.test_stat)
            %Precomputed test statistics 
            ind=ind_uplo(STATS.test_stat(i,:)>Threshold); 
        else

        end
        if Intensity 
            %Compute a test statistic matrix
            test_stat_mat=zeros(N,N); 
            if ~isempty(STATS.test_stat)
                test_stat_mat(ind_uplo)=STATS.test_stat(i,:)-Threshold;
            else

            end    
        end
        adj_perm=spalloc(N,N,length(ind));
        adj_perm(ind)=1;
        if bgl==1
            [a,sz]=components((adj_perm+adj_perm')/2); 
        else
            [a,sz]=get_components_dNBS((adj_perm+adj_perm')/2); 
        end
        ind_sz=find(sz>1);
        max_sz_perm=0; 
        for j=1:length(ind_sz)
            nodes=find(ind_sz(j)==a);
            if Intensity
                tmp=sum(sum(adj_perm(nodes,nodes).*test_stat_mat(nodes,nodes))); %/2;
            else
                tmp=sum(sum(adj_perm(nodes,nodes))); %/2;
            end
            if tmp>max_sz_perm
                max_sz_perm=full(tmp);
            end   
        end
        null_dist(i-1)=max_sz_perm; 
        if max_sz_perm>=max_sz
            p_approx=p_approx+1;
        end
       % str=sprintf('|   %5d/%5d |     %4d |     %4d |   %0.3f |',...
       %v1.1.2 Changed to %6.0f to %6.1f to allow fractional component sizes
       %that arise when component size is measured with intensity. 
       str=sprintf('| %5d/%5d |   %6.1f |   %6.1f |  %0.10f |',...
                i-1,K,max_sz_perm,max_sz,p_approx/(i-1));
            %Display no mare than nDisp most recent permutations
            new_str=[str,{new_str{1:min(nDisp,length(new_str))}}]';
            try set(H,'string',[new_str;pre_str]); drawnow; 
                catch;  fprintf([str,'\n']); end 
    end
    str1='| Permutation | Max Size | Max Size | Lowest  |';
    str2='|             |  Random  |  Actual  | p-value |';
    try tmp=get(H,'string'); set(H,'string',[{str1};{str2};tmp]); drawnow;
    catch;  fprintf([str1,'\n',str2,'\n']); end 
    
    
            % keyboard;


    test_stat=zeros(2,DIMS.nodes*(DIMS.nodes-1)/2);
    test_stat=STATS.test_stat(1,:);
    if isempty(STATS.test_stat)
        K=nbs.GLM.perms;
        %Temporarily set to 1 to save computation
        nbs.GLM.perms=1;
        test_stat=dNBSglm(nbs.GLM);
        %Set back to original value
        nbs.GLM.perms=K;
    else
        test_stat=STATS.test_stat(1,:);
    end

    nbs.NBS.test_stat=zeros(N,N);
    nbs.NBS.test_stat(ind_uplo)=test_stat(:,:); 
    %nbs.NBS.test_stat=nbs.NBS.test_stat+nbs.NBS.test_stat';

    %Determine components satisfying alpha significance threshold
    n_cnt=0;
    nbs.NBS.con_mat = {};
    nbs.NBS.n=0;
    for i=1:length(sz_links)
        tmp=sum(null_dist>=sz_links(i))/K;
        if tmp<=STATS.alpha
            n_cnt=n_cnt+1;
            ind=find(adj==i);
            con_mat{n_cnt}=spalloc(N,N,length(ind)*2);
            con_mat{n_cnt}(ind)=1; 
            con_mat{n_cnt}=con_mat{n_cnt};
            pval(n_cnt)=tmp;
            nbs.NBS.pval = pval;
            nbs.NBS.con_mat = con_mat;
            nbs.NBS.n = n_cnt;
       end
    end
    if n_cnt==0
        pval=[]; con_mat=[]; 
    end
    %Display significant results with dNBSview only if node coordinates provided
    if nbs.NBS.n>0
       %  keyboard;
    end
    
    % NBS.node_coor=[];
    if nbs.NBS.n>0 && length(NBS.node_coor) > 1



        FolderName =  strcat(num2str(Threshold,'%2.2f'),'-Treshold-',GLM.test, '-Test-',num2str(GLM.perms),'-Perms-',num2str(pval),'-pValue') ;
        %currentFolder = pwd;
        mkdir( fullfile(currentFolder,'Result', strFolder,FolderName) );
        cd (fullfile(currentFolder,'Result', strFolder, FolderName))
        
        % nbs.NBS.ResultPath = fullfile(currentFolder,'Result', strFolder, FolderName);
        
        % let us change that:
        nbs.NBS.ResultPath = fullfile(currentFolder,'Result', strFolder, FolderName);
        % nbs.NBS.ResultPath = sprintf('%s/../../Results/%s', NBS.FolderConnMatFile, FolderName);
        
        
        nbs.NBS.node_label = NBS.node_label;
        nbs.NBS.node_coor=NBS.node_coor;
        
        %
        %
        %
        %
        % keyboard;
                dNBSview(nbs.NBS);
        %
        %
        %
        %
        save my_NBS.mat nbs
        
                
        %savefig(gcf,'NBSnetwork');NBS

       %saveas(gcf,'NBSnetwork','jpg'); %Save as jpeg image
        nbs.NBS.ResultPath = '';
        
       if STATS.thresh(1) ~= STATS.thresh(3)
       close(gcf);
       end
       adj=[];
       x2=[];
       % keyboard;
       
       
        for x=1:length(nbs.NBS.con_mat)
            
            
            
            adj=nbs.NBS.con_mat{x};
           %networkName = concenate( num2str(x,'%02d') +'.txt');
             networkName =  strcat('adjacency_matrix_network_binary_',num2str(x),'.txt') ;
            dlmwrite(networkName,full(adj),'delimiter',' ','precision','%d');
            
            x2=[];
            y2=[];
            z2=[];
            x2=adj;
             adj=nbs.NBS.test_stat;
           %networkName = concenate( num2str(x,'%02d') +'.txt');
            networkName =  strcat('Network_weighted_',num2str(x),'.txt') ;
            dlmwrite(networkName,full(adj),'delimiter',' ','precision','%d');
            y2=adj;
           
            %Weighted with zeros
            %networkName = concenate( num2str(x,'%02d') +'.txt');
            networkName =  strcat('adjacency_matrix_network_weighted_',num2str(x),'.txt') ;
            idxMatrix = find(abs(x2)<1);
            z2= y2;
            z2(idxMatrix)=0;
            dlmwrite(networkName,full(z2),'delimiter',' ','precision','%d');
            y2=adj;
            
            %end Weighted with zeros
            
            [i,j]=find(nbs.NBS.con_mat{x});
            StatkName =  strcat('test_stat_',num2str(x),'.txt' );
            
            fid = fopen(StatkName, 'wt');
           %Cohens d
            TStatTotal   = 0;
            Cohansd = 0;
            AmountSubjects = 0;
            %Cohens d
            for n=1:length(i)        
                i_lab=nbs.NBS.node_label{i(n)};
                j_lab=nbs.NBS.node_label{j(n)};
                stat=nbs.NBS.test_stat(i(n),j(n));
                %fprintf('%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat);

                fprintf(fid,'%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat);
                TStatTotal= TStatTotal + stat;
                %dlmwrite(StatkName,'%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat,'delimiter',' ','precision','%d');
            end
            AmountSubjects= size(C,3);
            TStatTotal = TStatTotal / length(i);
        %     Cohansd = TStatTotal/(sqrt(AmountSubjects));
            Cohansd = (2*TStatTotal) / (sqrt(AmountSubjects-2));
            fprintf(fid,'****************************');
            fprintf(fid,'\n');
            %fprintf(fid,...
            %        'Cohens d =  %6.2f',...
            %        Cohansd);
            %fprintf(fid,'\n');
            fclose(fid);
        end
        % cd(currentFolder);
     elseif nbs.NBS.n>0 && ~UI.node_coor.ok
         str='Significant result - specify Node Coordinates to view';
     
     else
         str=strcat('No significant result for Threshold: ', num2str(Threshold,'%2.2f'));
    end
    if ~isempty(str)
        try tmp=get(S.OUT.ls,'string'); set(S.OUT.ls,'string',[{str};tmp]); drawnow;
        catch;  end%fprintf([str,'\n']); end
    end
    if nbs.NBS.n>0 
        str=strcat('Significant result at Threshold:  ', num2str(Threshold,'%2.2f'),'  p-value= ', num2str(pval,'%0.4f')) ;
    end
    clearvars nbs;
    %Check if singel Threshold or threshold Range
    %If range clear pval and con mat for next calculation
    %if Singel keep the Results to show in dNBSview to interact
    if  STATS.thresh(1) ~= STATS.thresh(3)
        pval=[]; con_mat=[];
    end
    
    disp(str);

 end
 
end
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read node coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node_coor,ok]=read_node_coor(Name,DIMS)
ok=1;
data=readUI(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data);
    if nr==DIMS.nodes && nc==3 && ns==1 && isnumeric(data)
        node_coor=data; 
    else
        ok=0; node_coor=[];
    end
else
    ok=0; node_coor=[];
end        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read node labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node_label,ok]=read_node_label(Name,DIMS)
ok=1;
data=importdata(Name);
if ~isempty(data)
    [nr,nc,ns]=size(data);
    if nr==DIMS.nodes && nc==1 && ns==1
        node_label=data; 
    else
        ok=0; node_label=[];
    end
else
    ok=0; node_label=[]; 
end
end
